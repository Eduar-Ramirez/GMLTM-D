#' Conditional Reliability based on Test Information Function (TIF)
#'
#' @description
#' Calculates conditional reliability using Test Information Function
#' for 3-parameter MLTM-D models. This approach is more precise than
#' quantile-based partitioning methods.
#'
#' @param fit Fitted model with \eqn{\alpha}, \eqn{\beta}, guessing parameters
#' @param theta_range Vector of \eqn{\theta} values where to evaluate reliability
#' @param component Integer or character. Specific component to analyze
#' @param n_samples Integer. Number of posterior samples to use
#'
#' @return A list of class \code{"conditional_reliability_tif"} with elements:
#'   \describe{
#'     \item{\code{theta}}{Numeric vector of theta values.}
#'     \item{\code{reliability}}{Numeric vector of reliability estimates
#'       at each theta value.}
#'     \item{\code{information}}{Numeric vector of test information
#'       values at each theta value.}
#'     \item{\code{component}}{Integer indicating the model component.}
#'     \item{\code{fit}}{The original fitted model object.}
#'   }
#'
#' @examples
#' \donttest{
#'   if (!requireNamespace("rstan", quietly = TRUE)) return()
#'   data(analogy)
#'   Q <- structure(
#'     c(0,0,1,0,1,0,1,0,1,1,0,1,1,1,0,1,1,1,0,1,0,1,0,0,1,0,1,
#'       1,0,0,0,0,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,0,
#'       1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,1,0,1,0,1,1,0,1,0,0,0,
#'       0,0,0,0,0,0,1,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,0,
#'       1,0,0,0,1,1,0,1,1,1,1,1,1,1,0,1,1,0,1,1,1,1,0,1),
#'     dim = c(27L, 5L),
#'     dimnames = list(NULL, c("rot_fig","rot_trap","reflection",
#'                             "subt_seg","mov_point")))
#'   components <- list(global = c(1, 2, 3), local = c(4, 5))
#'   fit <- GMLTM(analogy, Q, components,
#'                iters = 200, iter_warmup = 100, chains = 1, cores = 1)
#'   cond_rel <- conditional_reliability_tif(fit, theta_range = seq(-2, 2, 0.5))
#'   cond_rel$global$optimal_theta
#' }
#'
#' @family conditional reliability functions
#' @export
conditional_reliability_tif <- function(fit,
                                        theta_range = seq(-3, 3, 0.2),
                                        component = NULL,
                                        n_samples = 1000) {

  # Input validation
  if (!inherits(fit, c("GMLTM", "LLTM", "MLTM"))) {
    stop("Incompatible model type")
  }

  # Show available parameters for debugging
  available_params <- names(fit$posterior)
  message("Available posterior parameters: ", paste(available_params, collapse = ", "))

  # Check for gamma/guessing parameters
  gamma_param_name <- NULL
  if ("gamma" %in% available_params) {
    gamma_param_name <- "gamma"
  } else if ("guessing" %in% available_params) {
    gamma_param_name <- "guessing"
  } else if ("c" %in% available_params) {
    gamma_param_name <- "c"
  } else {
    stop("Guessing parameter not found (gamma/guessing/c). Available: ",
         paste(available_params, collapse = ", "))
  }

  required_params <- c("alpha", "beta", gamma_param_name)
  if (!all(required_params %in% available_params)) {
    stop("Missing model parameters (alpha, beta, gamma). Available: ",
         paste(available_params, collapse = ", "))
  }

  # Extract posterior parameters
  alpha_samples <- fit$posterior$alpha  # [iter, item, component]
  beta_samples <- fit$posterior$beta    # [iter, item, component]
  gamma_samples <- fit$posterior[[gamma_param_name]]  # [iter, item] or [iter, item, component]

  # Show dimensions for debugging
  message("Model dimensions: ", dim(alpha_samples)[1], " iterations, ",
      dim(alpha_samples)[2], " items, ", dim(alpha_samples)[3], " components")

  # Handle data structure for LLTM models
  if (inherits(fit, "LLTM")) {
    # Convert to 3D structure if necessary
    alpha_samples <- array(alpha_samples, dim = c(dim(alpha_samples), 1))
    beta_samples <- array(beta_samples, dim = c(dim(beta_samples), 1))
    if (length(dim(gamma_samples)) == 2) {
      gamma_samples <- array(gamma_samples, dim = c(dim(gamma_samples), 1))
    }
  }

  n_iter <- dim(alpha_samples)[1]
  n_items <- dim(alpha_samples)[2]
  n_comp <- dim(alpha_samples)[3]

  # Handle gamma_samples if it's 2D (shared across components)
  if (length(dim(gamma_samples)) == 2) {
    # gamma is shared across components - expand to 3D
    gamma_samples_3d <- array(NA, dim = c(n_iter, n_items, n_comp))
    for (comp in 1:n_comp) {
      gamma_samples_3d[, , comp] <- gamma_samples
    }
    gamma_samples <- gamma_samples_3d
  }

  # Select component(s) to analyze
  if (is.null(component)) {
    components_to_analyze <- 1:n_comp
  } else if (is.numeric(component)) {
    components_to_analyze <- component
  } else {
    # If character, search by name
    comp_names <- colnames(fit$EAP$theta) %||% paste0("Component_", 1:n_comp)
    components_to_analyze <- which(comp_names == component)
    if (length(components_to_analyze) == 0) {
      stop("Component '", component, "' not found. Available: ", paste(comp_names, collapse = ", "))
    }
  }

  # Limit samples for speed
  if (n_samples < n_iter) {
    sample_idx <- sample(1:n_iter, n_samples)
    alpha_samples <- alpha_samples[sample_idx, , , drop = FALSE]
    beta_samples <- beta_samples[sample_idx, , , drop = FALSE]
    gamma_samples <- gamma_samples[sample_idx, , , drop = FALSE]
    message(sprintf("Using %d out of %d posterior samples", n_samples, n_iter))
  }

  results <- list()

  for (comp in components_to_analyze) {
    message(sprintf("Processing component %d/%d...", comp, n_comp))

    # Information matrix for each theta and posterior sample
    information_matrix <- array(NA, dim = c(length(theta_range), nrow(alpha_samples)))

    for (s in 1:nrow(alpha_samples)) {
      # Parameters for this posterior sample
      alpha_s <- alpha_samples[s, , comp]
      beta_s <- beta_samples[s, , comp]
      gamma_s <- gamma_samples[s, , comp]

      for (t in seq_along(theta_range)) {
        theta <- theta_range[t]

        # Information for each item at this theta
        item_info <- numeric(n_items)

        for (i in 1:n_items) {
          # 3-parameter model: P(theta) = gamma + (1-gamma) * exp(alpha*(theta-beta)) / (1 + exp(alpha*(theta-beta)))
          z <- alpha_s[i] * (theta - beta_s[i])

          # Avoid numerical overflow
          if (z > 700) z <- 700
          if (z < -700) z <- -700

          exp_z <- exp(z)

          # Probability of correct response
          P <- gamma_s[i] + (1 - gamma_s[i]) * exp_z / (1 + exp_z)

          # Derivative of P with respect to theta
          dP_dtheta <- (1 - gamma_s[i]) * alpha_s[i] * exp_z / ((1 + exp_z)^2)

          # Item information: I(theta) = [dP/dtheta]^2 / [P(1-P)]
          if (P > 0.001 && P < 0.999) {  # Avoid division by zero
            item_info[i] <- (dP_dtheta^2) / (P * (1 - P))
          } else {
            item_info[i] <- 0
          }
        }

        # Total test information at this theta
        information_matrix[t, s] <- sum(item_info, na.rm = TRUE)
      }
    }

    # Conditional reliability: R(theta) = I(theta) / (1 + I(theta))
    # where I(theta) is test information
    reliability_matrix <- information_matrix / (1 + information_matrix)

    # Summary statistics
    comp_name <- if (!is.null(colnames(fit$EAP$theta))) {
      colnames(fit$EAP$theta)[comp]
    } else {
      paste0("Component_", comp)
    }

    results[[comp_name]] <- list(
      theta_values = theta_range,
      reliability_mean = apply(reliability_matrix, 1, mean, na.rm = TRUE),
      reliability_sd = apply(reliability_matrix, 1, sd, na.rm = TRUE),
      reliability_q025 = apply(reliability_matrix, 1, quantile, 0.025, na.rm = TRUE),
      reliability_q975 = apply(reliability_matrix, 1, quantile, 0.975, na.rm = TRUE),
      reliability_q05 = apply(reliability_matrix, 1, quantile, 0.05, na.rm = TRUE),
      reliability_q95 = apply(reliability_matrix, 1, quantile, 0.95, na.rm = TRUE),
      information_mean = apply(information_matrix, 1, mean, na.rm = TRUE),
      information_sd = apply(information_matrix, 1, sd, na.rm = TRUE),
      information_q025 = apply(information_matrix, 1, quantile, 0.025, na.rm = TRUE),
      information_q975 = apply(information_matrix, 1, quantile, 0.975, na.rm = TRUE),
      optimal_theta = theta_range[which.max(apply(reliability_matrix, 1, mean, na.rm = TRUE))],
      max_reliability = max(apply(reliability_matrix, 1, mean, na.rm = TRUE), na.rm = TRUE),
      samples = reliability_matrix  # For additional analyses
    )
  }

  class(results) <- "conditional_reliability_tif"
  attr(results, "theta_range") <- theta_range
  attr(results, "model_type") <- class(fit)[1]

  return(results)
}

#' Plot Conditional Reliability Results
#'
#' @description
#' Creates plots for conditional reliability analysis with multiple visualization options
#'
#' @param results Object of class conditional_reliability_tif
#' @param component Integer or character. Component to plot (NULL for first component)
#' @param plot_type Character. Type of plot: "reliability", "information", "both", "comparison"
#' @param include_ci Logical. Include confidence intervals
#' @param ci_level Numeric. Confidence level (0.90 or 0.95)
#' @param color_scheme Character. Color scheme: "blue", "viridis", "custom"
#' @param add_reference_lines Logical. Add reference lines for reliability levels
#' @param save_plot Logical. Save plot to file
#' @param filename Character. Filename if saving plot. Default \code{NULL}
#'   (no file written). Only writes when both \code{save_plot = TRUE} and
#'   \code{filename} is non-\code{NULL}.
#' @param ... Additional plotting parameters
#'
#' @return Invisibly returns \code{NULL}. Called for its side effect of
#'   producing one or two plots (reliability curve, test information
#'   function, or both) depending on \code{plot_type}.
#'
#' @examples
#' \donttest{
#'   if (!requireNamespace("rstan", quietly = TRUE)) return()
#'   data(analogy)
#'   Q <- structure(
#'     c(0,0,1,0,1,0,1,0,1,1,0,1,1,1,0,1,1,1,0,1,0,1,0,0,1,0,1,
#'       1,0,0,0,0,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,0,
#'       1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,1,0,1,0,1,1,0,1,0,0,0,
#'       0,0,0,0,0,0,1,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,0,
#'       1,0,0,0,1,1,0,1,1,1,1,1,1,1,0,1,1,0,1,1,1,1,0,1),
#'     dim = c(27L, 5L),
#'     dimnames = list(NULL, c("rot_fig","rot_trap","reflection",
#'                             "subt_seg","mov_point")))
#'   components <- list(global = c(1, 2, 3), local = c(4, 5))
#'   fit <- GMLTM(analogy, Q, components,
#'                iters = 200, iter_warmup = 100, chains = 1, cores = 1)
#'   cond_rel <- conditional_reliability_tif(fit, theta_range = seq(-2, 2, 0.5))
#'   plot_conditional_reliability(cond_rel, component = 1)
#' }
#'
#' @family conditional reliability functions
#' @export
plot_conditional_reliability <- function(results,
                                         component = NULL,
                                         plot_type = "both",
                                         include_ci = TRUE,
                                         ci_level = 0.95,
                                         color_scheme = "blue",
                                         add_reference_lines = TRUE,
                                         save_plot = FALSE,
                                         filename = NULL,
                                         ...) {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)

  # Check if it's the right class
  if (!inherits(results, "conditional_reliability_tif")) {
    stop("Object must be of class 'conditional_reliability_tif'")
  }

  # Select component
  if (is.null(component)) {
    comp_data <- results[[1]]
    comp_name <- names(results)[1]
  } else if (is.numeric(component)) {
    if (component > length(results)) {
      stop("Component ", component, " not found. Available components: 1 to ", length(results))
    }
    comp_data <- results[[component]]
    comp_name <- names(results)[component]
  } else {
    if (!component %in% names(results)) {
      stop("Component '", component, "' not found. Available: ", paste(names(results), collapse = ", "))
    }
    comp_data <- results[[component]]
    comp_name <- component
  }

  # Extract data
  theta_vals <- comp_data$theta_values
  rel_mean <- comp_data$reliability_mean
  info_mean <- comp_data$information_mean

  # Confidence intervals
  if (ci_level == 0.95) {
    rel_lower <- comp_data$reliability_q025
    rel_upper <- comp_data$reliability_q975
    info_lower <- comp_data$information_q025
    info_upper <- comp_data$information_q975
  } else {
    rel_lower <- comp_data$reliability_q05 %||% (rel_mean - comp_data$reliability_sd)
    rel_upper <- comp_data$reliability_q95 %||% (rel_mean + comp_data$reliability_sd)
    info_lower <- (info_mean - comp_data$information_sd)
    info_upper <- (info_mean + comp_data$information_sd)
  }

  # Color schemes
  colors <- switch(color_scheme,
                   "blue" = list(main = "darkblue", secondary = "lightblue", info = "darkgreen"),
                   "viridis" = list(main = "#440154", secondary = "#21908C", info = "#FDE725"),
                   "custom" = list(main = "#2E86AB", secondary = "#A23B72", info = "#F18F01"),
                   list(main = "darkblue", secondary = "lightblue", info = "darkgreen")  # default
  )

  # Setup plot
  if (save_plot && !is.null(filename)) {
    png(filename, width = 12, height = 8, units = "in", res = 300)
  }

  # Configure plot layout
  if (plot_type == "both") {
    par(mfrow = c(2, 1), mar = c(4, 4, 3, 2), oma = c(0, 0, 2, 0))
  } else if (plot_type == "comparison" && length(results) > 1) {
    par(mfrow = c(1, 1), mar = c(4, 4, 3, 2))
  } else {
    par(mfrow = c(1, 1), mar = c(4, 4, 3, 2))
  }

  # Plot 1: Conditional Reliability
  if (plot_type %in% c("reliability", "both")) {
    plot(theta_vals, rel_mean,
         type = "l", lwd = 3, col = colors$main,
         ylim = c(0, 1),
         xlab = expression(theta ~ "(Ability Level)"),
         ylab = "Conditional Reliability",
         main = paste("Conditional Reliability -", comp_name),
         ...)

    # Confidence band
    if (include_ci && !is.null(rel_lower) && !is.null(rel_upper)) {
      polygon(c(theta_vals, rev(theta_vals)),
              c(rel_lower, rev(rel_upper)),
              col = adjustcolor(colors$secondary, alpha.f = 0.3),
              border = NA)

      # Redraw main line
      lines(theta_vals, rel_mean, lwd = 3, col = colors$main)
    }

    # Reference lines
    if (add_reference_lines) {
      ref_levels <- c(0.6, 0.7, 0.8, 0.9)
      ref_labels <- c("Poor", "Questionable", "Acceptable", "Good")
      ref_colors <- c("red", "orange", "gold", "darkgreen")

      for (i in seq_along(ref_levels)) {
        abline(h = ref_levels[i], lty = 2, col = ref_colors[i], lwd = 1)
      }

      # Legend for reference lines
      legend("topright",
             legend = paste(ref_labels, "(>=", ref_levels, ")"),
             col = ref_colors, lty = 2, cex = 0.7, bty = "n")
    }

    # Optimal theta line
    abline(v = comp_data$optimal_theta, lty = 2, col = "red", lwd = 2)
    text(comp_data$optimal_theta, 0.95,
         sprintf("Optimal theta = %.2f\nMax Rel = %.3f",
                 comp_data$optimal_theta, comp_data$max_reliability),
         adj = c(0.5, 1), col = "red", cex = 0.8, font = 2)

    grid(col = "gray90", lty = 3)
  }

  # Plot 2: Test Information Function
  if (plot_type %in% c("information", "both")) {
    plot(theta_vals, info_mean,
         type = "l", lwd = 3, col = colors$info,
         xlab = expression(theta ~ "(Ability Level)"),
         ylab = "Test Information",
         main = paste("Test Information Function -", comp_name),
         ...)

    # Confidence band for information
    if (include_ci && !is.null(info_lower) && !is.null(info_upper)) {
      polygon(c(theta_vals, rev(theta_vals)),
              c(info_lower, rev(info_upper)),
              col = adjustcolor(colors$info, alpha.f = 0.2),
              border = NA)

      # Redraw main line
      lines(theta_vals, info_mean, lwd = 3, col = colors$info)
    }

    # Optimal theta line
    abline(v = comp_data$optimal_theta, lty = 2, col = "red", lwd = 2)

    # Maximum information point
    max_info_idx <- which.max(info_mean)
    max_info_theta <- theta_vals[max_info_idx]
    max_info_value <- info_mean[max_info_idx]

    points(max_info_theta, max_info_value, pch = 19, col = "red", cex = 1.5)
    text(max_info_theta, max_info_value,
         sprintf("Max Info = %.2f\nat theta = %.2f", max_info_value, max_info_theta),
         adj = c(0.5, -0.2), col = "red", cex = 0.8, font = 2)

    grid(col = "gray90", lty = 3)
  }

  # Plot 3: Comparison between components
  if (plot_type == "comparison" && length(results) > 1) {
    plot_comparison_conditional(results, include_ci = include_ci,
                                color_scheme = color_scheme, ...)
  }

  # Overall title
  if (plot_type == "both") {
    mtext(paste("Conditional Reliability Analysis -", comp_name),
          outer = TRUE, cex = 1.2, font = 2)
  }

  if (save_plot && !is.null(filename)) {
    dev.off()
    message("Plot saved as: ", filename)
  }

  invisible(NULL)
}

#' Plot Method for conditional_reliability_tif Objects
#'
#' @param x An object of class \code{conditional_reliability_tif}.
#' @param ... Additional arguments passed to \code{plot_conditional_reliability}.
#' @return Invisibly returns \code{NULL}. Called for its side effect of
#'   producing a reliability and/or information plot via
#'   \code{plot_conditional_reliability()}.
#' @export
plot.conditional_reliability_tif <- function(x, ...) {
  plot_conditional_reliability(x, ...)
}

#' Plot Comparison of Conditional Reliability Between Components
#'
#' @keywords internal
plot_comparison_conditional <- function(x, include_ci = TRUE,
                                        color_scheme = "blue", ...) {

  n_comp <- length(x)
  comp_names <- names(x)
  theta_vals <- x[[1]]$theta_values

  # Color palette
  if (n_comp <= 8) {
    colors <- rainbow(n_comp)
  } else {
    colors <- rainbow(n_comp)
  }

  # Initialize plot
  plot(range(theta_vals), c(0, 1),
       type = "n",
       xlab = expression(theta ~ "(Ability Level)"),
       ylab = "Conditional Reliability",
       main = "Conditional Reliability Comparison - All Components",
       ...)

  # Plot each component
  for (i in 1:n_comp) {
    comp_data <- x[[i]]
    rel_mean <- comp_data$reliability_mean

    # Confidence band
    if (include_ci) {
      rel_lower <- comp_data$reliability_q025
      rel_upper <- comp_data$reliability_q975

      polygon(c(theta_vals, rev(theta_vals)),
              c(rel_lower, rev(rel_upper)),
              col = adjustcolor(colors[i], alpha.f = 0.2),
              border = NA)
    }

    # Main line
    lines(theta_vals, rel_mean, lwd = 3, col = colors[i])

    # Optimal theta point
    points(comp_data$optimal_theta, comp_data$max_reliability,
           pch = 19, col = colors[i], cex = 1.2)
  }

  # Reference lines
  abline(h = c(0.6, 0.7, 0.8, 0.9), lty = 2, col = "gray60")

  # Legend
  legend("topright",
         legend = comp_names,
         col = colors, lwd = 3, bty = "n", cex = 0.9)

  grid(col = "gray90", lty = 3)
}

#' Summary Method for Conditional Reliability Analysis
#'
#' @param object An object of class \code{conditional_reliability_tif}.
#' @param digits Integer. Number of decimal places to display. Default is \code{3}.
#' @param ... Currently unused.
#' @return Invisibly returns \code{object}. Called for its side effect of
#'   printing a formatted summary of conditional reliability statistics to the
#'   console.
#' @export
summary.conditional_reliability_tif <- function(object, digits = 3, ...) {

  cat("CONDITIONAL RELIABILITY ANALYSIS (Test Information Function)\n")
  cat("============================================================\n\n")

  for (i in seq_along(object)) {
    comp_name <- names(object)[i]
    comp_data <- object[[i]]

    cat(sprintf("Component: %s\n", comp_name))
    cat(paste(rep("-", nchar(comp_name) + 11), collapse = ""), "\n")

    cat(sprintf("Optimal theta: %.*f (Max Reliability: %.*f)\n",
                digits, comp_data$optimal_theta, digits, comp_data$max_reliability))

    # Reliability ranges
    excellent_rel <- comp_data$theta_values[comp_data$reliability_mean >= 0.9]
    good_rel <- comp_data$theta_values[comp_data$reliability_mean >= 0.8]
    acceptable_rel <- comp_data$theta_values[comp_data$reliability_mean >= 0.7]

    if (length(excellent_rel) > 0) {
      cat(sprintf("Excellent reliability range (>=0.9): theta in [%.*f, %.*f]\n",
                  digits, min(excellent_rel), digits, max(excellent_rel)))
    }

    if (length(good_rel) > 0) {
      cat(sprintf("Good reliability range (>=0.8): theta in [%.*f, %.*f]\n",
                  digits, min(good_rel), digits, max(good_rel)))
    } else {
      cat("No range with good reliability (>=0.8)\n")
    }

    if (length(acceptable_rel) > 0) {
      cat(sprintf("Acceptable reliability range (>=0.7): theta in [%.*f, %.*f]\n",
                  digits, min(acceptable_rel), digits, max(acceptable_rel)))
    }

    # Information statistics
    cat(sprintf("Average test information: %.*f (SD: %.*f)\n",
                digits, mean(comp_data$information_mean), digits, mean(comp_data$information_sd)))

    # Reliability at key theta values
    key_thetas <- c(-2, -1, 0, 1, 2)
    cat("Reliability at key theta values:\n")
    for (theta_val in key_thetas) {
      idx <- which.min(abs(comp_data$theta_values - theta_val))
      if (length(idx) > 0) {
        rel_val <- comp_data$reliability_mean[idx]
        actual_theta <- comp_data$theta_values[idx]
        cat(sprintf("  theta = %.*f: Reliability = %.*f\n",
                    digits, actual_theta, digits, rel_val))
      }
    }

    cat("\n")
  }

  invisible(object)
}
