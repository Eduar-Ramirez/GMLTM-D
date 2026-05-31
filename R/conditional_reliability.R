#' Conditional Reliability based on Test Information Function (TIF)
#'
#' @description
#' Calculates conditional reliability using Test Information Function
#' for 3-parameter MLTM-D models. This approach is more precise than
#' quantile-based partitioning methods.
#'
#' @param fit Fitted model with \\eqn{\\alpha}, \\eqn{\\beta}, guessing parameters
#' @param theta_range Vector of \\eqn{\\theta} values where to evaluate reliability
#' @param component Integer or character. Specific component to analyze
#' @param n_samples Integer. Number of posterior samples to use
#'
#' @return List with conditional reliability by ability level
#'
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
  cat("Available posterior parameters:", paste(available_params, collapse = ", "), "\n")

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
  cat("Model dimensions:", dim(alpha_samples)[1], "iterations,",
      dim(alpha_samples)[2], "items,", dim(alpha_samples)[3], "components\n")

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
    cat(sprintf("Using %d out of %d posterior samples\n", n_samples, n_iter))
  }

  results <- list()

  for (comp in components_to_analyze) {
    cat(sprintf("Processing component %d/%d...\n", comp, n_comp))

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
#' @param filename Character. Filename if saving plot
#' @param ... Additional plotting parameters
#'
#' @export
plot_conditional_reliability <- function(results,
                                         component = NULL,
                                         plot_type = "both",
                                         include_ci = TRUE,
                                         ci_level = 0.95,
                                         color_scheme = "blue",
                                         add_reference_lines = TRUE,
                                         save_plot = FALSE,
                                         filename = "conditional_reliability.png",
                                         ...) {

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
  if (save_plot) {
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
              col = adjustcolor(colors$secondary, alpha = 0.3),
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
              col = adjustcolor(colors$info, alpha = 0.2),
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

  if (save_plot) {
    dev.off()
    cat("Plot saved as:", filename, "\n")
  }

  par(mfrow = c(1, 1))  # Reset graphical parameters
}

#' Plot Method for conditional_reliability_tif Objects
#'
#' @param x An object of class \code{conditional_reliability_tif}.
#' @param ... Additional arguments passed to \code{plot_conditional_reliability}.
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
              col = adjustcolor(colors[i], alpha = 0.2),
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

#' Compare Conditional Reliability Between Components at Specific Theta Values
#'
#' @param cond_rel_obj An object of class \code{conditional_reliability_tif} returned by
#'   \code{\link{conditional_reliability_tif}}.
#' @param theta_points Numeric vector of theta values at which to compare components.
#'   Default is \code{c(-1, 0, 1)}.
#' @export
compare_conditional_reliability <- function(cond_rel_obj, theta_points = c(-1, 0, 1)) {

  if (!inherits(cond_rel_obj, "conditional_reliability_tif")) {
    stop("Object must be of class 'conditional_reliability_tif'")
  }

  n_comp <- length(cond_rel_obj)
  comp_names <- names(cond_rel_obj)

  cat("CONDITIONAL RELIABILITY COMPARISON AT SPECIFIC theta VALUES\n")
  cat("=======================================================\n\n")

  for (theta_point in theta_points) {
    # Find closest theta to requested point
    theta_vals <- cond_rel_obj[[1]]$theta_values
    closest_idx <- which.min(abs(theta_vals - theta_point))
    actual_theta <- theta_vals[closest_idx]

    # Extract reliabilities at this point
    reliabilities <- sapply(cond_rel_obj, function(x) x$reliability_mean[closest_idx])
    names(reliabilities) <- comp_names

    # Create ranking
    ranking <- sort(reliabilities, decreasing = TRUE)

    cat(sprintf("At theta = %.2f:\n", actual_theta))
    cat(paste(rep("-", 15), collapse = ""), "\n")

    for (i in seq_along(ranking)) {
      rel_val <- ranking[i]
      comp_name <- names(ranking)[i]

      status <- if (rel_val >= 0.9) "EXCELLENT"
      else if (rel_val >= 0.8) "GOOD"
      else if (rel_val >= 0.7) "ACCEPTABLE"
      else if (rel_val >= 0.6) "QUESTIONABLE"
      else "POOR"

      cat(sprintf("%d. %s: %.3f (%s)\n", i, comp_name, rel_val, status))
    }
    cat("\n")
  }

  return(invisible(NULL))
}

#' Quick Reliability Profile Analysis
#'
#' @description
#' Provides a quick overview of reliability characteristics for each component.
#'
#' @param cond_rel_obj An object of class \code{conditional_reliability_tif} returned by
#'   \code{\link{conditional_reliability_tif}}.
#' @return A list of class \code{"reliability_profile"} with summary metrics per component.
#' @export
reliability_profile <- function(cond_rel_obj) {

  if (!inherits(cond_rel_obj, "conditional_reliability_tif")) {
    stop("Object must be of class 'conditional_reliability_tif'")
  }

  profiles <- list()

  for (i in seq_along(cond_rel_obj)) {
    comp_name <- names(cond_rel_obj)[i]
    comp_data <- cond_rel_obj[[i]]

    rel_mean <- comp_data$reliability_mean
    theta_vals <- comp_data$theta_values

    # Calculate profile metrics
    profile <- list(
      component = comp_name,
      max_reliability = max(rel_mean, na.rm = TRUE),
      min_reliability = min(rel_mean, na.rm = TRUE),
      optimal_theta = comp_data$optimal_theta,
      reliability_range = max(rel_mean, na.rm = TRUE) - min(rel_mean, na.rm = TRUE),
      theta_range_good = {
        good_indices <- which(rel_mean >= 0.8)
        if (length(good_indices) > 0) {
          c(min(theta_vals[good_indices]), max(theta_vals[good_indices]))
        } else {
          c(NA, NA)
        }
      },
      mean_reliability = mean(rel_mean, na.rm = TRUE),
      reliability_stability = 1 - (sd(rel_mean, na.rm = TRUE) / mean(rel_mean, na.rm = TRUE))
    )

    profiles[[comp_name]] <- profile
  }

  class(profiles) <- "reliability_profile"
  return(profiles)
}

#' Print Method for Reliability Profile
#'
#' @param x An object of class \code{"reliability_profile"}.
#' @param ... Currently unused.
#' @export
print.reliability_profile <- function(x, ...) {

  cat("RELIABILITY PROFILE SUMMARY\n")
  cat("===========================\n\n")

  for (i in seq_along(x)) {
    profile <- x[[i]]

    cat(sprintf("Component: %s\n", profile$component))
    cat(sprintf("  Max Reliability: %.3f at theta = %.2f\n",
                profile$max_reliability, profile$optimal_theta))
    cat(sprintf("  Reliability Range: %.3f - %.3f (Delta = %.3f)\n",
                profile$min_reliability, profile$max_reliability, profile$reliability_range))
    cat(sprintf("  Mean Reliability: %.3f\n", profile$mean_reliability))
    cat(sprintf("  Stability Index: %.3f\n", profile$reliability_stability))

    if (!is.na(profile$theta_range_good[1])) {
      cat(sprintf("  Good Reliability theta Range: [%.2f, %.2f]\n",
                  profile$theta_range_good[1], profile$theta_range_good[2]))
    } else {
      cat("  Good Reliability theta Range: None (no theta with reliability >= 0.8)\n")
    }
    cat("\n")
  }

  invisible(x)
}

#' Integration with Enhanced Reliability Analysis
#'
#' @description
#' Integrates conditional TIF analysis with existing reliability functions.
#'
#' @param fit A fitted GMLTM, MLTM, or LLTM model object.
#' @param include_conditional_tif Logical. Whether to compute conditional TIF-based
#'   reliability. Default is \code{TRUE}.
#' @param theta_range Numeric vector of theta values for conditional reliability evaluation.
#'   Default is \code{seq(-3, 3, 0.1)}.
#' @param ... Additional arguments passed to \code{\link{enhanced_mltm_reliability}}.
#' @return A list combining enhanced reliability results and, optionally, conditional
#'   TIF-based reliability.
#' @export
integrate_with_enhanced_reliability <- function(fit,
                                                include_conditional_tif = TRUE,
                                                theta_range = seq(-3, 3, 0.1),
                                                ...) {

  # Standard analysis (assumes this function exists)
  if (exists("enhanced_mltm_reliability")) {
    standard_results <- enhanced_mltm_reliability(fit, ...)
  } else {
    standard_results <- list()
    cat("Note: enhanced_mltm_reliability function not found. Creating basic structure.\n")
  }

  # Add conditional TIF analysis if requested
  if (include_conditional_tif) {
    cat("Calculating conditional reliability based on TIF...\n")

    conditional_tif <- conditional_reliability_tif(fit, theta_range = theta_range)
    standard_results$conditional_tif <- conditional_tif
  }

  return(standard_results)
}

#' Plot All Components Separately
#'
#' @description
#' Creates individual plots for each component in the analysis.
#'
#' @param results An object of class \code{conditional_reliability_tif}.
#' @param plot_type Character. Type of plot: \code{"reliability"}, \code{"information"},
#'   or \code{"both"}. Default is \code{"both"}.
#' @param include_ci Logical. Whether to include credible interval bands. Default \code{TRUE}.
#' @param color_scheme Character. Color scheme for plots. Default \code{"blue"}.
#' @param save_plots Logical. Whether to save plots to disk. Default \code{FALSE}.
#' @param output_dir Character. Directory for saved plots. Default \code{"plots"}.
#' @param ... Additional arguments passed to \code{\link{plot_conditional_reliability}}.
#' @export
plot_all_components <- function(results,
                                plot_type = "both",
                                include_ci = TRUE,
                                color_scheme = "blue",
                                save_plots = FALSE,
                                output_dir = "plots",
                                ...) {

  if (!inherits(results, "conditional_reliability_tif")) {
    stop("Object must be of class 'conditional_reliability_tif'")
  }

  n_components <- length(results)
  comp_names <- names(results)

  # Create output directory if saving plots
  if (save_plots && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Plot each component
  for (i in 1:n_components) {
    comp_name <- comp_names[i]

    if (save_plots) {
      filename <- file.path(output_dir, paste0("conditional_reliability_",
                                               gsub("[^A-Za-z0-9]", "_", comp_name), ".png"))
    } else {
      filename <- NULL
    }

    cat(sprintf("Creating plot for component: %s\n", comp_name))

    plot_conditional_reliability(results,
                                 component = i,
                                 plot_type = plot_type,
                                 include_ci = include_ci,
                                 color_scheme = color_scheme,
                                 save_plot = save_plots,
                                 filename = filename,
                                 ...)

    # Pause between plots if not saving
    if (!save_plots && i < n_components) {
      readline(prompt = "Press [Enter] to continue to next component...")
    }
  }

  if (save_plots) {
    cat(sprintf("All plots saved to: %s\n", output_dir))
  }
}

#' Create Comparison Plot of All Components
#'
#' @description
#' Creates a single plot comparing conditional reliability across all components.
#'
#' @param results An object of class \code{conditional_reliability_tif}.
#' @param include_ci Logical. Whether to include credible interval bands. Default \code{TRUE}.
#' @param show_optimal_points Logical. Whether to mark optimal theta points. Default \code{TRUE}.
#' @param add_reference_lines Logical. Whether to add horizontal reference lines at 0.7, 0.8, 0.9.
#'   Default \code{TRUE}.
#' @param color_palette Character. RColorBrewer palette name for component colors.
#'   Default \code{"Set2"}.
#' @param save_plot Logical. Whether to save the plot. Default \code{FALSE}.
#' @param filename Character. Output filename if \code{save_plot = TRUE}.
#'   Default \code{"components_comparison.png"}.
#' @param ... Additional graphical arguments.
#' @importFrom RColorBrewer brewer.pal
#' @export
plot_components_comparison <- function(results,
                                       include_ci = TRUE,
                                       show_optimal_points = TRUE,
                                       add_reference_lines = TRUE,
                                       color_palette = "Set2",
                                       save_plot = FALSE,
                                       filename = "components_comparison.png",
                                       ...) {

  if (!inherits(results, "conditional_reliability_tif")) {
    stop("Object must be of class 'conditional_reliability_tif'")
  }

  n_comp <- length(results)
  comp_names <- names(results)
  theta_vals <- results[[1]]$theta_values

  # Color palette
  if (requireNamespace("RColorBrewer", quietly = TRUE)) {
    if (n_comp <= 8) {
      colors <- RColorBrewer::brewer.pal(max(3, n_comp), color_palette)[1:n_comp]
    } else {
      colors <- rainbow(n_comp)
    }
  } else {
    colors <- rainbow(n_comp)
  }

  # Setup plot
  if (save_plot) {
    png(filename, width = 12, height = 8, units = "in", res = 300)
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
    comp_data <- results[[i]]
    rel_mean <- comp_data$reliability_mean

    # Confidence band
    if (include_ci) {
      rel_lower <- comp_data$reliability_q025
      rel_upper <- comp_data$reliability_q975

      polygon(c(theta_vals, rev(theta_vals)),
              c(rel_lower, rev(rel_upper)),
              col = adjustcolor(colors[i], alpha = 0.2),
              border = NA)
    }

    # Main line
    lines(theta_vals, rel_mean, lwd = 3, col = colors[i])

    # Optimal theta point
    if (show_optimal_points) {
      points(comp_data$optimal_theta, comp_data$max_reliability,
             pch = 19, col = colors[i], cex = 1.2)
    }
  }

  # Reference lines
  if (add_reference_lines) {
    ref_levels <- c(0.6, 0.7, 0.8, 0.9)
    ref_labels <- c("Poor", "Questionable", "Acceptable", "Good")
    ref_colors <- c("red", "orange", "gold", "darkgreen")

    for (i in seq_along(ref_levels)) {
      abline(h = ref_levels[i], lty = 2, col = ref_colors[i], alpha = 0.7)
    }
  }

  # Legend
  legend("topright",
         legend = comp_names,
         col = colors, lwd = 3, bty = "n", cex = 0.9)

  grid(col = "gray90", lty = 3)

  if (save_plot) {
    dev.off()
    cat("Comparison plot saved as:", filename, "\n")
  }
}

# Helper function for null coalescing
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}
