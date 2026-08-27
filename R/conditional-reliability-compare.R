#' Compare Conditional Reliability Between Components at Specific Theta Values
#'
#' @param cond_rel_obj An object of class \code{conditional_reliability_tif} returned by
#'   \code{\link{conditional_reliability_tif}}.
#' @param theta_points Numeric vector of theta values at which to compare components.
#'   Default is \code{c(-1, 0, 1)}.
#' @return Invisibly returns a data frame with columns \code{theta},
#'   \code{reliability}, and \code{se} evaluated at the requested
#'   \code{theta_points}. Called primarily for its side effect of
#'   printing a formatted summary table to the console.
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
#'   compare_conditional_reliability(cond_rel, theta_points = c(-1, 0, 1))
#' }
#'
#' @family conditional reliability functions
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
#' @return A list of class \code{"reliability_profile"} with elements:
#'   \describe{
#'     \item{\code{summary}}{Data frame with columns \code{theta},
#'       \code{reliability}, \code{lower}, and \code{upper}.}
#'     \item{\code{component}}{Integer indicating which model component
#'       was analysed.}
#'   }
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
#'   reliability_profile(cond_rel)
#' }
#'
#' @family conditional reliability functions
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
#' @return Invisibly returns \code{x}. Called for its side effect of printing
#'   a formatted reliability profile summary to the console.
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
#'   result <- integrate_with_enhanced_reliability(fit,
#'                theta_range = seq(-2, 2, 0.5))
#'   result$conditional_tif
#' }
#'
#' @family conditional reliability functions
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
    message("Note: enhanced_mltm_reliability function not found. Creating basic structure.")
  }

  # Add conditional TIF analysis if requested
  if (include_conditional_tif) {
    message("Calculating conditional reliability based on TIF...")

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
#' @param output_dir Character. Directory for saved plots. Default \code{NULL}
#'   (no files written). Only creates directory and saves when both
#'   \code{save_plots = TRUE} and \code{output_dir} is non-\code{NULL}.
#' @param ... Additional arguments passed to \code{\link{plot_conditional_reliability}}.
#' @return Invisibly returns \code{NULL}. Called for its side effect of
#'   producing reliability plots for all model components.
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
#'   plot_all_components(cond_rel, save_plots = TRUE, output_dir = tempdir())
#' }
#'
#' @family conditional reliability functions
#' @export
plot_all_components <- function(results,
                                plot_type = "both",
                                include_ci = TRUE,
                                color_scheme = "blue",
                                save_plots = FALSE,
                                output_dir = NULL,
                                ...) {

  if (!inherits(results, "conditional_reliability_tif")) {
    stop("Object must be of class 'conditional_reliability_tif'")
  }

  n_components <- length(results)
  comp_names <- names(results)

  # Create output directory if saving plots
  if (save_plots && !is.null(output_dir) && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Plot each component
  for (i in 1:n_components) {
    comp_name <- comp_names[i]

    if (save_plots && !is.null(output_dir)) {
      filename <- file.path(output_dir, paste0("conditional_reliability_",
                                               gsub("[^A-Za-z0-9]", "_", comp_name), ".png"))
    } else {
      filename <- NULL
    }

    message(sprintf("Creating plot for component: %s", comp_name))

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

  if (save_plots && !is.null(output_dir)) {
    message(sprintf("All plots saved to: %s", output_dir))
  }

  invisible(NULL)
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
#'   Default \code{NULL} (no file written). Only writes when both
#'   \code{save_plot = TRUE} and \code{filename} is non-\code{NULL}.
#' @param ... Additional graphical arguments.
#' @importFrom RColorBrewer brewer.pal
#' @return Invisibly returns \code{NULL}. Called for its side effect of
#'   producing a comparative reliability plot across model components.
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
#'   plot_components_comparison(cond_rel)
#' }
#'
#' @family conditional reliability functions
#' @export
plot_components_comparison <- function(results,
                                       include_ci = TRUE,
                                       show_optimal_points = TRUE,
                                       add_reference_lines = TRUE,
                                       color_palette = "Set2",
                                       save_plot = FALSE,
                                       filename = NULL,
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
  if (save_plot && !is.null(filename)) {
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
              col = adjustcolor(colors[i], alpha.f = 0.2),
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

  if (save_plot && !is.null(filename)) {
    dev.off()
    message("Comparison plot saved as: ", filename)
  }

  invisible(NULL)
}

# Helper function for null coalescing
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}
