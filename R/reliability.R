#' @title
#' Marginal reliability
#' @description
#'
#' Estimate the the marginal reliability of the GMLTM.
#'
#' @usage
#'
#' reliability(fit)
#'
#' @param fit MLTM object.
#'
#' @details \code{reliability} estimates a ...
#'
#' @return A number denoting the reliability estimate.
#'
#' @references
#'
#' Ramírez, E.S.; Jiménez, M.; Franco, V.R.; Alvarado, J.M. Delving into the Complexity of Analogical Reasoning: A Detailed Exploration with the Generalized Multicomponent Latent Trait Model for Diagnosis. \emph{J. Intell.} 2024, 12, 67. https://doi.org/10.3390/jintelligence12070067
#'
#' @export
reliability <- function(fit) {

  if (inherits(fit, "LLTM")) {
    EAP_theta <- matrix(fit$EAP$theta)
    VAR_theta <- matrix(apply(fit$posterior$theta, MARGIN = 2, FUN = var))
  } else {
    EAP_theta <- fit$EAP$theta
    VAR_theta <- apply(fit$posterior$theta, MARGIN = c(2, 3), FUN = var)
  }
  n <- nrow(EAP_theta)
  q <- ncol(EAP_theta)
  COV <- cov(EAP_theta)
  e <- colMeans(VAR_theta)
  s <- diag(COV)
  rxx <- s/(s+e)

  return(rxx)

}


#' Enhanced Reliability Analysis for GMLTM-D Models
#'
#' @description
#' Provides comprehensive reliability analysis for General Multicomponent Latent Trait Models
#' for Diagnosis (GMLTM-D) using Bayesian posterior distributions. Optimized for speed
#' and minimal dependencies.
#'
#' @param fit A fitted GMLTM, MLTM, or LLTM model object containing posterior samples
#'   of theta parameters.
#' @param include_conditional Logical. Whether to compute conditional reliability
#'   estimates across ability levels. Default is FALSE for speed.
#' @param include_hierarchical Logical. Whether to compute hierarchical reliability
#'   for the general factor. Default is TRUE.
#' @param include_comparisons Logical. Whether to perform Bayesian comparisons
#'   between components. Default is TRUE.
#' @param n_samples Integer. Number of posterior samples to use (for speed control).
#'   If NULL, uses all available samples.
#'
#' @return An object of class \code{enhanced_mltm_reliability}.
#'
#' @export
enhanced_mltm_reliability <- function(fit,
                                      include_conditional = FALSE,  # Default FALSE for speed
                                      include_hierarchical = TRUE,
                                      include_comparisons = TRUE,
                                      n_samples = NULL) {

  # Input validation (simplified)
  if (!inherits(fit, c("GMLTM", "LLTM", "MLTM"))) {
    stop("The input must be a fitted GMLTM, MLTM, or LLTM model.")
  }

  if (is.null(fit$posterior$theta) || is.null(fit$EAP$theta)) {
    stop("The model does not contain valid posterior estimates for theta.")
  }

  # Handle different model structures
  if (inherits(fit, "LLTM")) {
    theta_samples <- array(fit$posterior$theta,
                           dim = c(nrow(fit$posterior$theta),
                                   ncol(fit$posterior$theta), 1))
    EAP_theta <- matrix(fit$EAP$theta, ncol = 1)
  } else {
    theta_samples <- fit$posterior$theta
    EAP_theta <- fit$EAP$theta
  }

  # Subsample for speed if requested
  if (!is.null(n_samples) && n_samples < dim(theta_samples)[1]) {
    sample_indices <- sample(1:dim(theta_samples)[1], n_samples)
    theta_samples <- theta_samples[sample_indices, , , drop = FALSE]
  }

  results <- list()

  # Store model info
  results$model_info <- list(
    model_class = class(fit)[1],
    n_persons = dim(theta_samples)[2],
    n_components = ifelse(length(dim(theta_samples)) == 3,
                          dim(theta_samples)[3], 1),
    n_iterations = dim(theta_samples)[1]
  )

  # 1. Marginal reliability (optimized)
  results$marginal_robust <- bayesian_reliability_fast(theta_samples, EAP_theta)

  # 2. Hierarchical reliability
  if (include_hierarchical && results$model_info$n_components > 1) {
    results$hierarchical <- hierarchical_reliability_fast(theta_samples)
  }

  # 3. Component comparisons
  if (include_comparisons && results$model_info$n_components > 1) {
    results$comparisons <- reliability_comparison_fast(results$marginal_robust)
  }

  # 4. Conditional reliability (optional, computationally expensive)
  if (include_conditional && results$model_info$n_persons >= 20) {
    results$conditional <- conditional_reliability_fast(theta_samples, EAP_theta)
  }

  # 5. Basic diagnostics (without coda dependency)
  results$diagnostics <- compute_basic_diagnostics(theta_samples)

  # Set class
  class(results) <- "enhanced_mltm_reliability"
  attr(results, "call") <- match.call()

  return(results)
}

#' Fast Bayesian Marginal Reliability
#'
#' @description
#' Optimized computation of marginal reliability with uncertainty quantification.
#'
#' @param theta_samples Array of posterior samples [iter, person, component].
#' @param EAP_theta Matrix of EAP estimates [person, component].
#'
#' @return List of reliability estimates per component.
#'
#' @keywords internal
bayesian_reliability_fast <- function(theta_samples, EAP_theta) {

  n_iter <- dim(theta_samples)[1]
  n_person <- dim(theta_samples)[2]
  n_comp <- dim(theta_samples)[3]

  # Pre-allocate results
  reliability_samples <- matrix(NA, nrow = n_iter, ncol = n_comp)

  # Vectorized computation across iterations
  for (c in 1:n_comp) {
    # Between-person variance for each iteration (vectorized)
    var_between <- apply(theta_samples[, , c], 1, var)

    # Within-person error variance (pre-computed)
    var_within_person <- apply(theta_samples[, , c], 2, var)
    var_within <- rep(mean(var_within_person), n_iter)

    # Reliability computation (vectorized)
    total_var <- var_between + var_within
    reliability_samples[, c] <- ifelse(total_var > 0, var_between / total_var, NA)
  }

  # Summarize results
  results <- list()
  comp_names <- if (!is.null(colnames(EAP_theta))) {
    colnames(EAP_theta)
  } else {
    paste0("Component_", 1:n_comp)
  }

  for (c in 1:n_comp) {
    valid_samples <- reliability_samples[!is.na(reliability_samples[, c]), c]

    if (length(valid_samples) > 0) {
      # Basic quantiles (fast)
      quants <- quantile(valid_samples, c(0.025, 0.25, 0.5, 0.75, 0.975))

      # HDI approximation (fast)
      hdi_95 <- c(quants[1], quants[5])  # Approximation using quantiles
      hdi_90 <- quantile(valid_samples, c(0.05, 0.95))

      results[[comp_names[c]]] <- list(
        mean = mean(valid_samples),
        sd = sd(valid_samples),
        median = quants[3],
        quantiles = quants,
        hdi_90 = hdi_90,
        hdi_95 = hdi_95,
        samples = valid_samples,
        n_samples = length(valid_samples)
      )
    }
  }

  return(results)
}

#' Fast Hierarchical Reliability
#'
#' @description
#' Optimized hierarchical reliability computation.
#'
#' @param theta_samples Array of posterior samples.
#'
#' @return List with component and general factor reliability.
#'
#' @keywords internal
hierarchical_reliability_fast <- function(theta_samples) {

  n_iter <- dim(theta_samples)[1]
  n_person <- dim(theta_samples)[2]
  n_comp <- dim(theta_samples)[3]

  # General factor computation (vectorized)
  general_theta <- apply(theta_samples, c(1, 2), mean)  # Average across components

  # General factor reliability
  var_between_general <- apply(general_theta, 1, var)

  # Error propagation (simplified)
  component_vars <- apply(theta_samples, c(1, 2), var)  # Variance across components per person
  var_error_general <- apply(component_vars, 1, mean) / n_comp

  # General reliability samples
  total_var_general <- var_between_general + var_error_general
  general_reliability_samples <- ifelse(total_var_general > 0,
                                        var_between_general / total_var_general, NA)

  # Use existing component reliability + add general factor
  component_results <- bayesian_reliability_fast(theta_samples,
                                                 matrix(apply(theta_samples, c(2, 3), mean),
                                                        ncol = n_comp))

  # Add general factor
  valid_general <- general_reliability_samples[!is.na(general_reliability_samples)]

  if (length(valid_general) > 0) {
    quants <- quantile(valid_general, c(0.025, 0.25, 0.5, 0.75, 0.975))

    component_results$General_Factor <- list(
      mean = mean(valid_general),
      sd = sd(valid_general),
      median = quants[3],
      quantiles = quants,
      hdi_90 = quantile(valid_general, c(0.05, 0.95)),
      hdi_95 = c(quants[1], quants[5]),
      samples = valid_general,
      n_samples = length(valid_general)
    )
  }

  return(component_results)
}

#' Fast Reliability Comparisons
#'
#' @description
#' Optimized pairwise comparisons of reliability between components.
#'
#' @param marginal_results List of marginal reliability results.
#'
#' @return List of comparison results.
#'
#' @keywords internal
reliability_comparison_fast <- function(marginal_results) {

  n_comp <- length(marginal_results)
  if (n_comp < 2) return(NULL)

  # Extract samples (ensure equal length)
  reliability_samples <- lapply(marginal_results, function(x) x$samples)
  min_length <- min(sapply(reliability_samples, length))
  reliability_matrix <- sapply(reliability_samples, function(x) x[1:min_length])

  comp_names <- names(marginal_results)
  comparisons <- list()

  # Pairwise comparisons (optimized)
  for (i in 1:(n_comp - 1)) {
    for (j in (i + 1):n_comp) {
      diff_samples <- reliability_matrix[, i] - reliability_matrix[, j]

      comparison_name <- paste0(comp_names[i], "_vs_", comp_names[j])

      # Basic statistics (fast)
      mean_diff <- mean(diff_samples)
      sd_diff <- sd(diff_samples)
      prob_greater <- mean(diff_samples > 0)

      # HDI approximation
      hdi_diff <- quantile(diff_samples, c(0.025, 0.975))

      comparisons[[comparison_name]] <- list(
        component_1 = comp_names[i],
        component_2 = comp_names[j],
        mean_difference = mean_diff,
        sd_difference = sd_diff,
        prob_greater = prob_greater,
        prob_substantially_greater = mean(diff_samples > 0.05),
        hdi_difference = hdi_diff,
        effect_size = mean_diff / sd_diff
      )
    }
  }

  return(comparisons)
}

#' Fast Conditional Reliability (Simplified)
#'
#' @description
#' Simplified conditional reliability using quantile-based approach.
#'
#' @param theta_samples Array of posterior samples.
#' @param EAP_theta Matrix of EAP estimates.
#' @param n_bins Integer. Number of ability bins (default 5 for speed).
#'
#' @return List of conditional reliability by bins.
#'
#' @keywords internal
conditional_reliability_fast <- function(theta_samples, EAP_theta, n_bins = 5) {

  n_comp <- dim(theta_samples)[3]
  results <- list()

  for (c in 1:n_comp) {
    # Create ability bins based on EAP estimates
    theta_eap <- EAP_theta[, c]
    breaks <- quantile(theta_eap, seq(0, 1, length.out = n_bins + 1))
    bins <- cut(theta_eap, breaks, include.lowest = TRUE)

    bin_results <- list()

    for (b in 1:n_bins) {
      bin_indices <- which(bins == levels(bins)[b])

      if (length(bin_indices) >= 5) {  # Minimum for stable estimates
        # Subsample theta for this bin
        theta_bin <- theta_samples[, bin_indices, c]

        # Reliability for this bin (simplified)
        bin_reliability <- apply(theta_bin, 1, function(theta_iter) {
          var_between <- var(theta_iter)
          var_within <- mean(apply(theta_samples[, bin_indices, c], 2, var))

          if (var_between + var_within > 0) {
            return(var_between / (var_between + var_within))
          } else {
            return(NA)
          }
        })

        valid_rel <- bin_reliability[!is.na(bin_reliability)]

        if (length(valid_rel) > 0) {
          bin_results[[b]] <- list(
            theta_range = range(theta_eap[bin_indices]),
            theta_center = mean(theta_eap[bin_indices]),
            mean = mean(valid_rel),
            sd = sd(valid_rel),
            quantiles = quantile(valid_rel, c(0.025, 0.5, 0.975)),
            n_persons = length(bin_indices)
          )
        }
      }
    }

    comp_name <- ifelse(!is.null(colnames(EAP_theta)),
                        colnames(EAP_theta)[c],
                        paste0("Component_", c))
    results[[comp_name]] <- bin_results
  }

  return(results)
}

#' Basic Diagnostics (No External Dependencies)
#'
#' @description
#' Computes basic MCMC diagnostics without external packages.
#'
#' @param theta_samples Array of posterior samples.
#'
#' @return List of diagnostic statistics.
#'
#' @keywords internal
compute_basic_diagnostics <- function(theta_samples) {

  diagnostics <- list()

  # Effective sample size approximation (simple autocorrelation-based)
  n_eff_approx <- apply(theta_samples, c(2, 3), function(x) {
    n <- length(x)
    if (n <= 1) return(NA)

    # Lag-1 autocorrelation
    rho1 <- cor(x[-1], x[-n], use = "complete.obs")
    if (is.na(rho1)) return(n)

    # Simple effective size approximation
    eff_size <- n / (1 + 2 * pmax(0, rho1))
    return(eff_size)
  })

  diagnostics$n_effective_approx <- n_eff_approx

  # Monte Carlo Standard Error
  diagnostics$mcse <- apply(theta_samples, c(2, 3),
                            function(x) sd(x, na.rm = TRUE) / sqrt(length(x)))

  # Autocorrelation diagnostics
  diagnostics$autocorr_lag1 <- apply(theta_samples, c(2, 3), function(x) {
    n <- length(x)
    if (n <= 1) return(NA)
    cor(x[-1], x[-n], use = "complete.obs")
  })

  # Geweke diagnostic approximation (split-chain comparison)
  diagnostics$geweke_approx <- apply(theta_samples, c(2, 3), function(x) {
    n <- length(x)
    if (n < 20) return(NA)

    # First 10% vs last 50%
    first_part <- x[1:floor(0.1 * n)]
    last_part <- x[ceiling(0.5 * n):n]

    mean_diff <- mean(first_part) - mean(last_part)
    se_diff <- sqrt(var(first_part)/length(first_part) + var(last_part)/length(last_part))

    if (se_diff > 0) {
      z_score <- abs(mean_diff / se_diff)
      return(z_score)
    } else {
      return(0)
    }
  })

  return(diagnostics)
}

#' Print Method for Enhanced MLTM Reliability
#'
#' @param x An object of class \code{"enhanced_mltm_reliability"}.
#' @param digits Integer. Number of decimal places to display. Default is \code{3}.
#' @param ... Currently unused.
#' @export
print.enhanced_mltm_reliability <- function(x, digits = 3, ...) {

  cat("Enhanced GMLTM-D Reliability Analysis\n")
  cat("====================================\n\n")

  # Model information
  cat("Model Information:\n")
  cat(sprintf("  Model class: %s\n", x$model_info$model_class))
  cat(sprintf("  Persons: %d\n", x$model_info$n_persons))
  cat(sprintf("  Components: %d\n", x$model_info$n_components))
  cat(sprintf("  MCMC iterations used: %d\n\n", x$model_info$n_iterations))

  # Marginal reliability
  if (!is.null(x$marginal_robust)) {
    cat("Marginal Reliability (Bayesian):\n")
    cat("---------------------------------\n")

    for (i in 1:length(x$marginal_robust)) {
      rel <- x$marginal_robust[[i]]
      cat(sprintf("  %s:\n", names(x$marginal_robust)[i]))
      cat(sprintf("    Mean: %.*f (SD: %.*f)\n", digits, rel$mean, digits, rel$sd))
      cat(sprintf("    95%% CI: [%.*f, %.*f]\n",
                  digits, rel$hdi_95[1], digits, rel$hdi_95[2]))
      cat(sprintf("    Samples: %d\n\n", rel$n_samples))
    }
  }

  # General factor
  if (!is.null(x$hierarchical) && "General_Factor" %in% names(x$hierarchical)) {
    cat("General Factor Reliability:\n")
    cat("---------------------------\n")
    rel_gen <- x$hierarchical$General_Factor
    cat(sprintf("  Mean: %.*f (SD: %.*f)\n", digits, rel_gen$mean, digits, rel_gen$sd))
    cat(sprintf("  95%% CI: [%.*f, %.*f]\n\n",
                digits, rel_gen$hdi_95[1], digits, rel_gen$hdi_95[2]))
  }

  # Comparisons (top 3 only for brevity)
  if (!is.null(x$comparisons)) {
    cat("Component Comparisons (Top 3):\n")
    cat("------------------------------\n")

    # Sort by absolute difference
    comp_sorted <- x$comparisons[order(sapply(x$comparisons, function(y) abs(y$mean_difference)),
                                       decreasing = TRUE)]

    for (i in 1:min(3, length(comp_sorted))) {
      comp_name <- names(comp_sorted)[i]
      diff <- comp_sorted[[i]]
      cat(sprintf("  %s vs %s:\n", diff$component_1, diff$component_2))
      cat(sprintf("    Difference: %.*f, P(C1>C2): %.*f\n",
                  digits, diff$mean_difference, digits, diff$prob_greater))

      # Simple interpretation
      if (diff$prob_greater > 0.95) {
        cat("    -> Strong evidence C1 more reliable\n")
      } else if (diff$prob_greater < 0.05) {
        cat("    -> Strong evidence C2 more reliable\n")
      } else {
        cat("    -> Inconclusive evidence\n")
      }
      cat("\n")
    }

    if (length(x$comparisons) > 3) {
      cat(sprintf("  ... and %d more comparisons\n\n", length(x$comparisons) - 3))
    }
  }

  # Diagnostic warnings
  if (!is.null(x$diagnostics$n_effective_approx)) {
    min_eff <- min(x$diagnostics$n_effective_approx, na.rm = TRUE)
    if (min_eff < 50) {
      cat("[WARNING]  WARNING: Low effective sample size detected (<50)\n")
      cat("   Consider running longer MCMC chains\n\n")
    }
  }

  if (!is.null(x$diagnostics$geweke_approx)) {
    max_geweke <- max(x$diagnostics$geweke_approx, na.rm = TRUE)
    if (max_geweke > 2) {
      cat("[WARNING]  WARNING: Potential convergence issues detected\n")
      cat("   Consider checking chain mixing\n\n")
    }
  }

  invisible(x)
}

#' Simple Plot Method (Minimal Dependencies)
#'
#' @description
#' Creates basic plots using base R (no ggplot2 dependency).
#'
#' @param x Object of class \code{enhanced_mltm_reliability}.
#' @param type Character. Type of plot: "marginal", "conditional", "comparison".
#' @param component Integer or character. Specific component for conditional plots.
#' @param ... Additional plotting parameters.
#'
#' @export
plot.enhanced_mltm_reliability <- function(x, type = "marginal", component = NULL, ...) {

  if (type == "marginal" && !is.null(x$marginal_robust)) {
    plot_marginal_base(x$marginal_robust, ...)

  } else if (type == "conditional" && !is.null(x$conditional)) {
    plot_conditional_base(x$conditional, component, ...)

  } else if (type == "comparison" && !is.null(x$comparisons)) {
    plot_comparison_base(x$comparisons, ...)

  } else {
    cat("Plot type not available or no data for requested plot.\n")
    cat("Available types:",
        paste(c("marginal"[!is.null(x$marginal_robust)],
                "conditional"[!is.null(x$conditional)],
                "comparison"[!is.null(x$comparisons)]), collapse = ", "))
  }
}

#' Plot Marginal Reliability (Base R)
#'
#' @keywords internal
plot_marginal_base <- function(marginal_data, ...) {

  comp_names <- names(marginal_data)
  n_comp <- length(comp_names)

  # Extract data
  means <- sapply(marginal_data, function(x) x$mean)
  lower_95 <- sapply(marginal_data, function(x) x$hdi_95[1])
  upper_95 <- sapply(marginal_data, function(x) x$hdi_95[2])
  lower_90 <- sapply(marginal_data, function(x) x$hdi_90[1])
  upper_90 <- sapply(marginal_data, function(x) x$hdi_90[2])

  # Create plot
  x_pos <- 1:n_comp

  plot(x_pos, means,
       ylim = c(0, 1),
       xlim = c(0.5, n_comp + 0.5),
       xlab = "Component",
       ylab = "Reliability",
       main = "Bayesian Marginal Reliability",
       pch = 19,
       cex = 1.5,
       col = "darkblue",
       xaxt = "n",
       ...)

  # Add error bars
  segments(x_pos, lower_95, x_pos, upper_95, col = "darkblue", lwd = 1)
  segments(x_pos, lower_90, x_pos, upper_90, col = "blue", lwd = 3)

  # Add reference lines
  abline(h = c(0.6, 0.7, 0.8, 0.9), lty = 2, col = "gray60")

  # Labels
  axis(1, at = x_pos, labels = comp_names, las = 2)

  # Legend
  legend("bottomright",
         legend = c("Mean", "90% CI", "95% CI"),
         col = c("darkblue", "blue", "darkblue"),
         lwd = c(0, 3, 1),
         pch = c(19, NA, NA),
         bty = "n")
}

#' Plot Conditional Reliability (Base R)
#'
#' @keywords internal
plot_conditional_base <- function(conditional_data, component = NULL, ...) {

  # Select component
  if (!is.null(component)) {
    if (is.numeric(component)) {
      comp_data <- conditional_data[[component]]
      comp_name <- names(conditional_data)[component]
    } else {
      comp_data <- conditional_data[[component]]
      comp_name <- component
    }
  } else {
    comp_data <- conditional_data[[1]]
    comp_name <- names(conditional_data)[1]
  }

  # Extract valid data
  valid_data <- comp_data[!sapply(comp_data, is.null)]
  if (length(valid_data) == 0) {
    cat("No valid conditional data for plotting.\n")
    return(invisible())
  }

  # Extract plotting data
  theta_centers <- sapply(valid_data, function(x) x$theta_center)
  means <- sapply(valid_data, function(x) x$mean)
  lower <- sapply(valid_data, function(x) x$quantiles[1])
  upper <- sapply(valid_data, function(x) x$quantiles[3])

  # Sort by theta
  ord <- order(theta_centers)
  theta_centers <- theta_centers[ord]
  means <- means[ord]
  lower <- lower[ord]
  upper <- upper[ord]

  # Plot
  plot(theta_centers, means,
       type = "l",
       ylim = c(0, 1),
       xlab = expression(theta ~ "(Ability Level)"),
       ylab = "Conditional Reliability",
       main = paste("Conditional Reliability -", comp_name),
       col = "darkblue",
       lwd = 2,
       ...)

  # Add confidence bands
  polygon(c(theta_centers, rev(theta_centers)),
          c(lower, rev(upper)),
          col = adjustcolor("lightblue", alpha = 0.3),
          border = NA)

  # Re-draw main line
  lines(theta_centers, means, col = "darkblue", lwd = 2)
  points(theta_centers, means, pch = 19, col = "darkblue")

  # Reference lines
  abline(h = c(0.6, 0.7, 0.8, 0.9), lty = 2, col = "gray60")

  # Grid
  grid(col = "gray90")
}

#' Plot Comparisons (Base R)
#'
#' @keywords internal
plot_comparison_base <- function(comparison_data, ...) {

  n_comp <- length(comparison_data)
  if (n_comp == 0) return(invisible())

  # Extract data
  comp_names <- names(comparison_data)
  mean_diffs <- sapply(comparison_data, function(x) x$mean_difference)
  lower_diffs <- sapply(comparison_data, function(x) x$hdi_difference[1])
  upper_diffs <- sapply(comparison_data, function(x) x$hdi_difference[2])

  # Colors based on significance
  colors <- ifelse(lower_diffs > 0, "darkgreen",
                   ifelse(upper_diffs < 0, "darkred", "darkblue"))

  # Horizontal plot
  y_pos <- 1:n_comp

  plot(mean_diffs, y_pos,
       xlim = range(c(lower_diffs, upper_diffs)) * 1.1,
       ylim = c(0.5, n_comp + 0.5),
       xlab = "Reliability Difference",
       ylab = "",
       main = "Component Reliability Comparisons",
       pch = 19,
       cex = 1.5,
       col = colors,
       yaxt = "n",
       ...)

  # Error bars
  segments(lower_diffs, y_pos, upper_diffs, y_pos, col = colors, lwd = 2)

  # Zero line
  abline(v = 0, lty = 1, col = "black", lwd = 1)

  # Labels
  axis(2, at = y_pos, labels = comp_names, las = 2)

  # Grid
  grid(col = "gray90")

  # Legend
  legend("topright",
         legend = c("C1 > C2", "C1 < C2", "Uncertain"),
         col = c("darkgreen", "darkred", "darkblue"),
         pch = 19,
         bty = "n")
}

#' Quick Reliability Check
#'
#' @description
#' Ultra-fast reliability check for initial assessment.
#'
#' @param fit Fitted model object.
#' @param n_samples_quick Integer. Number of samples for quick analysis (default 500).
#'
#' @return Named vector of reliability estimates.
#'
#' @examples
#' \dontrun{
#' # Very fast initial check
#' quick_rel <- quick_reliability_check(fit)
#' print(quick_rel)
#' }
#'
#' @export
quick_reliability_check <- function(fit, n_samples_quick = 500) {

  # Input validation
  if (!inherits(fit, c("GMLTM", "LLTM", "MLTM"))) {
    stop("Invalid model type.")
  }

  # Handle model structure
  if (inherits(fit, "LLTM")) {
    theta_samples <- array(fit$posterior$theta[1:min(n_samples_quick, nrow(fit$posterior$theta)), ],
                           dim = c(min(n_samples_quick, nrow(fit$posterior$theta)),
                                   ncol(fit$posterior$theta), 1))
    EAP_theta <- matrix(fit$EAP$theta, ncol = 1)
  } else {
    n_use <- min(n_samples_quick, dim(fit$posterior$theta)[1])
    theta_samples <- fit$posterior$theta[1:n_use, , , drop = FALSE]
    EAP_theta <- fit$EAP$theta
  }

  n_comp <- dim(theta_samples)[3]

  # Very simple reliability computation
  reliability_means <- numeric(n_comp)

  for (c in 1:n_comp) {
    # Simple variance decomposition
    var_between <- var(apply(theta_samples[, , c], 2, mean))  # Between-person variance
    var_within <- mean(apply(theta_samples[, , c], 2, var))   # Average within-person variance

    if (var_between + var_within > 0) {
      reliability_means[c] <- var_between / (var_between + var_within)
    } else {
      reliability_means[c] <- NA
    }
  }

  # Component names
  comp_names <- if (!is.null(colnames(EAP_theta))) {
    colnames(EAP_theta)
  } else {
    paste0("Component_", 1:n_comp)
  }

  names(reliability_means) <- comp_names

  return(reliability_means)
}

# ================================================================================
# SIMPLIFIED USE EXAMPLE
# ================================================================================

#' Step-by-step Usage Example
#'
#' @description
#' Demonstration function to show how to use the optimized functions.
#'
#' @param fit Fitted model
#'
#' @examples
#' \dontrun{
#' # Basic and fast usage
#' demo_reliability_analysis(fit)
#' }
#'
#' @export
demo_reliability_analysis <- function(fit) {

  cat("=== STEP-BY-STEP RELIABILITY ANALYSIS ===\n\n")

  # Step 1: Quick check
  cat("1. QUICK CHECK (ultra-fast)\n")
  cat("-------------------------------------\n")

  quick_rel <- quick_reliability_check(fit)
  print(round(quick_rel, 3))

  cat("\nQuick interpretation:\n")
  for (i in seq_along(quick_rel)) {
    rel_val <- quick_rel[i]
    comp_name <- names(quick_rel)[i]

    if (rel_val >= 0.8) {
      cat(sprintf("  %s: GOOD (%.3f)\n", comp_name, rel_val))
    } else if (rel_val >= 0.7) {
      cat(sprintf("  %s: ACCEPTABLE (%.3f)\n", comp_name, rel_val))
    } else if (rel_val >= 0.6) {
      cat(sprintf("  %s: QUESTIONABLE (%.3f)\n", comp_name, rel_val))
    } else {
      cat(sprintf("  %s: POOR (%.3f)\n", comp_name, rel_val))
    }
  }

  cat("\n" %+% paste(rep("=", 50), collapse = "") %+% "\n\n")

  # Step 2: Detailed analysis (optimized)
  cat("2. DETAILED ANALYSIS (optimized)\n")
  cat("----------------------------------\n")

  rel_detailed <- enhanced_mltm_reliability(
    fit,
    include_conditional = FALSE,    # Fast
    include_hierarchical = TRUE,
    include_comparisons = TRUE,
    n_samples = 1000               # Limit samples for speed
  )

  print(rel_detailed)

  cat("\n" %+% paste(rep("=", 50), collapse = "") %+% "\n\n")

  # Step 3: Basic visualization
  cat("3. VISUALIZATION\n")
  cat("----------------\n")
  cat("Creating marginal reliability plot...\n")

  plot(rel_detailed, type = "marginal")

  if (!is.null(rel_detailed$comparisons)) {
    cat("Press Enter to see component comparisons...")
    readline()
    plot(rel_detailed, type = "comparison")
  }

  return(rel_detailed)
}

# ================================================================================
# ADDITIONAL UTILITY FUNCTIONS
# ================================================================================

#' Check Data Quality for Reliability Analysis
#'
#' @description
#' Verifies if the data is suitable for robust reliability analysis.
#'
#' @param fit Fitted model
#'
#' @return List with data quality diagnostics
#'
#' @export
check_reliability_data_quality <- function(fit) {

  # Basic checks
  checks <- list()

  # 1. Data structure
  if (inherits(fit, "LLTM")) {
    theta_samples <- array(fit$posterior$theta,
                           dim = c(nrow(fit$posterior$theta),
                                   ncol(fit$posterior$theta), 1))
  } else {
    theta_samples <- fit$posterior$theta
  }

  checks$n_iterations <- dim(theta_samples)[1]
  checks$n_persons <- dim(theta_samples)[2]
  checks$n_components <- dim(theta_samples)[3]

  # 2. Sample size
  checks$sample_size_adequate <- checks$n_persons >= 100
  checks$mcmc_length_adequate <- checks$n_iterations >= 1000

  # 3. Data variability
  eap_theta <- if (inherits(fit, "LLTM")) matrix(fit$EAP$theta, ncol = 1) else fit$EAP$theta

  checks$theta_variability <- apply(eap_theta, 2, function(x) {
    list(
      range = diff(range(x, na.rm = TRUE)),
      sd = sd(x, na.rm = TRUE),
      adequate_variability = sd(x, na.rm = TRUE) > 0.5
    )
  })

  # 4. Basic convergence
  autocorrs <- apply(theta_samples, c(2, 3), function(x) {
    if (length(x) > 1) {
      cor(x[-1], x[-length(x)], use = "complete.obs")
    } else {
      NA
    }
  })

  checks$high_autocorrelation <- mean(autocorrs > 0.8, na.rm = TRUE) > 0.1

  # 5. Extreme values
  checks$extreme_values <- apply(eap_theta, 2, function(x) {
    q99 <- quantile(x, 0.99, na.rm = TRUE)
    q01 <- quantile(x, 0.01, na.rm = TRUE)
    sum(x > q99 | x < q01, na.rm = TRUE) / length(x)
  })

  # 6. Recommendations
  checks$recommendations <- list()

  if (!checks$sample_size_adequate) {
    checks$recommendations <- c(checks$recommendations,
                                "Small sample size (<100). Interpret results with caution.")
  }

  if (!checks$mcmc_length_adequate) {
    checks$recommendations <- c(checks$recommendations,
                                "Short MCMC chains (<1000). Consider more iterations.")
  }

  if (checks$high_autocorrelation) {
    checks$recommendations <- c(checks$recommendations,
                                "High autocorrelation detected. Check MCMC convergence.")
  }

  if (any(sapply(checks$theta_variability, function(x) !x$adequate_variability))) {
    checks$recommendations <- c(checks$recommendations,
                                "Low variability in some components. Check ability range.")
  }

  class(checks) <- "reliability_data_quality"
  return(checks)
}

#' Print Method for Data Quality Diagnostics
#'
#' @param x An object of class \code{"reliability_data_quality"}.
#' @param ... Currently unused.
#' @export
print.reliability_data_quality <- function(x, ...) {

  cat("DATA QUALITY DIAGNOSTICS FOR RELIABILITY\n")
  cat("=======================================\n\n")

  cat("Data Characteristics:\n")
  cat(sprintf("  Persons: %d %s\n", x$n_persons,
              ifelse(x$sample_size_adequate, "[OK]", "[WARNING]")))
  cat(sprintf("  Components: %d\n", x$n_components))
  cat(sprintf("  MCMC iterations: %d %s\n", x$n_iterations,
              ifelse(x$mcmc_length_adequate, "[OK]", "[WARNING]")))

  cat("\nVariability by Component:\n")
  for (i in seq_along(x$theta_variability)) {
    var_info <- x$theta_variability[[i]]
    comp_name <- names(x$theta_variability)[i] %||% paste("Component", i)
    cat(sprintf("  %s: Range=%.2f, SD=%.2f %s\n",
                comp_name, var_info$range, var_info$sd,
                ifelse(var_info$adequate_variability, "[OK]", "[WARNING]")))
  }

  if (x$high_autocorrelation) {
    cat("\n[WARNING]  WARNING: High autocorrelation detected\n")
  }

  if (length(x$recommendations) > 0) {
    cat("\nRECOMMENDATIONS:\n")
    for (i in seq_along(x$recommendations)) {
      cat(sprintf("%d. %s\n", i, x$recommendations[[i]]))
    }
  } else {
    cat("\n[OK] Data suitable for robust reliability analysis\n")
  }

  cat("\n")
  invisible(x)
}

#' @noRd
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' @noRd
`%+%` <- function(x, y) {
  paste0(x, y)
}

# ================================================================================
# STEP-BY-STEP USAGE INSTRUCTIONS
# ================================================================================

#' Usage Instructions
#'
#' @description
#' Prints detailed usage instructions.
#'
#' @export
reliability_usage_instructions <- function() {

  cat("USAGE INSTRUCTIONS - OPTIMIZED RELIABILITY ANALYSIS\n")
  cat("===================================================\n\n")

  cat("STEP 1: Quick check (5-10 seconds)\n")
  cat("-----------------------------------\n")
  cat("quick_rel <- quick_reliability_check(your_model)\n")
  cat("print(quick_rel)\n\n")

  cat("STEP 2: Check data quality\n")
  cat("---------------------------\n")
  cat("data_quality <- check_reliability_data_quality(your_model)\n")
  cat("print(data_quality)\n\n")

  cat("STEP 3: Optimized full analysis (30-60 seconds)\n")
  cat("----------------------------------------------\n")
  cat("rel_analysis <- enhanced_mltm_reliability(your_model,\n")
  cat("                include_conditional = FALSE,  # Fast\n")
  cat("                n_samples = 1000)            # Limit samples\n")
  cat("print(rel_analysis)\n\n")

  cat("STEP 4: Basic visualization\n")
  cat("----------------------------\n")
  cat("plot(rel_analysis, type = 'marginal')\n")
  cat("plot(rel_analysis, type = 'comparison')  # If multiple components\n\n")

  cat("STEP 5: Conditional analysis (OPTIONAL - slower)\n")
  cat("--------------------------------------------------\n")
  cat("rel_conditional <- enhanced_mltm_reliability(your_model,\n")
  cat("                   include_conditional = TRUE,\n")
  cat("                   n_samples = 500)  # Fewer samples for speed\n")
  cat("plot(rel_conditional, type = 'conditional')\n\n")

  cat("AUTOMATIC DEMONSTRATION:\n")
  cat("------------------------\n")
  cat("demo_reliability_analysis(your_model)  # Step-by-step guide\n\n")

  cat("TROUBLESHOOTING:\n")
  cat("----------------\n")
  cat("1. Error 'could not find function': All dependencies are eliminated\n")
  cat("2. High time: Use n_samples to limit MCMC samples\n")
  cat("3. No conditional data: include_conditional = FALSE by default\n")
  cat("4. Insufficient memory: Use quick_reliability_check() for basic analysis\n\n")

  invisible()
}

# ================================================================================
# EXPORT FUNCTIONS
# ================================================================================

#' Export Reliability Results
#'
#' @description
#' Exports results in tabular format for publications.
#'
#' @param reliability_obj Object of class enhanced_mltm_reliability
#' @param file_name File name (optional)
#'
#' @return data.frame with tabulated results
#'
#' @export
export_reliability_results <- function(reliability_obj, file_name = NULL) {

  if (!inherits(reliability_obj, "enhanced_mltm_reliability")) {
    stop("The object must be of class 'enhanced_mltm_reliability'")
  }

  if (is.null(reliability_obj$marginal_robust)) {
    stop("No marginal reliability results available")
  }

  # Create results table
  results_table <- data.frame(
    Component = names(reliability_obj$marginal_robust),
    Mean = sapply(reliability_obj$marginal_robust, function(x) round(x$mean, 3)),
    SD = sapply(reliability_obj$marginal_robust, function(x) round(x$sd, 3)),
    CI_Lower = sapply(reliability_obj$marginal_robust, function(x) round(x$hdi_95[1], 3)),
    CI_Upper = sapply(reliability_obj$marginal_robust, function(x) round(x$hdi_95[2], 3)),
    N_Samples = sapply(reliability_obj$marginal_robust, function(x) x$n_samples),
    Quality = sapply(reliability_obj$marginal_robust, function(x) {
      if (x$mean >= 0.8) "Good"
      else if (x$mean >= 0.7) "Acceptable"
      else if (x$mean >= 0.6) "Questionable"
      else "Poor"
    }),
    stringsAsFactors = FALSE
  )

  # Add general factor if exists
  if (!is.null(reliability_obj$hierarchical) &&
      "General_Factor" %in% names(reliability_obj$hierarchical)) {

    gen_factor <- reliability_obj$hierarchical$General_Factor

    gen_row <- data.frame(
      Component = "General_Factor",
      Mean = round(gen_factor$mean, 3),
      SD = round(gen_factor$sd, 3),
      CI_Lower = round(gen_factor$hdi_95[1], 3),
      CI_Upper = round(gen_factor$hdi_95[2], 3),
      N_Samples = gen_factor$n_samples,
      Quality = if (gen_factor$mean >= 0.8) "Good"
      else if (gen_factor$mean >= 0.7) "Acceptable"
      else if (gen_factor$mean >= 0.6) "Questionable"
      else "Poor",
      stringsAsFactors = FALSE
    )

    results_table <- rbind(results_table, gen_row)
  }

  # Export if file name is provided
  if (!is.null(file_name)) {
    write.csv(results_table, file_name, row.names = FALSE)
    cat("Results exported to:", file_name, "\n")
  }

  return(results_table)
}
