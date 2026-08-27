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
#'   fit <- LLTM(analogy, Q, iters = 200, iter_warmup = 100, chains = 1, cores = 1)
#'   quick_rel <- quick_reliability_check(fit)
#'   print(quick_rel)
#' }
#'
#' @family reliability functions
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
#' @return Invisibly returns a list with the reliability estimates computed
#'   at each step of the analysis. Called primarily for its side effect of
#'   printing a step-by-step explanation to the console.
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
#'   fit <- GMLTM(analogy, Q, components, iters = 200, iter_warmup = 100, chains = 1, cores = 1)
#'   demo_reliability_analysis(fit)
#' }
#'
#' @family reliability functions
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
#'   dq <- check_reliability_data_quality(fit)
#'   print(dq)
#' }
#'
#' @family reliability functions
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
#' @return Invisibly returns \code{x}. Called for its side effect of printing
#'   the data quality diagnostics to the console.
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
#' @return Invisibly returns \code{NULL}. Called for its side effect of
#'   printing step-by-step usage instructions to the console.
#' @examples
#' reliability_usage_instructions()
#'
#' @family reliability functions
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
#'   rel <- enhanced_mltm_reliability(fit, n_samples = 100)
#'   export_reliability_results(rel)
#' }
#'
#' @family reliability functions
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
    message("Results exported to: ", file_name)
  }

  return(results_table)
}
