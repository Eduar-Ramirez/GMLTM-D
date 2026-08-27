#' Compute LOO and WAIC for GMLTM models
#'
#' This function extracts the log-likelihood from a GMLTM model
#' and computes the Leave-One-Out Cross-Validation (LOO) and
#' the Widely Applicable Information Criterion (WAIC).
#' LOO is a Bayesian model comparison metric based on Pareto-smoothed importance
#' sampling, while WAIC is a fully Bayesian criterion that estimates predictive
#' accuracy.
#'
#' @param fit A fitted GMLTM model or a list of fitted models. Models fitted with
#' \code{\link{GMLTM_corr}} (correlated \eqn{\Sigma}), \code{\link{GMLTM}}
#' (\eqn{\Sigma = \sigma^2 I}), \code{\link{MLTM}}, and \code{\link{LLTM}} can be
#' freely mixed in the same list; each is labeled by model class in the output.
#' @return If a single model is provided, returns a list with LOO and WAIC results.
#' If multiple models are provided, returns a summary table with key LOO and WAIC
#' indices, labeled by model type (e.g. \code{"GMLTM-D (Sigma libre)"} for
#' \code{GMLTM_corr} fits, \code{"GMLTM-D (Sigma=I)"} for \code{GMLTM} fits).
#' @references
#' Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian model
#' evaluation using LOO-CV and WAIC. Statistics and Computing, 27(5), 1413--1432.
#' \doi{10.1007/s11222-016-9696-4}
#' @importFrom loo loo waic loo_compare
#' @family model validation functions
#' @export
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
#'   fit1 <- GMLTM(data = analogy, Q = Q, components = components,
#'                 iters = 200, iter_warmup = 100, chains = 1, cores = 1)
#'   compute_model_validation(fit1)
#' }
compute_model_validation <- function(fit) {
  model_label <- function(m, i) {
    if (inherits(m, "GMLTM_corr")) return("GMLTM-D (Sigma libre)")
    if (inherits(m, "GMLTM"))      return("GMLTM-D (Sigma=I)")
    if (inherits(m, "MLTM"))       return("MLTM-D")
    if (inherits(m, "LLTM"))       return("LLTM")
    paste0("Model", i)
  }

  if (is.list(fit) && all(sapply(fit, function(m) "posterior" %in% names(m) && "loglik" %in% names(m$posterior)))) {
    log_lik_list <- lapply(fit, function(m) {
      dims <- dim(m$posterior$loglik)
      matrix(m$posterior$loglik, nrow = dims[1], ncol = prod(dims[-1]))
    })

    labels <- mapply(model_label, fit, seq_along(fit))
    labels <- make.unique(unname(labels))

    loo_results  <- lapply(log_lik_list, loo::loo)
    waic_results <- lapply(log_lik_list, loo::waic)
    names(loo_results)  <- labels
    names(waic_results) <- labels

    loo_summary <- do.call(rbind, lapply(loo_results, function(x) data.frame(
      elpd_loo = x$estimates["elpd_loo", "Estimate"],
      p_loo    = x$estimates["p_loo",    "Estimate"],
      looic    = x$estimates["looic",    "Estimate"]
    )))

    waic_summary <- do.call(rbind, lapply(waic_results, function(x) data.frame(
      elpd_waic = x$estimates["elpd_waic", "Estimate"],
      p_waic    = x$estimates["p_waic",    "Estimate"],
      waic      = x$estimates["waic",      "Estimate"]
    )))

    summary_table <- cbind(loo_summary, waic_summary[, -1])
    rownames(summary_table) <- labels
    comparison <- loo::loo_compare(loo_results)

    return(list(LOO = loo_results, WAIC = waic_results,
                Summary = summary_table, Comparison = comparison))
  }

  if (!"posterior" %in% names(fit) || !"loglik" %in% names(fit$posterior)) {
    stop("The model does not contain log-likelihood in fit$posterior$loglik.")
  }

  dims <- dim(fit$posterior$loglik)
  log_lik_matrix <- matrix(fit$posterior$loglik, nrow = dims[1], ncol = prod(dims[-1]))

  loo_result  <- loo::loo(log_lik_matrix)
  waic_result <- loo::waic(log_lik_matrix)

  return(list(LOO = loo_result, WAIC = waic_result))
}
