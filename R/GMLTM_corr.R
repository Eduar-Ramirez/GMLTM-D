# =====================================================================
# GMLTM_corr(): variant of the GMLTM-D with a CORRELATED multivariate
# latent structure (general Sigma, estimated via a Cholesky factor and
# an LKJ prior), instead of the Sigma = sigma^2 * I used by GMLTM().
#
# Requires the helper functions get_eta(), get_C(), locate_eta(),
# locate_alpha() from R/GMLTM.R (identical, unmodified) and the
# GMLTM_corr.stan file (same inst/ directory as GMLTM.stan).
# =====================================================================

#' GMLTM-D with correlated latent components (general Sigma)
#'
#' @description
#' Variant of \code{\link{GMLTM}} that replaces the independent prior on
#' \eqn{\theta} (\eqn{\Sigma = \sigma^2 I}) with a multivariate Normal prior
#' with a general correlation matrix \eqn{\Sigma}, estimated from the data via
#' its Cholesky factor (non-centered parameterization with an LKJ prior on the
#' correlation). This corresponds to innovation (iii) of the GMLTM-D as
#' formulated in Ramirez et al. (2024): allowing correlations between
#' components to reflect the natural interdependence between cognitive
#' abilities, instead of assuming independent latent traits.
#'
#' @inheritParams GMLTM
#' @param lkj_eta Concentration parameter of the LKJ prior on the correlation
#'   matrix \eqn{\Sigma}. \code{lkj_eta = 1} (default) is uniform over the
#'   space of correlation matrices; \code{lkj_eta > 1} concentrates the prior
#'   towards \eqn{\Sigma = I} (low correlations); \code{lkj_eta < 1} towards
#'   more extreme correlations.
#'
#' @details
#' \code{GMLTM_corr} estimates the same likelihood as \code{\link{GMLTM}}
#' (three IRT parameters per item/component: discrimination \eqn{\alpha},
#' guessing \eqn{c}, and rule difficulty \eqn{\eta}/\eqn{\beta}), but replaces
#' the prior on the subject-level ability matrix \eqn{\theta} with
#' \eqn{\theta_i \sim N_M(0, \Sigma)}, where \eqn{\Sigma} is a correlation
#' matrix among the \code{M} cognitive components. \eqn{\Sigma} is
#' parameterized through its Cholesky factor \code{L_Omega} (a
#' \code{cholesky_factor_corr}), with \code{theta = t(L_Omega \%*\% theta_raw)}
#' and \code{theta_raw ~ N(0, 1)} i.i.d. (non-centered parameterization, which
#' improves posterior geometry relative to sampling \eqn{\theta} directly from
#' the multivariate Normal). \code{priors$theta} is accepted for interface
#' compatibility with \code{\link{GMLTM}} but is not used: the scale of
#' \eqn{\theta} is fixed to 1 by construction (diagonal of \eqn{\Sigma}), and
#' its prior is governed entirely by \code{lkj_eta}.
#'
#' @return A list of class \code{c("GMLTM_corr", "GMLTM")} with elements:
#' \describe{
#'   \item{\code{EAP}}{Posterior mean estimates: \code{theta}, \code{alpha},
#'     \code{eta}, \code{beta}, \code{guessing}, and \code{Sigma} (the
#'     estimated \code{M x M} correlation matrix among components).}
#'   \item{\code{quantiles}}{Posterior credible intervals for each parameter,
#'     including \code{Sigma}.}
#'   \item{\code{posterior}}{Full posterior samples and derived quantities
#'     (including \code{loglik}, used by \code{\link{compute_model_validation}}).}
#'   \item{\code{fit}}{The \code{stanfit} object from \code{rstan::sampling}.}
#'   \item{\code{data}}{The original data matrix.}
#'   \item{\code{priors}}{The prior hyperparameters used.}
#'   \item{\code{lkj_eta}}{The LKJ concentration parameter used.}
#' }
#'
#' @references
#' Ramirez, E. S., Jimenez, M., Franco, V. R., & Alvarado, J. M. (2024).
#' Delving into the complexity of analogical reasoning: A detailed exploration
#' with the Generalized Multicomponent Latent Trait Model for Diagnosis.
#' \emph{Journal of Intelligence}, \bold{12}, 67.
#' \doi{10.3390/jintelligence12070067}
#'
#' Embretson, S. E., & Yang, X. (2013). A multicomponent latent trait model
#' for diagnosis. \emph{Psychometrika}, \bold{78}, 14--36.
#' \doi{10.1007/s11336-012-9296-y}
#'
#' Lewandowski, D., Kurowicka, D., & Joe, H. (2009). Generating random
#' correlation matrices based on vines and extended onion method.
#' \emph{Journal of Multivariate Analysis}, 100(9), 1989--2001.
#' \doi{10.1016/j.jmva.2009.04.008}
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
#'   fit <- GMLTM_corr(analogy, Q, components,
#'                      iters = 200, iter_warmup = 100, chains = 1, cores = 1)
#'   fit$EAP$eta
#'   fit$EAP$Sigma
#' }
#'
#' @family model fitting functions
#' @export
GMLTM_corr <- function(data, Q, components,
                        iters = 2000, chains = 2,
                        iter_warmup = 1000,
                        quantiles = c(0.025, 0.50, 0.975),
                        cores = parallel::detectCores() - 1,
                        priors = list(
                          theta = list(mu = 0, sigma = 1),
                          eta   = list(mu = 0, sigma = 1),
                          alpha = list(mu = 0, sigma = 1, family = "normal"),
                          c     = list(shape1 = 3, shape2 = 20)
                        ),
                        lkj_eta = 1,
                        ...) {

  default_priors <- list(
    theta = list(mu = 0, sigma = 1),
    eta   = list(mu = 0, sigma = 1),
    alpha = list(mu = 0, sigma = 1, family = "normal"),
    c     = list(shape1 = 3, shape2 = 20)
  )
  priors <- utils::modifyList(default_priors, priors)

  for (par in c("theta", "eta", "alpha")) {
    if (!is.numeric(priors[[par]]$sigma) || priors[[par]]$sigma <= 0)
      stop(sprintf("'priors$%s$sigma' must be a positive number.", par))
  }
  if (!is.numeric(priors$c$shape1) || priors$c$shape1 <= 0)
    stop("'priors$c$shape1' must be a positive number.")
  if (!is.numeric(priors$c$shape2) || priors$c$shape2 <= 0)
    stop("'priors$c$shape2' must be a positive number.")
  if (!is.character(priors$alpha$family) || length(priors$alpha$family) != 1 ||
      !priors$alpha$family %in% c("normal", "lognormal"))
    stop("'priors$alpha$family' must be either \"normal\" or \"lognormal\".")
  if (!is.numeric(lkj_eta) || length(lkj_eta) != 1 || is.na(lkj_eta) || lkj_eta <= 0)
    stop("'lkj_eta' must be a positive number.")

  binary_eta    <- get_eta(components)
  C             <- get_C(Q, components)
  indexes_eta   <- locate_eta(binary_eta)
  indexes_alpha <- locate_alpha(Q, components)
  M      <- length(components)
  K      <- ncol(Q)
  y      <- as.vector(as.matrix(data))
  N_subj <- nrow(data)
  N_item <- ncol(data)
  ID     <- rep(1:N_subj, times = N_item)
  item   <- rep(1:N_item, each  = N_subj)
  n_alpha <- nrow(indexes_alpha)
  n_eta   <- nrow(indexes_eta)

  data_list <- list(
    N_subj = N_subj, N_item = N_item,
    ID = ID, item = item, K = K, M = M,
    indexes_eta   = indexes_eta,   n_eta   = n_eta,
    indexes_alpha = indexes_alpha, n_alpha = n_alpha,
    Q = Q, C = C, y = y, cf = 1L,
    prior_theta_mu    = priors$theta$mu,
    prior_theta_sigma = priors$theta$sigma,
    prior_eta_mu      = priors$eta$mu,
    prior_eta_sigma   = priors$eta$sigma,
    prior_alpha_mu    = priors$alpha$mu,
    prior_alpha_sigma = priors$alpha$sigma,
    prior_c_shape1    = priors$c$shape1,
    prior_c_shape2    = priors$c$shape2,
    prior_theta_lkj_eta = lkj_eta,
    alpha_lognormal   = as.integer(priors$alpha$family == "lognormal")
  )

  stan_file <- system.file("GMLTM_corr.stan", package = "GMLTM")
  model <- rstan::stan_model(stan_file)

  fit <- rstan::sampling(
    model,
    data   = data_list,
    chains = chains,
    iter   = iters + iter_warmup,
    warmup = iter_warmup,
    cores  = cores,
    ...
  )

  theta    <- as.matrix(fit, pars = "theta")
  alpha    <- as.matrix(fit, pars = "alpha")
  eta      <- as.matrix(fit, pars = "eta")
  guessing <- as.matrix(fit, pars = "c")
  Sigma    <- as.matrix(fit, pars = "Sigma")

  summary_theta    <- colMeans(theta)
  summary_alpha    <- colMeans(alpha)
  summary_eta      <- colMeans(eta)
  summary_guessing <- colMeans(guessing)
  summary_Sigma    <- colMeans(Sigma)

  if (is.null(colnames(data))) {
    item_names <- paste0("item", seq_len(N_item))
  } else {
    item_names <- colnames(data)
  }
  if (is.null(rownames(data))) {
    subject_names <- seq_len(N_subj)
  } else {
    subject_names <- rownames(data)
  }
  if (is.null(colnames(Q))) {
    rule_names <- paste0("rule", seq_len(K))
  } else {
    rule_names <- colnames(Q)
  }
  rownames(Q) <- item_names
  comp_names  <- names(components)
  if (is.null(comp_names)) comp_names <- paste0("Component", seq_len(M))
  J <- length(quantiles)

  theta_EAP <- matrix(summary_theta, nrow = N_subj, ncol = M)
  rownames(theta_EAP) <- as.character(subject_names)
  colnames(theta_EAP) <- comp_names
  quant_theta     <- t(apply(theta, MARGIN = 2, quantile, probs = quantiles))
  quantiles_theta <- list()
  for (j in seq_len(J)) {
    temp <- matrix(quant_theta[, j], nrow = N_subj, ncol = M)
    rownames(temp) <- as.character(subject_names)
    colnames(temp) <- comp_names
    quantiles_theta[[j]] <- temp
  }
  names(quantiles_theta) <- as.character(quantiles)

  alpha_EAP <- matrix(0, N_item, M)
  for (i in seq_len(n_alpha)) {
    alpha_EAP[indexes_alpha[i, 1], indexes_alpha[i, 2]] <- summary_alpha[indexes_alpha[i, 3]]
  }
  rownames(alpha_EAP) <- item_names
  colnames(alpha_EAP) <- comp_names
  quant_alpha     <- t(apply(alpha, MARGIN = 2, quantile, probs = quantiles))
  quantiles_alpha <- list()
  for (j in seq_len(J)) {
    temp <- matrix(0, nrow = N_item, ncol = M)
    for (i in seq_len(n_alpha)) {
      temp[indexes_alpha[i, 1], indexes_alpha[i, 2]] <- quant_alpha[indexes_alpha[i, 3], j]
    }
    rownames(temp) <- item_names
    colnames(temp) <- comp_names
    quantiles_alpha[[j]] <- temp
  }
  names(quantiles_alpha) <- as.character(quantiles)

  eta_EAP <- matrix(0, K, M)
  for (i in seq_len(n_eta)) {
    eta_EAP[indexes_eta[i, 1], indexes_eta[i, 2]] <- summary_eta[i]
  }
  rownames(eta_EAP) <- rule_names
  colnames(eta_EAP) <- comp_names
  beta_EAP <- Q %*% eta_EAP
  quant_eta <- t(apply(eta, MARGIN = 2, quantile, probs = quantiles))
  quantiles_eta <- quantiles_beta <- list()
  for (j in seq_len(J)) {
    temp <- matrix(0, nrow = K, ncol = M)
    for (i in seq_len(n_eta)) {
      temp[indexes_eta[i, 1], indexes_eta[i, 2]] <- quant_eta[i, j]
    }
    rownames(temp) <- rule_names
    colnames(temp) <- comp_names
    quantiles_eta[[j]]  <- temp
    quantiles_beta[[j]] <- Q %*% temp
  }
  names(quantiles_eta)  <- as.character(quantiles)
  names(quantiles_beta) <- as.character(quantiles)

  guessing_EAP <- summary_guessing
  names(guessing_EAP) <- item_names
  quant_guessing     <- t(apply(guessing, MARGIN = 2, quantile, probs = quantiles))
  quantiles_guessing <- list()
  for (j in seq_len(J)) {
    temp <- quant_guessing[, j]
    names(temp) <- item_names
    quantiles_guessing[[j]] <- temp
  }
  names(quantiles_guessing) <- as.character(quantiles)

  Sigma_EAP <- matrix(summary_Sigma, M, M)
  rownames(Sigma_EAP) <- colnames(Sigma_EAP) <- comp_names
  quant_Sigma     <- t(apply(Sigma, MARGIN = 2, quantile, probs = quantiles))
  quantiles_Sigma <- list()
  for (j in seq_len(J)) {
    temp <- matrix(quant_Sigma[, j], M, M)
    rownames(temp) <- colnames(temp) <- comp_names
    quantiles_Sigma[[j]] <- temp
  }
  names(quantiles_Sigma) <- as.character(quantiles)

  n_draws <- nrow(theta)
  y_vec   <- as.vector(as.matrix(data))
  p <- loglik <- array(NA, dim = c(n_draws, N_subj, N_item))
  theta_matrix <- array(theta, dim = c(n_draws, N_subj, M))
  alpha_matrix <- array(0,     dim = c(n_draws, N_item, M))
  eta_matrix   <- array(0,     dim = c(n_draws, K, M))
  beta_matrix  <- array(0,     dim = c(n_draws, N_item, M))

  for (i in seq_len(n_draws)) {
    for (j in seq_len(n_alpha)) {
      alpha_matrix[i, indexes_alpha[j, 1], indexes_alpha[j, 2]] <- alpha[i, indexes_alpha[j, 3]]
    }
    for (j in seq_len(n_eta)) {
      eta_matrix[i, indexes_eta[j, 1], indexes_eta[j, 2]] <- eta[i, j]
    }
    beta_matrix[i, , ] <- Q %*% eta_matrix[i, , ]
    mu <- plogis(alpha_matrix[i, item, ] * (theta_matrix[i, ID, ] - beta_matrix[i, item, ])) ^ C[item, ]
    p[i, , ]      <- guessing[i, item] + (1 - guessing[i, item]) * apply(mu, MARGIN = 1, prod)
    loglik[i, , ] <- dbinom(y_vec, size = 1, prob = c(p[i, , ]), log = TRUE)
  }

  EAP <- list(theta = theta_EAP, alpha = alpha_EAP, eta = eta_EAP,
              beta = beta_EAP, guessing = guessing_EAP, Sigma = Sigma_EAP)
  quantiles_out <- list(theta = quantiles_theta, alpha = quantiles_alpha,
                        eta = quantiles_eta, beta = quantiles_beta,
                        guessing = quantiles_guessing, Sigma = quantiles_Sigma)
  posterior <- list(theta = theta_matrix, alpha = alpha_matrix, eta = eta_matrix,
                    beta = beta_matrix, guessing = guessing, Sigma = Sigma,
                    probabilities = p, loglik = loglik)

  result <- list(EAP = EAP, quantiles = quantiles_out, posterior = posterior,
                 fit = fit, data = data, priors = priors, lkj_eta = lkj_eta)
  class(result) <- c("GMLTM_corr", "GMLTM")
  return(result)
}

#' Print method for GMLTM_corr objects
#'
#' @param x An object of class \code{c("GMLTM_corr", "GMLTM")}, as returned by
#'   \code{\link{GMLTM_corr}}.
#' @param digits Number of decimal digits to display. Default is 3.
#' @param ... Currently ignored.
#'
#' @return Invisibly returns \code{x}.
#' @export
print.GMLTM_corr <- function(x, digits = 3, ...) {
  cat("GMLTM-D model fit (correlated components)\n")
  cat("------------------------------------------\n")
  cat(sprintf("Subjects: %d | Items: %d | Components: %d\n",
              nrow(x$EAP$theta), nrow(x$EAP$alpha), ncol(x$EAP$theta)))
  cat("\nPosterior mean correlation matrix (Sigma):\n")
  print(round(x$EAP$Sigma, digits))
  cat("\nPosterior mean (EAP) rule difficulty (eta):\n")
  print(round(x$EAP$eta, digits))
  cat("\nUse extract_correlation(x) for pairwise credible intervals,\n")
  cat("or summary(x$fit)/x$quantiles for full posterior summaries.\n")
  invisible(x)
}
