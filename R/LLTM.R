#' The Linear Logistic Test Model
#'
#' @description
#' Estimate the parameters of the LLTM via Bayesian Hamiltonian Monte Carlo.
#'
#' @param data An \eqn{n \times p} matrix or data.frame of binary responses
#'   (rows = subjects, columns = items).
#' @param Q A \eqn{p \times K} matrix specifying which cognitive rules each item
#'   requires (Q-matrix).
#' @param iters Number of post-warmup MCMC iterations per chain. Default is 2000.
#' @param chains Number of Markov chains. Default is 2.
#' @param iter_warmup Number of warmup iterations per chain. Default is 1000.
#' @param quantiles Numeric vector of probabilities for posterior quantiles.
#'   Default is \code{c(0.025, 0.50, 0.975)}.
#' @param cores Number of CPU cores for parallel chains.
#'   Default is \code{parallel::detectCores() - 1}.
#' @param priors A named list of prior hyperparameters. Each element is a named
#'   list with \code{mu} and \code{sigma}. Available parameters: \code{theta}
#'   (ability, Normal prior) and \code{eta} (rule difficulty, Normal prior).
#'   Unspecified elements retain defaults. Example:
#'   \code{priors = list(eta = list(sigma = 3))}.
#' @param ... Additional arguments passed to \code{rstan::sampling()}.
#'
#' @details
#' \code{LLTM} estimates the Bayesian version of the Linear Logistic Test Model
#' (Fischer, 1973), which extends the Rasch model by decomposing item difficulty
#' into cognitive rules. Item difficulty is expressed as
#' \eqn{\beta_i = \mathbf{q}_i^\top \boldsymbol{\eta}},
#' where \eqn{\mathbf{q}_i} is the i-th row of Q and \eqn{\boldsymbol{\eta}} is the
#' vector of rule difficulty parameters.
#'
#' \strong{Prior distributions:} Ability (\eqn{\theta}) and rule difficulty
#' (\eqn{\eta}) receive Normal priors. Prior sensitivity analysis is recommended.
#'
#' @return A list of class \code{"LLTM"} with elements:
#' \describe{
#'   \item{\code{EAP}}{Posterior mean estimates: \code{theta}, \code{eta}, \code{beta}.}
#'   \item{\code{quantiles}}{Posterior credible intervals for each parameter.}
#'   \item{\code{posterior}}{Full posterior samples and derived quantities.}
#'   \item{\code{fit}}{The \code{stanfit} object from \code{rstan::sampling}.}
#'   \item{\code{data}}{The original data matrix.}
#'   \item{\code{priors}}{The prior hyperparameters used.}
#' }
#'
#' @references
#' Fischer, G. H. (1973). The linear logistic test model as an instrument
#' in educational research. \emph{Acta Psychologica}, \bold{37}(6), 359--374.
#' \doi{10.1016/0001-6918(73)90003-6}
#'
#' Ramirez, E. S., Jimenez, M., Franco, V. R., & Alvarado, J. M. (2024).
#' Delving into the complexity of analogical reasoning: A detailed exploration
#' with the Generalized Multicomponent Latent Trait Model for Diagnosis.
#' \emph{Journal of Intelligence}, \bold{12}, 67.
#' \doi{10.3390/jintelligence12070067}
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
#'   fit$EAP$eta
#'   reliability(fit)
#' }
#'
#' @family model fitting functions
#' @export
LLTM <- function(data, Q,
                 iters = 2000, chains = 2,
                 iter_warmup = 1000,
                 quantiles = c(0.025, 0.50, 0.975),
                 cores = parallel::detectCores() - 1,
                 priors = list(
                   theta = list(mu = 0, sigma = 1),
                   eta   = list(mu = 0, sigma = 1)
                 ),
                 ...) {

  default_priors <- list(
    theta = list(mu = 0, sigma = 1),
    eta   = list(mu = 0, sigma = 1)
  )
  priors <- utils::modifyList(default_priors, priors)

  if (!is.numeric(priors$theta$sigma) || priors$theta$sigma <= 0)
    stop("'priors$theta$sigma' must be a positive number.")
  if (!is.numeric(priors$eta$sigma) || priors$eta$sigma <= 0)
    stop("'priors$eta$sigma' must be a positive number.")

  K      <- ncol(Q)
  y      <- as.vector(as.matrix(data))
  N_subj <- nrow(data)
  N_item <- ncol(data)
  ID     <- rep(1:N_subj, times = N_item)
  item   <- rep(1:N_item, each  = N_subj)

  data_list <- list(
    N_subj = N_subj, N_item = N_item,
    ID = ID, item = item, K = K, Q = Q, y = y,
    prior_theta_mu    = priors$theta$mu,
    prior_theta_sigma = priors$theta$sigma,
    prior_eta_mu      = priors$eta$mu,
    prior_eta_sigma   = priors$eta$sigma
  )

  stan_file <- system.file("LLTM.stan", package = "GMLTM")
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

  theta <- as.matrix(fit, pars = "theta")
  eta   <- as.matrix(fit, pars = "eta")
  beta  <- as.matrix(fit, pars = "beta")

  summary_theta <- colMeans(theta)
  summary_eta   <- colMeans(eta)
  summary_beta  <- colMeans(beta)

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
  J <- length(quantiles)

  colnames(theta) <- as.character(subject_names)
  theta_EAP <- summary_theta
  names(theta_EAP) <- as.character(subject_names)
  quantiles_theta <- t(apply(theta, MARGIN = 2, quantile, probs = quantiles))
  rownames(quantiles_theta) <- as.character(subject_names)
  colnames(quantiles_theta) <- as.character(quantiles)

  colnames(eta)  <- rule_names
  colnames(beta) <- item_names
  eta_EAP  <- summary_eta
  names(eta_EAP) <- rule_names
  beta_EAP <- Q %*% eta_EAP
  quantiles_eta  <- t(apply(eta,  MARGIN = 2, quantile, probs = quantiles))
  quantiles_beta <- t(apply(beta, MARGIN = 2, quantile, probs = quantiles))
  rownames(quantiles_eta)  <- rule_names
  colnames(quantiles_eta)  <- as.character(quantiles)
  rownames(quantiles_beta) <- item_names
  colnames(quantiles_beta) <- as.character(quantiles)

  n_draws <- nrow(theta)
  y_vec   <- as.vector(as.matrix(data))
  p <- loglik <- array(NA, dim = c(n_draws, N_subj, N_item))
  for (i in seq_len(n_draws)) {
    p[i, , ]      <- plogis(theta[i, ID] - beta[i, item])
    loglik[i, , ] <- dbinom(y_vec, size = 1, prob = c(p[i, , ]), log = TRUE)
  }

  EAP      <- list(theta = theta_EAP, eta = eta_EAP, beta = beta_EAP)
  quantiles_out <- list(theta = quantiles_theta, eta = quantiles_eta, beta = quantiles_beta)
  posterior <- list(theta = theta, eta = eta, beta = beta,
                    probabilities = p, loglik = loglik)

  result <- list(EAP = EAP, quantiles = quantiles_out, posterior = posterior,
                 fit = fit, data = data, priors = priors)
  class(result) <- "LLTM"
  return(result)
}

#' Print method for LLTM objects
#'
#' @param x An object of class \code{"LLTM"}, as returned by \code{\link{LLTM}}.
#' @param digits Number of decimal digits to display. Default is 3.
#' @param ... Currently ignored.
#'
#' @return Invisibly returns \code{x}.
#' @export
print.LLTM <- function(x, digits = 3, ...) {
  cat("LLTM model fit\n")
  cat("--------------\n")
  cat(sprintf("Subjects: %d | Items: %d | Rules: %d\n",
              length(x$EAP$theta), length(x$EAP$beta), length(x$EAP$eta)))
  cat("\nPosterior mean (EAP) rule difficulty (eta):\n")
  print(round(x$EAP$eta, digits))
  cat("\nUse summary(x$fit) or x$quantiles for credible intervals.\n")
  invisible(x)
}
