#' The Multicomponent Latent Trait Model for Diagnosis
#'
#' @description
#' Estimate the parameters of the MLTM-D via Bayesian Hamiltonian Monte Carlo.
#'
#' @param data An \eqn{n \times p} matrix or data.frame of binary responses
#'   (rows = subjects, columns = items).
#' @param Q A \eqn{p \times K} matrix specifying which cognitive rules each item
#'   requires (Q-matrix).
#' @param components A named list grouping rules into components. Each element is
#'   a numeric vector of rule indices belonging to that component.
#'   Example: \code{list(global = c(1,2,3), local = c(4,5))}.
#' @param iters Number of post-warmup MCMC iterations per chain. Default is 2000.
#' @param chains Number of Markov chains. Default is 2.
#' @param iter_warmup Number of warmup iterations per chain. Default is 1000.
#' @param quantiles Numeric vector of probabilities for posterior quantiles.
#'   Default is \code{c(0.025, 0.50, 0.975)}.
#' @param cores Number of CPU cores for parallel chains.
#'   Default is \code{parallel::detectCores() - 1}.
#' @param priors A named list of prior hyperparameters with elements \code{theta},
#'   \code{eta}, and \code{alpha}. Each is a list with \code{mu} and \code{sigma}.
#'   \code{alpha} uses a half-Normal prior (truncated at 0).
#'   Unspecified elements retain defaults.
#' @param ... Additional arguments passed to \code{rstan::sampling()}.
#'
#' @details
#' \code{MLTM} estimates the Bayesian version of the Multicomponent Latent Trait
#' Model for Diagnosis (MLTM-D; Embretson & Yang, 2013). This noncompensatory model
#' specifies a hierarchical relationship between components and rules.
#'
#' \strong{Prior distributions:} Ability (\eqn{\theta}) and rule difficulty (\eqn{\eta})
#' receive Normal priors. Discrimination (\eqn{\alpha}) receives a half-Normal prior.
#'
#' @return A list of class \code{"MLTM"} with elements:
#' \describe{
#'   \item{\code{EAP}}{Posterior mean estimates: \code{theta}, \code{alpha}, \code{eta}, \code{beta}.}
#'   \item{\code{quantiles}}{Posterior credible intervals for each parameter.}
#'   \item{\code{posterior}}{Full posterior samples and derived quantities.}
#'   \item{\code{fit}}{The \code{stanfit} object from \code{rstan::sampling}.}
#'   \item{\code{data}}{The original data matrix.}
#'   \item{\code{priors}}{The prior hyperparameters used.}
#' }
#'
#' @references
#' Embretson, S. E., & Yang, X. (2013). A multicomponent latent trait model
#' for diagnosis. \emph{Psychometrika}, \bold{78}, 14--36.
#' \doi{10.1007/s11336-012-9296-y}
#'
#' Embretson, S. E. (2019). Diagnostic modeling of skill hierarchies and
#' cognitive processes with MLTM-D. In M. von Davier & Y.-S. Lee (Eds.),
#' \emph{Handbook of Diagnostic Classification Models} (pp. 185--208).
#' Springer. \doi{10.1007/978-3-030-05584-4_9}
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
#'   components <- list(global = c(1, 2, 3), local = c(4, 5))
#'   fit <- MLTM(analogy, Q, components, iters = 200, iter_warmup = 100, chains = 1, cores = 1)
#'   fit$EAP$eta
#'   reliability(fit)
#' }
#'
#' @family model fitting functions
#' @export
MLTM <- function(data, Q, components,
                 iters = 2000, chains = 2,
                 iter_warmup = 1000,
                 quantiles = c(0.025, 0.50, 0.975),
                 cores = parallel::detectCores() - 1,
                 priors = list(
                   theta = list(mu = 0, sigma = 1),
                   eta   = list(mu = 0, sigma = 1),
                   alpha = list(mu = 0, sigma = 1)
                 ),
                 ...) {

  default_priors <- list(
    theta = list(mu = 0, sigma = 1),
    eta   = list(mu = 0, sigma = 1),
    alpha = list(mu = 0, sigma = 1)
  )
  priors <- utils::modifyList(default_priors, priors)

  for (par in c("theta", "eta", "alpha")) {
    if (!is.numeric(priors[[par]]$sigma) || priors[[par]]$sigma <= 0)
      stop(sprintf("'priors$%s$sigma' must be a positive number.", par))
  }

  binary_eta  <- get_eta(components)
  C           <- get_C(Q, components)
  indexes_eta <- locate_eta(binary_eta)
  M      <- length(components)
  K      <- ncol(Q)
  y      <- as.vector(as.matrix(data))
  N_subj <- nrow(data)
  N_item <- ncol(data)
  ID     <- rep(1:N_subj, times = N_item)
  item   <- rep(1:N_item, each  = N_subj)
  n_eta  <- nrow(indexes_eta)
  ones   <- rep(1L, N_subj * N_item)

  data_list <- list(
    N_subj = N_subj, N_item = N_item,
    ID = ID, item = item, K = K, M = M,
    indexes_eta = indexes_eta, n_eta = n_eta,
    Q = Q, C = C, y = y, ones = ones,
    prior_theta_mu    = priors$theta$mu,
    prior_theta_sigma = priors$theta$sigma,
    prior_eta_mu      = priors$eta$mu,
    prior_eta_sigma   = priors$eta$sigma,
    prior_alpha_mu    = priors$alpha$mu,
    prior_alpha_sigma = priors$alpha$sigma
  )

  stan_file <- system.file("MLTM.stan", package = "GMLTM")
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
  alpha <- as.matrix(fit, pars = "alpha")
  eta   <- as.matrix(fit, pars = "eta")

  summary_theta <- colMeans(theta)
  summary_alpha <- colMeans(alpha)
  summary_eta   <- colMeans(eta)

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
  quant_theta    <- t(apply(theta, MARGIN = 2, quantile, probs = quantiles))
  quantiles_theta <- list()
  for (j in seq_len(J)) {
    temp <- matrix(quant_theta[, j], nrow = N_subj, ncol = M)
    rownames(temp) <- as.character(subject_names)
    colnames(temp) <- comp_names
    quantiles_theta[[j]] <- temp
  }
  names(quantiles_theta) <- as.character(quantiles)

  colnames(alpha) <- comp_names
  alpha_EAP <- c(summary_alpha)
  names(alpha_EAP) <- comp_names
  quantiles_alpha <- t(apply(alpha, MARGIN = 2, quantile, probs = quantiles))
  colnames(quantiles_alpha) <- as.character(quantiles)

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

  n_draws <- nrow(theta)
  y_vec   <- as.vector(as.matrix(data))
  p <- loglik <- array(NA, dim = c(n_draws, N_subj, N_item))
  theta_matrix <- array(theta, dim = c(n_draws, N_subj, M))
  eta_matrix   <- array(0,     dim = c(n_draws, K, M))
  beta_matrix  <- array(0,     dim = c(n_draws, N_item, M))
  comp <- rep(1L, N_subj * N_item)
  for (i in seq_len(n_draws)) {
    for (j in seq_len(n_eta)) {
      eta_matrix[i, indexes_eta[j, 1], indexes_eta[j, 2]] <- eta[i, j]
    }
    beta_matrix[i, , ] <- Q %*% eta_matrix[i, , ]
    mu <- plogis(alpha[i, ][comp] * (theta_matrix[i, ID, ] - beta_matrix[i, item, ])) ^ C[item, ]
    p[i, , ]      <- apply(mu, MARGIN = 1, prod)
    loglik[i, , ] <- dbinom(y_vec, size = 1, prob = c(p[i, , ]), log = TRUE)
  }

  EAP       <- list(theta = theta_EAP, alpha = alpha_EAP, eta = eta_EAP, beta = beta_EAP)
  quantiles_out <- list(theta = quantiles_theta, alpha = quantiles_alpha,
                        eta = quantiles_eta, beta = quantiles_beta)
  posterior <- list(theta = theta_matrix, alpha = alpha, eta = eta_matrix,
                    beta = beta_matrix, probabilities = p, loglik = loglik)

  result <- list(EAP = EAP, quantiles = quantiles_out, posterior = posterior,
                 fit = fit, data = data, priors = priors)
  class(result) <- "MLTM"
  return(result)
}

#' Print method for MLTM objects
#'
#' @param x An object of class \code{"MLTM"}, as returned by \code{\link{MLTM}}.
#' @param digits Number of decimal digits to display. Default is 3.
#' @param ... Currently ignored.
#'
#' @return Invisibly returns \code{x}.
#' @export
print.MLTM <- function(x, digits = 3, ...) {
  cat("MLTM-D model fit\n")
  cat("----------------\n")
  cat(sprintf("Subjects: %d | Components: %d\n",
              nrow(x$EAP$theta), ncol(x$EAP$theta)))
  cat(sprintf("Components: %s\n\n", paste(colnames(x$EAP$theta), collapse = ", ")))
  cat("Posterior mean (EAP) rule difficulty (eta):\n")
  print(round(x$EAP$eta, digits))
  cat(sprintf("\nDiscrimination (alpha, EAP): %s\n",
              paste(round(x$EAP$alpha, digits), collapse = ", ")))
  cat("\nUse summary(x$fit) or x$quantiles for credible intervals.\n")
  invisible(x)
}
