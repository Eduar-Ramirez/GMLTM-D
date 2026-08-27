get_eta <- function(components) {
  M   <- length(components)
  K   <- length(unique(unlist(components)))
  eta <- matrix(0, K, M)
  for (m in seq_len(M)) eta[components[[m]], m] <- 1
  return(eta)
}

get_C <- function(Q, components) {
  M <- length(components)
  C <- matrix(0, nrow(Q), M)
  for (i in seq_len(M)) {
    sums1   <- rowSums(Q[, components[[i]], drop = FALSE])
    indexes <- which(sums1 > 0)
    C[indexes, i] <- 1
  }
  return(C)
}

locate_eta <- function(eta) {
  indexes <- which(ifelse(eta == 1, TRUE, FALSE), arr.ind = TRUE)
  return(indexes)
}

locate_alpha <- function(Q, components) {
  C        <- get_C(Q, components)
  indexes_2 <- locate_eta(C)
  M <- length(components)
  i <- 1
  x <- c()

  for (m in seq_len(M)) {
    X  <- Q[, components[[m]], drop = FALSE]
    D  <- cbind(C[, m], 0)
    Z  <- X %*% t(X) - rowSums(X)
    Z2 <- (t(Z) + Z) / 2
    nr <- seq_len(nrow(X))

    repeat {
      j       <- min(nr)
      indexes <- which(Z2[j, ] == 0)
      nr      <- nr[-match(indexes, nr)]
      D[indexes, 2] <- i
      if (length(nr) == 0) break
      i <- i + 1
    }
    x <- c(x, D[D[, 1] != 0, 2])
  }

  x         <- as.numeric(as.factor(x))
  indexes_2 <- cbind(indexes_2, x)
  return(indexes_2)
}

#' The Generalized Multicomponent Latent Trait Model for Diagnosis
#'
#' @description
#' Estimate the parameters of the GMLTM-D via Bayesian Hamiltonian Monte Carlo.
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
#'   \code{eta}, \code{alpha}, and \code{c}. For Normal parameters supply
#'   \code{mu} and \code{sigma}; for the guessing parameter supply \code{shape1}
#'   and \code{shape2} (Beta prior). \code{priors$alpha} additionally accepts
#'   \code{family}, either \code{"normal"} (default; a half-Normal, i.e. Normal
#'   truncated at 0) or \code{"lognormal"} (a Log-Normal prior on \eqn{\alpha},
#'   using the same \code{mu}/\code{sigma} as the location/scale of the
#'   underlying Normal on the log scale). Unspecified elements retain defaults.
#' @param ... Additional arguments passed to \code{rstan::sampling()}.
#'
#' @details
#' \code{GMLTM} estimates the Generalized Multicomponent Latent Trait Model for
#' Diagnosis (GMLTM-D; Ramirez et al., 2024) in its Bayesian version. This model
#' analyses items composed of cognitive rules or operations, incorporating three IRT
#' parameters. Rules can be grouped into distinct components.
#'
#' \strong{Prior distributions:} Ability (\eqn{\theta}) and rule difficulty
#' (\eqn{\eta}) receive Normal priors. Discrimination (\eqn{\alpha}) receives
#' either a half-Normal prior (\code{priors$alpha$family = "normal"}, default)
#' or a Log-Normal prior (\code{priors$alpha$family = "lognormal"}). Guessing
#' (\eqn{c}) receives a Beta prior.
#'
#' @return A list of class \code{"GMLTM"} with elements:
#' \describe{
#'   \item{\code{EAP}}{Posterior mean estimates: \code{theta}, \code{alpha},
#'     \code{eta}, \code{beta}, \code{guessing}.}
#'   \item{\code{quantiles}}{Posterior credible intervals for each parameter.}
#'   \item{\code{posterior}}{Full posterior samples and derived quantities.}
#'   \item{\code{fit}}{The \code{stanfit} object from \code{rstan::sampling}.}
#'   \item{\code{data}}{The original data matrix.}
#'   \item{\code{priors}}{The prior hyperparameters used.}
#' }
#'
#' @references
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
#'   fit <- GMLTM(analogy, Q, components, iters = 200, iter_warmup = 100, chains = 1, cores = 1)
#'   fit$EAP$eta
#'   reliability(fit)
#' }
#'
#' @family model fitting functions
#' @export
GMLTM <- function(data, Q, components,
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
    alpha_lognormal   = as.integer(priors$alpha$family == "lognormal")
  )

  stan_file <- system.file("GMLTM.stan", package = "GMLTM")
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

  summary_theta    <- colMeans(theta)
  summary_alpha    <- colMeans(alpha)
  summary_eta      <- colMeans(eta)
  summary_guessing <- colMeans(guessing)

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
              beta = beta_EAP, guessing = guessing_EAP)
  quantiles_out <- list(theta = quantiles_theta, alpha = quantiles_alpha,
                        eta = quantiles_eta, beta = quantiles_beta,
                        guessing = quantiles_guessing)
  posterior <- list(theta = theta_matrix, alpha = alpha_matrix, eta = eta_matrix,
                    beta = beta_matrix, guessing = guessing,
                    probabilities = p, loglik = loglik)

  result <- list(EAP = EAP, quantiles = quantiles_out, posterior = posterior,
                 fit = fit, data = data, priors = priors)
  class(result) <- "GMLTM"
  return(result)
}

#' Print method for GMLTM objects
#'
#' @param x An object of class \code{"GMLTM"}, as returned by \code{\link{GMLTM}}.
#' @param digits Number of decimal digits to display. Default is 3.
#' @param ... Currently ignored.
#'
#' @return Invisibly returns \code{x}.
#' @export
print.GMLTM <- function(x, digits = 3, ...) {
  cat("GMLTM-D model fit\n")
  cat("-----------------\n")
  cat(sprintf("Subjects: %d | Items: %d | Components: %d\n",
              nrow(x$EAP$theta), nrow(x$EAP$alpha), ncol(x$EAP$theta)))
  cat(sprintf("Components: %s\n\n", paste(colnames(x$EAP$theta), collapse = ", ")))
  cat("Posterior mean (EAP) rule difficulty (eta):\n")
  print(round(x$EAP$eta, digits))
  cat("\nUse summary(x$fit) or x$quantiles for credible intervals.\n")
  invisible(x)
}
