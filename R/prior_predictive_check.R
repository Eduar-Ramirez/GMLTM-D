# =====================================================================
# prior_predictive_check(): simulate parameters from the priors declared
# for GMLTM()/GMLTM_corr() and simulate response data from them, to
# inspect whether the resulting data are substantively plausible before
# fitting the model (Gelman et al., 2020, "Bayesian Workflow", Sect. 2.4).
#
# Pure-R simulation, mirroring the exact generative structure used
# internally by GMLTM() (R/GMLTM.R) and GMLTM.stan: eta ~ Normal per
# rule/component, beta = Q %*% eta, theta ~ Normal per subject/component,
# c ~ Beta per item, and alpha drawn from the SAME family/hyperparameters
# a user would pass to GMLTM() (half-Normal or Log-Normal), respecting the
# item-pattern grouping of locate_alpha() -- items that require the same
# combination of rules within a component share a single alpha draw, just
# as they share a single alpha parameter when the model is actually fit.
# =====================================================================

# Bare column name referenced inside aes() in plot_prior_predictive_check()
# below (NSE; not a missing binding).
utils::globalVariables("value")

# Inverse-CDF sampler for a Normal(mu, sigma) truncated to (0, Inf), i.e.
# the half-Normal family used for alpha in GMLTM()/GMLTM_corr() by default.
.gmltm_rtruncnorm0 <- function(n, mu, sigma) {
  p0 <- stats::pnorm(0, mean = mu, sd = sigma)
  u  <- stats::runif(n, min = p0, max = 1)
  stats::qnorm(u, mean = mu, sd = sigma)
}

#' Prior Predictive Check for the GMLTM-D
#'
#' @description
#' Simulates parameters directly from the priors declared for
#' \code{\link{GMLTM}}/\code{\link{GMLTM_corr}}, simulates response data from
#' those parameters under the GMLTM-D's conjunctive structure, and returns
#' the resulting distribution of simulated proportion-correct -- globally, per
#' item, and (optionally) per component -- across replicates. This lets the
#' analyst check whether a candidate prior implies substantively plausible
#' data \emph{before} spending compute on fitting the model to real data
#' (Gelman et al., 2020, Sect. 2.4).
#'
#' @param Q A \eqn{p \times K} Q-matrix specifying which cognitive rules each
#'   item requires, identical in role to the one passed to
#'   \code{\link{GMLTM}}.
#' @param components A named list grouping rules into components, identical
#'   in role to the one passed to \code{\link{GMLTM}}.
#' @param N Number of simulated examinees per replicate.
#' @param S Number of replicates. Default is \code{1000}, following the
#'   number of prior predictive draws used in Fig. 3 of Gelman et al. (2020);
#'   values from 500 to 1000 are typical.
#' @param priors A named list of prior hyperparameters with the same
#'   structure accepted by \code{\link{GMLTM}}'s \code{priors} argument
#'   (elements \code{theta}, \code{eta}, \code{alpha}, \code{c}), including
#'   \code{priors$alpha$family} (\code{"normal"} or \code{"lognormal"}).
#'   Parameters are simulated under exactly this prior, so the check reflects
#'   what will actually be fit -- not an approximation of it.
#' @param by_component Logical. If \code{TRUE} (default), also returns the
#'   simulated proportion-correct restricted to the items of each component.
#'
#' @details
#' In each of the \code{S} replicates: rule difficulties \eqn{\eta} are drawn
#' as \eqn{\eta_k \sim \text{Normal}(\mu_\eta, \sigma_\eta)} for each rule of
#' each component, and item difficulty is \eqn{\beta = Q \eta}; abilities
#' \eqn{\theta} are drawn i.i.d. \eqn{\text{Normal}(\mu_\theta, \sigma_\theta)}
#' for \code{N} simulated examinees; guessing \eqn{c} is drawn
#' \eqn{\text{Beta}(\text{shape1}, \text{shape2})} per item; and discrimination
#' \eqn{\alpha} is drawn from \code{priors$alpha$family} -- half-Normal
#' (\eqn{\text{Normal}(\mu_\alpha, \sigma_\alpha)} truncated to \eqn{(0,
#' \infty)}) or Log-Normal(\eqn{\mu_\alpha, \sigma_\alpha}) -- with one draw
#' per unique combination of rules that a component's items require (the same
#' grouping \code{\link{GMLTM}} itself uses, so \code{alpha} has as many free
#' values here as it would have as a Stan parameter). Response probabilities
#' follow the GMLTM-D's noncompensatory structure: for item \eqn{i} of
#' examinee \eqn{n}, \eqn{p_{ni} = c_i + (1 - c_i) \prod_{m}
#' \text{plogis}(\alpha_{im}(\theta_{nm} - \beta_{im}))^{C_{im}}}, and
#' \eqn{y_{ni} \sim \text{Bernoulli}(p_{ni})}.
#'
#' Use \code{\link{plot_prior_predictive_check}} to visualize the resulting
#' distribution against a substantively plausible range.
#'
#' @return A list of class \code{"GMLTM_prior_predictive_check"} with elements:
#' \describe{
#'   \item{\code{global}}{Numeric vector of length \code{S}: the overall
#'     simulated proportion correct in each replicate.}
#'   \item{\code{by_item}}{An \code{S x p} matrix of simulated proportion
#'     correct per item.}
#'   \item{\code{by_component}}{An \code{S x M} matrix of simulated proportion
#'     correct restricted to each component's items, or \code{NULL} if
#'     \code{by_component = FALSE}.}
#'   \item{\code{N}, \code{S}}{The arguments used.}
#'   \item{\code{Q}, \code{components}, \code{priors}}{The arguments used
#'     (\code{priors} after applying defaults).}
#' }
#'
#' @references
#' Gelman, A., Vehtari, A., Simpson, D., Margossian, C. C., Carpenter, B.,
#' Yao, Y., Kennedy, L., Gabry, J., Burkner, P.-C., & Modrak, M. (2020).
#' Bayesian workflow. \emph{arXiv}. \doi{10.48550/arXiv.2011.01808}
#'
#' Ramirez, E. S., Jimenez, M., Franco, V. R., & Alvarado, J. M. (2024).
#' Delving into the complexity of analogical reasoning: A detailed exploration
#' with the Generalized Multicomponent Latent Trait Model for Diagnosis.
#' \emph{Journal of Intelligence}, \bold{12}, 67.
#' \doi{10.3390/jintelligence12070067}
#'
#' @examples
#' Q <- structure(
#'   c(0,0,1,0,1,0,1,0,1,1,0,1,1,1,0,1,1,1,0,1,0,1,0,0,1,0,1,
#'     1,0,0,0,0,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,0,
#'     1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,1,0,1,0,1,1,0,1,0,0,0,
#'     0,0,0,0,0,0,1,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,0,
#'     1,0,0,0,1,1,0,1,1,1,1,1,1,1,0,1,1,0,1,1,1,1,0,1),
#'   dim = c(27L, 5L),
#'   dimnames = list(NULL, c("rot_fig","rot_trap","reflection",
#'                           "subt_seg","mov_point")))
#' components <- list(global = c(1, 2, 3), local = c(4, 5))
#'
#' ppc_normal <- prior_predictive_check(Q, components, N = 150, S = 200)
#' summary(ppc_normal$global)
#'
#' ppc_lognormal <- prior_predictive_check(
#'   Q, components, N = 150, S = 200,
#'   priors = list(alpha = list(mu = 0, sigma = 1, family = "lognormal"))
#' )
#' summary(ppc_lognormal$global)
#'
#' @family prior predictive check functions
#' @export
prior_predictive_check <- function(Q, components, N, S = 1000,
                                    priors = list(
                                      theta = list(mu = 0, sigma = 1),
                                      eta   = list(mu = 0, sigma = 1),
                                      alpha = list(mu = 0, sigma = 1, family = "normal"),
                                      c     = list(shape1 = 3, shape2 = 20)
                                    ),
                                    by_component = TRUE) {

  if (!is.numeric(N) || length(N) != 1 || N <= 0 || N != as.integer(N))
    stop("'N' must be a single positive integer.")
  if (!is.numeric(S) || length(S) != 1 || S <= 0 || S != as.integer(S))
    stop("'S' must be a single positive integer.")
  N <- as.integer(N)
  S <- as.integer(S)

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

  p <- nrow(Q)
  K <- ncol(Q)
  M <- length(components)
  comp_names <- names(components)
  if (is.null(comp_names)) comp_names <- paste0("Component", seq_len(M))
  item_names <- if (is.null(rownames(Q))) paste0("item", seq_len(p)) else rownames(Q)

  C <- get_C(Q, components)
  indexes_alpha <- locate_alpha(Q, components)
  n_alpha_unique <- max(indexes_alpha[, 3])

  global_prop <- numeric(S)
  item_prop   <- matrix(NA_real_, nrow = S, ncol = p, dimnames = list(NULL, item_names))
  comp_prop   <- if (by_component) {
    matrix(NA_real_, nrow = S, ncol = M, dimnames = list(NULL, comp_names))
  } else {
    NULL
  }

  for (s in seq_len(S)) {
    eta_matrix <- matrix(0, K, M)
    for (m in seq_len(M)) {
      rules_m <- components[[m]]
      eta_matrix[rules_m, m] <- stats::rnorm(length(rules_m), priors$eta$mu, priors$eta$sigma)
    }
    beta_matrix <- Q %*% eta_matrix

    alpha_vec <- if (priors$alpha$family == "lognormal") {
      stats::rlnorm(n_alpha_unique, meanlog = priors$alpha$mu, sdlog = priors$alpha$sigma)
    } else {
      .gmltm_rtruncnorm0(n_alpha_unique, priors$alpha$mu, priors$alpha$sigma)
    }
    alpha_matrix <- matrix(0, p, M)
    for (r in seq_len(nrow(indexes_alpha))) {
      alpha_matrix[indexes_alpha[r, 1], indexes_alpha[r, 2]] <- alpha_vec[indexes_alpha[r, 3]]
    }

    theta <- matrix(stats::rnorm(N * M, priors$theta$mu, priors$theta$sigma), nrow = N, ncol = M)
    guessing <- stats::rbeta(p, priors$c$shape1, priors$c$shape2)

    prodmu <- matrix(1, N, p)
    for (m in seq_len(M)) {
      active <- which(C[, m] == 1)
      if (length(active) == 0) next
      linear <- outer(theta[, m], beta_matrix[active, m], "-") *
        matrix(alpha_matrix[active, m], nrow = N, ncol = length(active), byrow = TRUE)
      prodmu[, active] <- prodmu[, active] * stats::plogis(linear)
    }

    guess_mat <- matrix(guessing, nrow = N, ncol = p, byrow = TRUE)
    p_mat <- guess_mat + (1 - guess_mat) * prodmu
    y <- matrix(stats::rbinom(N * p, size = 1, prob = as.vector(p_mat)), nrow = N, ncol = p)

    global_prop[s] <- mean(y)
    item_prop[s, ] <- colMeans(y)
    if (by_component) {
      for (m in seq_len(M)) {
        active <- which(C[, m] == 1)
        comp_prop[s, m] <- if (length(active) > 0) mean(y[, active, drop = FALSE]) else NA_real_
      }
    }
  }

  result <- list(
    global = global_prop, by_item = item_prop, by_component = comp_prop,
    N = N, S = S, Q = Q, components = components, priors = priors
  )
  class(result) <- "GMLTM_prior_predictive_check"
  result
}

#' Print method for GMLTM_prior_predictive_check objects
#'
#' @param x An object of class \code{"GMLTM_prior_predictive_check"}, as
#'   returned by \code{\link{prior_predictive_check}}.
#' @param digits Number of decimal digits to display. Default is 3.
#' @param ... Currently ignored.
#'
#' @return Invisibly returns \code{x}.
#' @export
print.GMLTM_prior_predictive_check <- function(x, digits = 3, ...) {
  cat("GMLTM-D prior predictive check\n")
  cat("-------------------------------\n")
  cat(sprintf("Simulations: %d | Subjects: %d | Items: %d\n",
              x$S, x$N, ncol(x$by_item)))
  cat(sprintf(
    "\nOverall proportion correct across simulations: mean = %.*f, [%.*f, %.*f] (2.5%%-97.5%%)\n",
    digits, mean(x$global), digits, stats::quantile(x$global, 0.025),
    digits, stats::quantile(x$global, 0.975)
  ))
  if (!is.null(x$by_component) && !all(is.na(x$by_component))) {
    cat("\nProportion correct by component (mean across simulations):\n")
    print(round(colMeans(x$by_component, na.rm = TRUE), digits))
  }
  cat("\nUse plot_prior_predictive_check(x) to visualize the distribution.\n")
  invisible(x)
}

#' Plot a Prior Predictive Check for the GMLTM-D
#'
#' @description
#' Draws a histogram of the simulated proportion-correct distribution from
#' \code{\link{prior_predictive_check}}, with configurable reference lines
#' marking a substantively plausible range, in the style of Fig. 3 of Gelman
#' et al. (2020): a prior that concentrates simulated data outside that range
#' is a sign the prior is too informative, too diffuse, or otherwise
#' misspecified.
#'
#' @param x An object of class \code{"GMLTM_prior_predictive_check"}, as
#'   returned by \code{\link{prior_predictive_check}}.
#' @param type Which simulated distribution to plot: \code{"global"}
#'   (default, one value per replicate), \code{"by_item"} (item-level values
#'   pooled across items and replicates), or \code{"by_component"}
#'   (component-level values pooled across components and replicates;
#'   requires \code{x} to have been computed with \code{by_component = TRUE}).
#' @param plausible_range Numeric vector of length 2 giving the substantively
#'   plausible range for the proportion correct, drawn as dashed reference
#'   lines. Default is \code{c(0.2, 0.9)}.
#' @param bins Number of histogram bins. Default is \code{30}.
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#' Q <- structure(
#'   c(0,0,1,0,1,0,1,0,1,1,0,1,1,1,0,1,1,1,0,1,0,1,0,0,1,0,1,
#'     1,0,0,0,0,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,0,
#'     1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,1,0,1,0,1,1,0,1,0,0,0,
#'     0,0,0,0,0,0,1,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,0,
#'     1,0,0,0,1,1,0,1,1,1,1,1,1,1,0,1,1,0,1,1,1,1,0,1),
#'   dim = c(27L, 5L),
#'   dimnames = list(NULL, c("rot_fig","rot_trap","reflection",
#'                           "subt_seg","mov_point")))
#' components <- list(global = c(1, 2, 3), local = c(4, 5))
#' ppc <- prior_predictive_check(Q, components, N = 150, S = 200)
#' plot_prior_predictive_check(ppc)
#'
#' @family prior predictive check functions
#' @export
plot_prior_predictive_check <- function(x, type = c("global", "by_item", "by_component"),
                                         plausible_range = c(0.2, 0.9), bins = 30) {
  if (!inherits(x, "GMLTM_prior_predictive_check"))
    stop("'x' must be an object of class 'GMLTM_prior_predictive_check', as returned by prior_predictive_check().")
  type <- match.arg(type)
  if (!is.numeric(plausible_range) || length(plausible_range) != 2)
    stop("'plausible_range' must be a numeric vector of length 2.")

  if (type == "global") {
    values <- x$global
    subtitle_unit <- "replicates"
    plot_title <- "Prior Predictive Check -- Global Proportion Correct"
  } else if (type == "by_item") {
    values <- as.vector(x$by_item)
    subtitle_unit <- "item x replicate values"
    plot_title <- "Prior Predictive Check -- Item Proportion Correct"
  } else {
    if (is.null(x$by_component))
      stop("'x' has no 'by_component' element (it was computed with by_component = FALSE).")
    values <- as.vector(x$by_component)
    subtitle_unit <- "component x replicate values"
    plot_title <- "Prior Predictive Check -- Component Proportion Correct"
  }

  in_range <- mean(values >= plausible_range[1] & values <= plausible_range[2], na.rm = TRUE)
  df <- data.frame(value = values)

  ggplot2::ggplot(df, ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(bins = bins, fill = "#4C72B0", color = "white", alpha = 0.85) +
    ggplot2::geom_vline(xintercept = plausible_range, linetype = "dashed",
                        color = "firebrick", linewidth = 0.8) +
    ggplot2::labs(
      title = plot_title,
      subtitle = sprintf("%.1f%% of %s fall inside the plausible range [%.2f, %.2f]",
                         100 * in_range, subtitle_unit, plausible_range[1], plausible_range[2]),
      x = "Simulated proportion correct", y = "Count"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = 15),
      plot.subtitle = ggplot2::element_text(size = 11, color = "gray30")
    )
}
