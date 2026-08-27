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
#' @param fit A fitted model object of class \code{"GMLTM"}, \code{"MLTM"}, or
#'   \code{"LLTM"}, as returned by \code{\link{GMLTM}}, \code{\link{MLTM}}, or
#'   \code{\link{LLTM}}.
#'
#' @details \code{reliability} computes the classical marginal reliability
#'   coefficient for each cognitive component's ability estimates, from the
#'   posterior mean and posterior variance of theta across subjects:
#'   \eqn{r_{xx} = s^2 / (s^2 + e)}, where \eqn{s^2} is the variance of the
#'   EAP theta estimates across subjects (true-score variance) and \eqn{e} is
#'   the average posterior variance of theta across subjects (error variance).
#'
#' @return A named numeric vector with one marginal reliability coefficient
#'   per cognitive component (length 1 for \code{LLTM} fits, which have a
#'   single ability dimension).
#'
#' @references
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
#'   fit <- GMLTM(analogy, Q, components,
#'                iters = 200, iter_warmup = 100, chains = 1, cores = 1)
#'   reliability(fit)
#' }
#'
#' @family reliability functions
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
