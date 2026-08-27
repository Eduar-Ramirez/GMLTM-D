#' Extract and Summarize the Latent Correlation Matrix from a GMLTM_corr Fit
#'
#' @description
#' Extracts the posterior correlation matrix \eqn{\Sigma} among cognitive
#' components from a model fitted with \code{\link{GMLTM_corr}}, together with
#' pairwise credible intervals and a heatmap for visual inspection.
#'
#' @param fit An object of class \code{"GMLTM_corr"}, as returned by
#'   \code{\link{GMLTM_corr}}.
#' @param credible Width of the credible interval for each pairwise
#'   correlation, e.g. \code{0.95} for a 95\% interval. Default is \code{0.95}.
#' @param plot Logical. If \code{TRUE} (default), draws a correlation heatmap
#'   of the posterior mean \eqn{\Sigma} on the active graphics device, using
#'   the same color scheme as the correlation diagnostics in
#'   \code{\link{generate_Q_with_interactions}}.
#'
#' @details
#' The pairwise credible intervals are computed from the posterior draws of
#' \code{Sigma} in \code{fit$fit}, extracted with \code{rstan::extract()} (not
#' from \code{fit$EAP$Sigma}, which only holds the posterior mean). A pair of
#' components is flagged as showing evidence of a real (non-zero) correlation
#' when its credible interval excludes 0.
#'
#' @return A list of class \code{"GMLTM_correlation"} with elements:
#' \describe{
#'   \item{\code{Sigma}}{The \code{M x M} posterior mean correlation matrix
#'     (\code{fit$EAP$Sigma}).}
#'   \item{\code{pairs}}{A data.frame with one row per pair of components:
#'     \code{component_1}, \code{component_2}, \code{estimate}, \code{lower},
#'     \code{upper}, and \code{excludes_zero} (logical).}
#'   \item{\code{credible}}{The credible-interval width used.}
#'   \item{\code{plot}}{A function that redraws the heatmap when called with
#'     no arguments, or \code{NULL} if \code{plot = FALSE}.}
#' }
#'
#' @references
#' Embretson, S. E., & Yang, X. (2013). A multicomponent latent trait model
#' for diagnosis. \emph{Psychometrika}, \bold{78}, 14--36.
#' \doi{10.1007/s11336-012-9296-y}
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
#'   corr <- extract_correlation(fit)
#'   corr$pairs
#' }
#'
#' @importFrom stats quantile
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics image contour axis legend par
#' @family diagnostic reporting functions
#' @export
extract_correlation <- function(fit, credible = 0.95, plot = TRUE) {
  if (!inherits(fit, "GMLTM_corr"))
    stop("'fit' must be an object of class 'GMLTM_corr', as returned by GMLTM_corr().")
  if (!is.numeric(credible) || length(credible) != 1 || credible <= 0 || credible >= 1)
    stop("'credible' must be a single number strictly between 0 and 1.")

  Sigma_EAP  <- fit$EAP$Sigma
  comp_names <- colnames(Sigma_EAP)
  M <- ncol(Sigma_EAP)

  Sigma_draws <- rstan::extract(fit$fit, pars = "Sigma")$Sigma
  alpha_level <- 1 - credible
  probs <- c(alpha_level / 2, 1 - alpha_level / 2)

  pairs <- data.frame(
    component_1 = character(0), component_2 = character(0),
    estimate = numeric(0), lower = numeric(0), upper = numeric(0)
  )
  for (i in seq_len(M - 1)) {
    for (j in seq(i + 1, M)) {
      draws_ij <- Sigma_draws[, i, j]
      q <- stats::quantile(draws_ij, probs = probs)
      pairs <- rbind(pairs, data.frame(
        component_1 = comp_names[i], component_2 = comp_names[j],
        estimate = mean(draws_ij), lower = q[1], upper = q[2]
      ))
    }
  }
  rownames(pairs) <- NULL
  pairs$excludes_zero <- pairs$lower > 0 | pairs$upper < 0

  draw_heatmap <- function() {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)
    par(mar = c(5, 5, 4, 9))

    colors <- colorRampPalette(c("blue", "white", "red"))(100)
    image(1:ncol(Sigma_EAP), 1:nrow(Sigma_EAP), as.matrix(Sigma_EAP),
          col = colors, zlim = c(-1, 1),
          main = "GMLTM-D: Latent Component Correlation (Sigma)",
          xlab = "", ylab = "", axes = FALSE)
    contour(1:ncol(Sigma_EAP), 1:nrow(Sigma_EAP), as.matrix(Sigma_EAP),
            levels = c(-0.8, -0.5, 0.5, 0.8), add = TRUE, col = "black", lwd = 1)
    axis(1, at = 1:ncol(Sigma_EAP), labels = comp_names, cex.axis = 0.9, las = 2)
    axis(2, at = 1:nrow(Sigma_EAP), labels = comp_names, cex.axis = 0.9, las = 2)
    legend(x = ncol(Sigma_EAP) + 0.7, y = nrow(Sigma_EAP),
           legend = c("1.0", "0.5", "0.0", "-0.5", "-1.0"),
           fill = colors[c(100, 75, 50, 25, 1)],
           title = "Correlation", xpd = TRUE, bty = "n", cex = 0.85)
  }

  plot_fn <- NULL
  if (isTRUE(plot)) {
    draw_heatmap()
    plot_fn <- draw_heatmap
  }

  result <- list(Sigma = Sigma_EAP, pairs = pairs, credible = credible, plot = plot_fn)
  class(result) <- "GMLTM_correlation"
  return(result)
}

#' Print method for GMLTM_correlation objects
#'
#' @param x An object of class \code{"GMLTM_correlation"}, as returned by
#'   \code{\link{extract_correlation}}.
#' @param digits Number of decimal digits to display. Default is 3.
#' @param ... Currently ignored.
#'
#' @return Invisibly returns \code{x}.
#' @export
print.GMLTM_correlation <- function(x, digits = 3, ...) {
  cat("GMLTM-D latent component correlation (posterior mean Sigma)\n")
  cat("--------------------------------------------------------\n")
  print(round(x$Sigma, digits))
  cat("\n")
  pct <- round(100 * x$credible)
  cat(sprintf("Pairwise correlations (%d%% credible interval):\n", pct))
  for (i in seq_len(nrow(x$pairs))) {
    row <- x$pairs[i, ]
    verdict <- if (isTRUE(row$excludes_zero)) {
      "interval excludes 0 -> evidence of a real correlation"
    } else {
      "interval includes 0 -> no clear evidence of correlation"
    }
    cat(sprintf(
      "  %s - %s: r = %.*f, [%.*f, %.*f] -- %s\n",
      row$component_1, row$component_2,
      digits, row$estimate, digits, row$lower, digits, row$upper, verdict
    ))
  }
  invisible(x)
}
