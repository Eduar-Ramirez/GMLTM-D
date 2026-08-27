#' @title
#' Posterior Predictive Checks (PPC) for Model Fit Evaluation
#'
#' @description
#'
#' This function generates posterior predictive checks by plotting the observed and simulated
#' total scores distribution, comparing empirical data against model predictions. It also
#' computes fitted values and prediction intervals.
#'
#' @usage
#' ppchecks(fit, nsim = 100, interval = 0.95, ...)
#'
#' @param fit A fitted MLTM object containing model results.
#' @param nsim Number of simulated posterior samples (default = 100).
#' @param interval Probability associated with the credibility intervals (default = 0.95).
#' @param ... Additional graphical parameters to customize the plot.
#'
#' @details
#' The function simulates multiple datasets from the posterior distribution and compares
#' the empirical distribution of total scores with the predicted distribution.
#' It overlays the observed and predicted distributions using a histogram with transparency.
#'
#' The fitted values, along with their credibility intervals, are computed and returned.
#'
#' @return
#' A list containing:
#' \itemize{
#'   \item \code{items}: A table of fitted values and prediction intervals for each item.
#'   \item \code{subjects}: Fitted and observed mean scores per subject.
#'   \item \code{ysim}: Simulated response matrices.
#' }
#'
#' The function also generates:
#' \itemize{
#'   \item A histogram comparing observed vs. predicted total scores.
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
#'   fit <- LLTM(analogy, Q, iters = 200, iter_warmup = 100, chains = 1, cores = 1)
#'   ppchecks(fit)
#' }
#'
#' @family posterior predictive check functions
#' @export
ppchecks <- function(fit, nsim = 100, interval = 0.95, ...) {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)

  # Extract data dimensions
  n <- nrow(fit$data)
  p <- ncol(fit$data)
  scores <- rowSums(fit$data)

  # Simulate response matrices from posterior predictive distribution
  ysim <- array(NA, dim = c(nsim, n, p))
  probs <- fit$posterior$probabilities

  for(i in 1:nsim) {
    ysim[i, , ] <- rbinom(n * p, size = 1, prob = c(probs[i, , ]))
  }

  scores_ysim <- apply(ysim, MARGIN = 2, rowSums)

  # Create histograms without plotting
  p1 <- hist(scores, breaks = 0:p, plot = FALSE)
  p2 <- hist(c(scores_ysim), breaks = 0:p, plot = FALSE)

  # Determine y-axis limits dynamically
  ymax <- max(c(p1$density, p2$density)) * 1.2

  # Compute median of the observed scores
  median_score <- median(scores)

  # Decide the legend position based on the median score
  legend_position <- ifelse(median_score < (p / 2), "topright", "topleft")

  # Configure the plot margins with extra space on the right
  par(mar = c(5, 5, 2, 4))

  # Plot the histogram of observed scores
  plot(p1, col = rgb(0, 0, 1, 1/4), freq = FALSE, ylim = c(0, ymax),
       xlab = "Total Scores", main = "", xlim = c(0, p),
       cex.axis = 1.2, cex.lab = 1.4, lwd = 2, ...)

  # Overlay the histogram of predicted scores
  plot(p2, col = rgb(1, 0, 0, 1/4), freq = FALSE, add = TRUE)

  # Add grid lines for better visualization
  abline(h = seq(0, ymax, length.out = 6), col = "gray90", lty = "dotted")

  # Define colors for the legend
  cols2 <- c(rgb(1, 0, 0, 1/4), rgb(0, 0, 1, 1/4), rgb(0.36, 0, 0.64, 0.5))

  # Place the legend dynamically based on the data distribution
  legend(legend_position, legend = c("Predicted", "Observed", "Overlap"),
         fill = cols2, cex = 1.2, bty = "n",
         title = "Posterior Predictive Check", title.font = 2)

  return(invisible(NULL))
}


utils::globalVariables(c("ID", "lower", "upper", "predicted",
                          "observed", "obs_color"))

#' @title
#' Marginal Proportions Predictive Checks
#'
#' @description
#'
#' Computes and visualizes marginal success proportions, including predicted values,
#' confidence intervals, RMSR, SRMR, and bias estimation.
#'
#' @usage
#'
#' marginal_Pchecks(fit, interval = 0.95)
#'
#' @param fit MLTM object containing model results.
#' @param interval Probability associated with the credibility intervals (default = 0.95).
#'
#' @details \code{marginal_Pchecks} calculates marginal prediction intervals and observed success proportions.
#' It prints a table with observed vs. predicted values, generates a forest plot for visualization, and
#' computes key fit indices: RMSR, SRMR, and bias.
#'
#' @return A list containing:
#' \itemize{
#'  \item \code{items}: A table of fitted values and prediction intervals for each item.
#'  \item \code{rmsr}: The Root Mean Square Residual (RMSR).
#'  \item \code{srmr}: The Standardized Root Mean Square Residual (SRMR).
#'  \item \code{bias}: The difference between the total observed and predicted proportions.
#' }
#'
#' The function also generates:
#' \itemize{
#'  \item A forest plot visualizing prediction intervals and observed success probabilities.
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
#'   fit <- LLTM(analogy, Q, iters = 200, iter_warmup = 100, chains = 1, cores = 1)
#'   marginal_Pchecks(fit)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_segment geom_point
#'   scale_color_identity labs scale_y_continuous coord_cartesian
#'   theme_minimal theme element_text unit element_line element_rect margin
#' @family posterior predictive check functions
#' @export
marginal_Pchecks <- function(fit, interval = 0.95) {

  if (!"posterior" %in% names(fit) || !"probabilities" %in% names(fit$posterior)) {
    stop("Error: The provided model object does not contain posterior probabilities.")
  }

  # Compute fitted and observed success probabilities
  p <- fit$posterior$probabilities
  EAP_p <- apply(p, MARGIN = c(2, 3), mean)  # Expected a posteriori (EAP) estimates
  p_obs <- colMeans(fit$data)  # Observed success proportions
  p_pred <- colMeans(EAP_p)  # Predicted success probabilities
  marginals <- apply(p, MARGIN = c(1, 3), mean)  # Item-level marginals

  # Compute confidence intervals
  CI_p <- t(apply(marginals, MARGIN = 2, quantile, probs = c((1 - interval) / 2, 1 - (1 - interval) / 2)))

  # Create result table
  fitted_items <- cbind(CI_p[, 1], p_pred, CI_p[, 2], p_obs)
  colnames(fitted_items) <- c("lower", "predicted", "upper", "observed")

  # Compute RMSR (Root Mean Square Residual)
  residuals <- fitted_items[, "observed"] - fitted_items[, "predicted"]
  rmsr <- sqrt(mean(residuals^2))

  # Compute SRMR (Standardized Root Mean Square Residual)
  srmr <- sqrt(mean((residuals / sd(residuals))^2))

  # Compute Bias (Sum of observed - sum of predicted)
  bias <- sum(fitted_items[, "observed"]) - sum(fitted_items[, "predicted"])

  # Print results
  print(fitted_items)
  cat("\nPosterior Predictive Discrepancy Measures for Marginal Proportions:\n")
  cat("RMSR:", round(rmsr, 6), "\n")
  cat("SRMR:", round(srmr, 6), "\n")
  cat("Bias:", round(bias, 6), "\n")

  # Create a data frame for plotting
  df <- as.data.frame(fitted_items)
  df$ID <- 1:nrow(df)
  ybreaks <- seq(1, nrow(df), 2)
  ylabels <- rownames(df)[c(TRUE, FALSE)]
  df$obs_color <- ifelse(df$observed >= df$lower & df$observed <= df$upper, "lightgreen", "salmon")

  # Create forest plot
  p_forest <- ggplot(df, aes(y = ID)) +
    geom_segment(aes(x = lower, xend = upper, y = ID, yend = ID), color = "black", linewidth = 1) +
    geom_point(aes(x = predicted), color = "black", size = 3) +
    scale_color_identity() +
    geom_point(aes(x = observed, color = obs_color), size = 3) +
    labs(x = "Success probability", y = "Item", title = "Marginal PPchecks (Items)") +
    scale_y_continuous(breaks = ybreaks, labels = ylabels) +
    coord_cartesian(xlim = c(min(df[, 1:4]), max(df[, 1:4]))) +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 12, margin = unit(c(3, 0, 0, 0), "mm")),
          axis.title.x = element_text(size = 14, margin = unit(c(5, 0, 0, 0), "mm")),
          axis.text.y = element_text(size = 12, margin = unit(c(0, 3, 0, 0), "mm")),
          axis.title.y = element_text(size = 14, margin = unit(c(0, 5, 0, 0), "mm")),
          axis.ticks.x = element_line(size = 0.6),
          axis.ticks.y = element_line(size = 0.6),
          panel.border = element_rect(fill = NA),
          axis.ticks.length = unit(0.2, "cm"),
          plot.margin = margin(0.5, 1, 0.2, 0.5, 'cm'),
          plot.title = element_text(hjust = 0.5, margin = unit(c(0, 0, 3, 0), "mm")),
          axis.line = element_line(size = 0.5, colour = "black"))

  print(p_forest)

  return(invisible(list(items = fitted_items, rmsr = rmsr, srmr = srmr, bias = bias)))
}


