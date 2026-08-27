#' @title Grouped Item Characteristic Curves (ICC) Plot
#'
#' @description
#' Generates Item Characteristic Curves (ICCs) for a group of items displayed in a
#' grid layout with a shared legend.
#'
#' @param fit A fitted GMLTM model object containing EAP parameter estimates.
#' @param Q The Q-matrix indicating the association between items and rules.
#' @param components A list where each element is a vector of rule indices per component.
#' @param page Integer specifying which page of items to display.
#' @param n_items_per_page Number of items to include per page. Default is 9.
#' @param ncol Number of columns in the layout grid. Default is 3.
#' @param nrow Number of rows in the layout grid. Default is 3.
#'
#' @details
#' Displays one legend shared across all plots and ensures consistency across theta
#' and probability axes. Ideal for publications or appendices.
#'
#' @return A composed ICC grid with one shared legend, plotted to the active device.
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
#'   plot_ICC_grouped(fit, Q, components, page = 1)
#' }
#'
#' @family item characteristic curve plots
#' @export
plot_ICC_grouped <- function(fit, Q, components, page = 1, n_items_per_page = 9,
                             ncol = 3, nrow = 3) {
  theta    <- seq(-4, 4, length.out = 100)
  num_items <- nrow(fit$EAP$alpha)
  cols2    <- c("#1f77b4", "#d62728")

  K <- ncol(Q)
  M <- length(components)
  eta_mat <- matrix(0, K, M)
  for (i in seq_along(components)) eta_mat[components[[i]], i] <- 1
  C <- (Q %*% eta_mat > 0) * 1

  start_item <- (page - 1) * n_items_per_page + 1
  end_item   <- min(page * n_items_per_page, num_items)

  plots <- list()
  p1 <- p2 <- NULL

  for (j in start_item:end_item) {
    mu1 <- plogis(fit$EAP$alpha[j, 1] * (theta - fit$EAP$beta[j, 1]))
    if (C[j, 1] == 0) mu1 <- rep(0, length(theta))
    p1 <- fit$EAP$guessing[j] + (1 - fit$EAP$guessing[j]) * mu1

    mu2 <- plogis(fit$EAP$alpha[j, 2] * (theta - fit$EAP$beta[j, 2]))
    if (C[j, 2] == 0) mu2 <- rep(0, length(theta))
    p2 <- fit$EAP$guessing[j] + (1 - fit$EAP$guessing[j]) * mu2

    df <- data.frame(theta = theta, p1 = p1, p2 = p2)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = theta)) +
      ggplot2::geom_line(ggplot2::aes(y = p1, color = "Component 1"), linewidth = 1.3) +
      ggplot2::geom_line(ggplot2::aes(y = p2, color = "Component 2"), linewidth = 1.3) +
      ggplot2::scale_color_manual(values = cols2, name = NULL) +
      ggplot2::scale_x_continuous(limits = c(-4, 4)) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(title = paste("Item", j),
                    x = expression(theta), y = expression(P(theta))) +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::theme(
        plot.title  = ggplot2::element_text(hjust = 0.5, face = "bold", size = 15),
        axis.title  = ggplot2::element_text(size = 13),
        axis.text   = ggplot2::element_text(size = 11),
        legend.position = "none"
      )

    plots[[j - start_item + 1]] <- p
  }

  df_legend <- data.frame(theta = theta, p1 = p1, p2 = p2)
  legend_plot <- ggplot2::ggplot(df_legend, ggplot2::aes(x = theta)) +
    ggplot2::geom_line(ggplot2::aes(y = p1, color = "Component 1"), linewidth = 1.3) +
    ggplot2::geom_line(ggplot2::aes(y = p2, color = "Component 2"), linewidth = 1.3) +
    ggplot2::scale_color_manual(values = cols2, name = NULL) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      legend.position  = "bottom",
      legend.key.size  = grid::unit(1.3, "lines"),
      legend.text      = ggplot2::element_text(size = 13),
      legend.background = ggplot2::element_rect(fill = "gray95", color = NA)
    )

  get_legend <- function(plot) {
    g <- ggplot2::ggplotGrob(plot)
    legend_index <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
    g$grobs[[legend_index]]
  }
  legend <- get_legend(legend_plot)

  gridExtra::grid.arrange(
    do.call(gridExtra::arrangeGrob, c(plots, ncol = ncol, nrow = nrow)),
    legend,
    nrow    = 2,
    heights = c(10, 1.3)
  )
}

#' @title Individual Item Characteristic Curves (ICC)
#'
#' @description Returns a list of individual ICC \code{ggplot2} plots (one per item).
#'
#' @param fit A fitted GMLTM model object.
#' @param Q The Q-matrix for rule-item associations.
#' @param components A list indicating rule groupings per component.
#'
#' @return A list of individual \code{ggplot} objects, one per item.
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
#'   plots <- plot_ICC_individual(fit, Q, components)
#'   print(plots[[1]])
#' }
#'
#' @family item characteristic curve plots
#' @export
plot_ICC_individual <- function(fit, Q, components) {
  theta     <- seq(-4, 4, length.out = 100)
  num_items <- nrow(fit$EAP$alpha)
  cols2     <- c("#1f77b4", "#d62728")

  K <- ncol(Q)
  M <- length(components)
  eta_mat <- matrix(0, K, M)
  for (i in seq_along(components)) eta_mat[components[[i]], i] <- 1
  C <- (Q %*% eta_mat > 0) * 1

  plots <- list()

  for (j in seq_len(num_items)) {
    mu1 <- plogis(fit$EAP$alpha[j, 1] * (theta - fit$EAP$beta[j, 1]))
    if (C[j, 1] == 0) mu1 <- rep(0, length(theta))
    p1 <- fit$EAP$guessing[j] + (1 - fit$EAP$guessing[j]) * mu1

    mu2 <- plogis(fit$EAP$alpha[j, 2] * (theta - fit$EAP$beta[j, 2]))
    if (C[j, 2] == 0) mu2 <- rep(0, length(theta))
    p2 <- fit$EAP$guessing[j] + (1 - fit$EAP$guessing[j]) * mu2

    df <- data.frame(theta = theta, p1 = p1, p2 = p2)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = theta)) +
      ggplot2::geom_line(ggplot2::aes(y = p1, color = "Component 1"), linewidth = 1) +
      ggplot2::geom_line(ggplot2::aes(y = p2, color = "Component 2"), linewidth = 1) +
      ggplot2::scale_color_manual(values = cols2, name = NULL) +
      ggplot2::labs(title = paste("Item", j),
                    x = expression(theta), y = expression(P(theta))) +
      ggplot2::scale_x_continuous(limits = c(-4, 4)) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::theme_minimal(base_size = 16) +
      ggplot2::theme(
        plot.title      = ggplot2::element_text(hjust = 0.5, face = "bold", size = 18),
        legend.position = "bottom"
      )

    plots[[j]] <- p
  }

  return(plots)
}
