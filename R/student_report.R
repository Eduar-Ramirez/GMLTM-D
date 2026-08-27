# =====================================================================
# student_report(): individual diagnostic report for a student from
# a GMLTM/GMLTM_corr object, following the mastery cutline and decision
# confidence procedure described in Embretson (2019, Chap. 9, Sects.
# 9.2.2 and 9.4.2) for the MLTM-D, and inspired by Figs. 4 and 5 of
# Embretson & Yang (2013) for the graphical representation of the rule
# difficulty continuum.
#
# In addition to the component-level cutline gamma_m, a rule-level
# cutline tau_km is computed for each rule k of component m, restricted
# to the items that require that specific rule within the component
# (Embretson, 2019, Sect. 9.2.2).
# =====================================================================

# The mean-probability-equals-y root-finding kernel used for both the
# component-level cutline gamma_m and the rule-level cutline tau_km.
# alpha_m/beta_m/guess_m are already the item subset the caller wants the
# cutline computed over (all of a component's items for gamma_m; only the
# items requiring one specific rule within the component for tau_km) --
# the kernel itself is agnostic to which grouping that subset represents.
# unit_label only changes the wording of the warnings below, so a
# rule-level call reads "Rule '...'" instead of "Component '...'".

# Bare column names referenced inside aes() in .gmltm_plot_component() and
# .gmltm_plot_badges() below (NSE; not missing bindings).
utils::globalVariables(c("difficulty", "rule", "theta", "label_y", "x", "y",
                          "label", "badge_label", "badge_color"))

.gmltm_component_gamma <- function(y, alpha_m, beta_m, guess_m,
                                    comp_name, lower = -10, upper = 10,
                                    unit_label = "Component") {
  f <- function(theta) {
    mean(guess_m + (1 - guess_m) * plogis(alpha_m * (theta - beta_m))) - y
  }
  f_lo <- f(lower)
  f_hi <- f(upper)
  if (f_lo >= 0) {
    warning(sprintf(
      "%s '%s': the mean probability is already >= y at theta = %d; gamma set to -Inf (mastery guaranteed).",
      unit_label, comp_name, lower))
    return(-Inf)
  }
  if (f_hi <= 0) {
    warning(sprintf(
      "%s '%s': the mean probability never reaches y in [%d, %d]; gamma set to Inf (mastery unreachable).",
      unit_label, comp_name, lower, upper))
    return(Inf)
  }
  stats::uniroot(f, lower = lower, upper = upper)$root
}

# Rule-level cutline tau_km: same root-finding kernel as
# .gmltm_component_gamma(), restricted to the items that both load on
# component 'comp_idx' and require rule 'rule_idx' specifically
# (which(Q[, rule_idx] == 1 & C[, comp_idx] == 1)). If no such item
# exists -- the rule never appears on its own within the component's
# items -- tau_km is undefined: a warning is issued and NA is returned,
# so the caller can omit that rule instead of failing.
.gmltm_rule_tau <- function(y, fit, Q, C, comp_idx, rule_idx, rule_name,
                             lower = -10, upper = 10) {
  items_km <- which(Q[, rule_idx] == 1 & C[, comp_idx] == 1)
  if (length(items_km) == 0) {
    warning(sprintf(
      "Rule '%s': no items in this component require this rule specifically; tau is undefined and will be omitted from the report.",
      rule_name))
    return(NA_real_)
  }
  .gmltm_component_gamma(
    y = y,
    alpha_m = fit$EAP$alpha[items_km, comp_idx],
    beta_m  = fit$EAP$beta[items_km, comp_idx],
    guess_m = fit$EAP$guessing[items_km],
    comp_name = rule_name,
    lower = lower, upper = upper,
    unit_label = "Rule"
  )
}

.gmltm_declutter <- function(vals, min_gap) {
  ord <- order(vals)
  v <- vals[ord]
  for (i in seq_along(v)[-1]) {
    if (v[i] - v[i - 1] < min_gap) v[i] <- v[i - 1] + min_gap
  }
  out <- numeric(length(vals))
  out[ord] <- v
  out
}

#' @keywords internal
.gmltm_plot_component <- function(comp_name, rule_diff, theta_eap, theta_lo, theta_hi,
                                   gamma, mastery, y, digits) {
  rules_df <- data.frame(rule = names(rule_diff), difficulty = as.numeric(rule_diff))

  y_vals <- c(rules_df$difficulty, theta_eap, theta_lo, theta_hi)
  if (is.finite(gamma)) y_vals <- c(y_vals, gamma)
  y_range <- range(y_vals, na.rm = TRUE) + c(-0.5, 0.5)

  min_gap <- diff(y_range) * 0.07
  rules_df$label_y <- .gmltm_declutter(rules_df$difficulty, min_gap)

  color     <- if (mastery == 1) "#2ca02c" else "#d62728"
  state_lab <- if (mastery == 1) "Mastery" else "Non-mastery"
  gamma_lab <- if (is.finite(gamma)) sprintf("%.*f", digits, gamma) else as.character(gamma)

  student_df <- data.frame(x = 1.5, theta = theta_eap, lower = theta_lo, upper = theta_hi)

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(ggplot2::aes(x = 1, xend = 1, y = y_range[1], yend = y_range[2]),
                          color = "gray70", linewidth = 0.6)

  if (is.finite(gamma)) {
    p <- p + ggplot2::geom_hline(yintercept = gamma, linetype = "dashed", color = "gray40")
  }

  p <- p +
    ggplot2::geom_segment(data = rules_df,
                          ggplot2::aes(x = 1, xend = 0.82, y = difficulty, yend = label_y),
                          color = "gray70", linewidth = 0.3) +
    ggplot2::geom_point(data = rules_df,
                        ggplot2::aes(x = 1, y = difficulty), shape = 45, size = 6) +
    ggplot2::geom_text(data = rules_df,
                       ggplot2::aes(x = 0.78, y = label_y, label = rule),
                       hjust = 1, size = 3.4) +
    ggplot2::geom_pointrange(data = student_df,
                             ggplot2::aes(x = x, y = theta, ymin = lower, ymax = upper),
                             color = color, linewidth = 0.9, size = 0.8)

  if (is.finite(gamma)) {
    p <- p + ggplot2::annotate("label", x = 1.5, y = gamma,
                               label = sprintf("cutline (y = %.2f)", y),
                               vjust = -0.4, size = 3.2, color = "gray30",
                               fill = "white", alpha = 0.85)
  }

  p +
    ggplot2::scale_x_continuous(limits = c(0.1, 1.9), breaks = NULL) +
    ggplot2::scale_y_continuous(limits = y_range) +
    ggplot2::labs(
      title = sprintf("Component: %s (%s)", comp_name, state_lab),
      subtitle = sprintf("theta EAP = %.*f  [%.*f, %.*f]   gamma = %s",
                         digits, theta_eap, digits, theta_lo, digits, theta_hi, gamma_lab),
      x = NULL, y = expression(theta)
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = 15),
      plot.subtitle = ggplot2::element_text(size = 11, color = "gray30")
    )
}

# A simpler, "at a glance" complement to .gmltm_plot_component(): one row
# per component and, indented below it, one row per rule with a valid
# tau_km, each with its theta/tau value and a green ("1"/mastery) or red
# ("0"/non-mastery) badge. Rules omitted from tau_list (no items of their
# own within the component) simply do not get a row.
#' @keywords internal
.gmltm_plot_badges <- function(comp_names, theta_eap, domain, tau_list,
                                domain_rules_list, digits) {
  rows <- list()
  y_pos <- 0
  for (m in seq_along(comp_names)) {
    y_pos <- y_pos + 1
    rows[[length(rows) + 1]] <- data.frame(
      y = y_pos, level = "component",
      label = comp_names[m],
      value = sprintf("theta = %.*f", digits, theta_eap[m]),
      badge = domain[m],
      stringsAsFactors = FALSE
    )

    rule_names_m <- names(tau_list[[comp_names[m]]])
    for (k in seq_along(rule_names_m)) {
      y_pos <- y_pos + 1
      rows[[length(rows) + 1]] <- data.frame(
        y = y_pos, level = "rule",
        label = paste0("    ", rule_names_m[k]),
        value = sprintf("tau = %.*f", digits, tau_list[[comp_names[m]]][k]),
        badge = domain_rules_list[[comp_names[m]]][k],
        stringsAsFactors = FALSE
      )
    }
  }

  df <- do.call(rbind, rows)
  df$y <- max(df$y) - df$y + 1  # first component/rule at the top
  df$badge_color <- ifelse(df$badge == 1, "#2ca02c", "#d62728")
  df$badge_label <- ifelse(df$badge == 1, "1", "0")

  comp_rows <- df[df$level == "component", ]
  rule_rows <- df[df$level == "rule", ]

  ggplot2::ggplot(df, ggplot2::aes(y = y)) +
    ggplot2::geom_point(ggplot2::aes(x = 0, fill = I(badge_color)),
                        size = 9, shape = 21, color = "white", stroke = 0.8) +
    ggplot2::geom_text(ggplot2::aes(x = 0, label = badge_label),
                       color = "white", size = 4, fontface = "bold") +
    ggplot2::geom_text(data = comp_rows, ggplot2::aes(x = 0.3, label = label),
                       hjust = 0, size = 4.2, fontface = "bold") +
    ggplot2::geom_text(data = rule_rows, ggplot2::aes(x = 0.3, label = label),
                       hjust = 0, size = 3.6, color = "gray20") +
    ggplot2::geom_text(ggplot2::aes(x = 3.6, label = value),
                       hjust = 1, size = 3.3, color = "gray40") +
    ggplot2::scale_x_continuous(limits = c(-0.3, 3.8)) +
    ggplot2::scale_y_continuous(limits = c(0.3, max(df$y) + 0.7)) +
    ggplot2::labs(title = "Mastery summary", x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      axis.text  = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 15)
    )
}

#' Individual Student Diagnostic Report for GMLTM-D Models
#'
#' @description
#' Builds a per-component and per-rule diagnostic report for a single
#' examinee from a fitted \code{\link{GMLTM}}/\code{\link{GMLTM_corr}}
#' model: posterior EAP and credible interval of \eqn{\theta}, a mastery
#' cutline \eqn{\gamma_m} per component and \eqn{\tau_{km}} per rule,
#' derived from a mastery probability \code{y}, binary mastery vectors, a
#' decision-confidence index at both levels, and the difficulty of the
#' rules/items involved in each component. Follows the mastery-diagnosis
#' procedure for the MLTM-D described in Embretson (2019, Chap. 9, Sects.
#' 9.2.2 and 9.4.2).
#'
#' @param fit An object of class \code{"GMLTM"} (as returned by
#'   \code{\link{GMLTM}} or \code{\link{GMLTM_corr}}).
#' @param Q A \eqn{p \times K} Q-matrix, identical to the one used to fit
#'   \code{fit}, specifying which rules each item requires.
#' @param components A named list grouping rules into components, identical
#'   to the one used to fit \code{fit}.
#' @param student_id Identifier of the examinee: either a row name of
#'   \code{fit$data} (character/numeric matching \code{rownames(fit$data)})
#'   or an integer row index in \code{1:nrow(fit$data)}.
#' @param y Mastery probability used to derive the cutlines \eqn{\gamma_m}
#'   and \eqn{\tau_{km}} (Embretson, 2019, Sect. 9.2.2). Default is
#'   \code{0.50}. Ignored for components overridden via \code{theta_cut}
#'   (\eqn{\tau_{km}} is always derived from \code{y}; see Details).
#' @param theta_cut Optional fixed cutline(s) on the \eqn{\theta} scale,
#'   used instead of deriving \eqn{\gamma_m} from \code{y}. Either a single
#'   number (applied to every component) or a numeric vector of length
#'   \code{length(components)} (optionally named with the component names).
#'   Default is \code{NULL} (derive \eqn{\gamma_m} from \code{y}). Does not
#'   affect the rule-level cutlines \eqn{\tau_{km}}.
#' @param credible Width of the credible interval for \eqn{\theta}, e.g.
#'   \code{0.95} for a 95\% interval. Default is \code{0.95}.
#' @param digits Number of decimal digits used in the printed summary and
#'   in the plot subtitles. Default is \code{3}.
#' @param include_plots Logical. If \code{TRUE} (default), builds the
#'   \code{ggplot2} objects described under \code{plots} below. If
#'   \code{FALSE}, \code{plots} is \code{NULL} and that (comparatively
#'   expensive) step is skipped, leaving only the numeric summary; this is
#'   used by \code{\link{student_report_batch}} to compute a lightweight
#'   summary for every requested student while only rendering plots for the
#'   current page.
#'
#' @details
#' For each component \eqn{m}, the mastery cutline \eqn{\gamma_m} is the
#' value of \eqn{\theta_m} at which the mean predicted probability of
#' solving the items that load on component \eqn{m},
#' \eqn{\bar{P}_m(\theta_m) = \text{mean}_{k}\left[c_k + (1-c_k)\,
#' \text{plogis}(\alpha_{km}(\theta_m - \beta_{km}))\right]}, equals
#' \code{y} (Embretson, 2019, Eq. in Sect. 9.2.2). It is found numerically
#' with \code{\link[stats]{uniroot}} using the item EAP estimates
#' \code{fit$EAP$alpha}, \code{fit$EAP$beta}, and \code{fit$EAP$guessing}.
#' If \code{theta_cut} is supplied, it is used directly as \eqn{\gamma_m}
#' and \code{y} is not used to derive it.
#'
#' In addition, a rule-level cutline \eqn{\tau_{km}} is computed for each
#' rule \eqn{k} belonging to component \eqn{m}, using the same root-finding
#' procedure as \eqn{\gamma_m} but restricted to the items that both load on
#' component \eqn{m} and require rule \eqn{k} specifically (Embretson, 2019,
#' Sect. 9.2.2). If no such item exists -- the rule never appears on its own
#' among the component's items -- \eqn{\tau_{km}} is undefined: a warning is
#' issued and that rule is omitted from \code{domain_rules} and
#' \code{confidence_rules} for that component. Rule-level cutlines are
#' always derived from \code{y} and are not affected by \code{theta_cut}.
#'
#' The student is classified as mastering component \eqn{m} (\code{domain
#' = 1}) if \eqn{\theta_{m}^{EAP} \geq \gamma_m}, and as non-mastering
#' (\code{domain = 0}) otherwise; the same comparison of \eqn{\theta_m^{EAP}}
#' against \eqn{\tau_{km}} classifies mastery of rule \eqn{k}
#' (\code{domain_rules}). Decision confidence is the proportion of
#' posterior draws of \eqn{\theta_m} for that student that agree with the
#' classification: the proportion \eqn{\geq \gamma_m} (or \eqn{\tau_{km}})
#' when classified as mastery, or the proportion \eqn{< \gamma_m} (or
#' \eqn{\tau_{km}}) when classified as non-mastery (Embretson, 2019, Sect.
#' 9.4.2).
#'
#' The difficulty of the rules assigned to each component is read from
#' \code{fit$EAP$eta} (rule-level difficulty, on the \eqn{\theta} scale),
#' and the difficulty of the items that load on each component is read
#' from \code{fit$EAP$beta} (item-level difficulty). \code{Q} and
#' \code{components} are required because \code{fit} does not store the
#' item/component structure matrix \code{C}; it is recomputed internally
#' the same way as in \code{\link{plot_ICC_individual}}.
#'
#' @return A list of class \code{"GMLTM_student_report"} with elements:
#' \describe{
#'   \item{\code{student}}{The resolved identifier of the examinee (as it
#'     appears in \code{rownames(fit$data)}).}
#'   \item{\code{theta}}{A list with \code{EAP}, \code{lower}, \code{upper}
#'     (named by component) and \code{credible}.}
#'   \item{\code{cutline}}{A list with \code{y}, \code{gamma} (named by
#'     component), \code{tau} (named list by component of named numeric
#'     vectors of \eqn{\tau_{km}} by rule, excluding any rule omitted for
#'     lack of items) and \code{fixed} (logical, \code{TRUE} if
#'     \code{theta_cut} was supplied).}
#'   \item{\code{domain}}{Named binary vector (\code{1} = mastery, \code{0}
#'     = non-mastery) per component.}
#'   \item{\code{domain_rules}}{Named list by component of named binary
#'     vectors (\code{1} = mastery, \code{0} = non-mastery) per rule,
#'     excluding any rule omitted for lack of items.}
#'   \item{\code{confidence}}{Named decision-confidence index per
#'     component.}
#'   \item{\code{confidence_rules}}{Named list by component of named
#'     decision-confidence indices per rule, same exclusions as
#'     \code{domain_rules}.}
#'   \item{\code{difficulty}}{A list with \code{rules} (named list of
#'     rule-level \eqn{\eta} difficulties per component) and \code{items}
#'     (named list of item-level \eqn{\beta} difficulties per component).}
#'   \item{\code{plots}}{Named list of \code{ggplot2} objects: one per
#'     component, showing the rule difficulty continuum, the mastery
#'     cutline, and the student's \eqn{\theta} estimate with its credible
#'     interval (colored green for mastery, red for non-mastery), inspired
#'     by Figs. 4 and 5 of Embretson & Yang (2013); plus \code{badges}, a
#'     single scannable summary of every component and rule with a green/red
#'     mastery badge, complementing the continuum plots. \code{NULL} if
#'     \code{include_plots = FALSE}.}
#' }
#'
#' @references
#' Embretson, S. E. (2019). Diagnostic modeling of skill hierarchies and
#' cognitive processes with MLTM-D. In M. von Davier & Y.-S. Lee (Eds.),
#' \emph{Handbook of Diagnostic Classification Models} (pp. 185--208).
#' Springer. \doi{10.1007/978-3-030-05584-4_9}
#'
#' Embretson, S. E., & Yang, X. (2013). A multicomponent latent trait model
#' for diagnosis. \emph{Psychometrika}, \bold{78}, 14--36.
#' \doi{10.1007/s11336-012-9296-y}
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
#'   fit <- GMLTM(analogy, Q, components, iters = 200, iter_warmup = 100, chains = 1, cores = 1)
#'   report <- student_report(fit, Q, components, student_id = 11)
#'   report$domain
#'   report$domain_rules
#'   report$plots$global
#'   report$plots$badges
#' }
#'
#' @family diagnostic reporting functions
#' @export
student_report <- function(fit, Q, components, student_id,
                            y = 0.50, theta_cut = NULL,
                            credible = 0.95, digits = 3,
                            include_plots = TRUE) {

  if (!inherits(fit, "GMLTM"))
    stop("'fit' must be an object of class 'GMLTM', as returned by GMLTM() or GMLTM_corr().")
  if (!is.numeric(y) || length(y) != 1 || y <= 0 || y >= 1)
    stop("'y' must be a single number strictly between 0 and 1.")
  if (!is.numeric(credible) || length(credible) != 1 || credible <= 0 || credible >= 1)
    stop("'credible' must be a single number strictly between 0 and 1.")

  M <- length(components)
  comp_names <- colnames(fit$EAP$theta)
  if (M != ncol(fit$EAP$theta))
    stop("'components' must have the same length as the number of components in 'fit'.")
  if (ncol(Q) != nrow(fit$EAP$eta))
    stop("'Q' must have the same number of columns (rules) as 'fit$EAP$eta'.")
  if (nrow(Q) != nrow(fit$EAP$alpha))
    stop("'Q' must have the same number of rows (items) as 'fit$EAP$alpha'.")

  if (!is.null(theta_cut)) {
    if (!is.numeric(theta_cut) || !(length(theta_cut) %in% c(1, M)))
      stop("'theta_cut' must be a single number or a numeric vector of length length(components).")
    if (length(theta_cut) == 1) {
      gamma_m <- rep(theta_cut, M)
    } else if (!is.null(names(theta_cut))) {
      gamma_m <- theta_cut[comp_names]
    } else {
      gamma_m <- theta_cut
    }
    names(gamma_m) <- comp_names
  }

  subject_names <- rownames(fit$EAP$theta)
  if (is.numeric(student_id)) {
    idx <- as.integer(student_id)
    if (length(idx) != 1 || idx < 1 || idx > nrow(fit$EAP$theta))
      stop("'student_id' as a numeric index must be a single integer in 1:nrow(fit$data).")
  } else {
    idx <- match(as.character(student_id), subject_names)
    if (is.na(idx))
      stop(sprintf("'student_id' = '%s' was not found in rownames(fit$data).", student_id))
  }
  student_label <- subject_names[idx]

  C <- get_C(Q, components)

  theta_EAP_student <- fit$EAP$theta[idx, ]
  theta_draws <- matrix(fit$posterior$theta[, idx, ], ncol = M)

  alpha_level <- 1 - credible
  probs <- c(alpha_level / 2, 1 - alpha_level / 2)
  ci <- apply(theta_draws, MARGIN = 2, FUN = quantile, probs = probs)
  theta_lower <- ci[1, ]
  theta_upper <- ci[2, ]
  names(theta_EAP_student) <- names(theta_lower) <- names(theta_upper) <- comp_names

  if (is.null(theta_cut)) {
    gamma_m <- vapply(seq_len(M), function(m) {
      items_m <- which(C[, m] == 1)
      .gmltm_component_gamma(
        y = y,
        alpha_m = fit$EAP$alpha[items_m, m],
        beta_m  = fit$EAP$beta[items_m, m],
        guess_m = fit$EAP$guessing[items_m],
        comp_name = comp_names[m]
      )
    }, numeric(1))
    names(gamma_m) <- comp_names
  }

  mastery <- as.integer(theta_EAP_student >= gamma_m)
  names(mastery) <- comp_names

  confidence_idx <- numeric(M)
  for (m in seq_len(M)) {
    draws_m <- theta_draws[, m]
    confidence_idx[m] <- if (mastery[m] == 1) mean(draws_m >= gamma_m[m]) else mean(draws_m < gamma_m[m])
  }
  names(confidence_idx) <- comp_names

  # Rule-level cutlines tau_km, mastery, and decision confidence: same
  # procedure as gamma_m/mastery/confidence_idx above, restricted to the
  # items requiring each specific rule within its component. Rules with no
  # items of their own in the component are omitted (with a warning).
  rule_names_all <- rownames(fit$EAP$eta)
  tau_m <- vector("list", M)
  domain_rules <- vector("list", M)
  confidence_rules <- vector("list", M)
  names(tau_m) <- names(domain_rules) <- names(confidence_rules) <- comp_names

  for (m in seq_len(M)) {
    rules_m <- components[[m]]
    draws_m <- theta_draws[, m]

    tau_vals <- numeric(0)
    dom_vals <- integer(0)
    conf_vals <- numeric(0)
    kept_names <- character(0)

    for (k in rules_m) {
      tau_k <- .gmltm_rule_tau(
        y = y, fit = fit, Q = Q, C = C,
        comp_idx = m, rule_idx = k,
        rule_name = rule_names_all[k]
      )
      if (is.na(tau_k)) next

      dom_k  <- as.integer(theta_EAP_student[m] >= tau_k)
      conf_k <- if (dom_k == 1) mean(draws_m >= tau_k) else mean(draws_m < tau_k)

      tau_vals   <- c(tau_vals, tau_k)
      dom_vals   <- c(dom_vals, dom_k)
      conf_vals  <- c(conf_vals, conf_k)
      kept_names <- c(kept_names, rule_names_all[k])
    }

    names(tau_vals) <- names(dom_vals) <- names(conf_vals) <- kept_names
    tau_m[[m]] <- tau_vals
    domain_rules[[m]] <- dom_vals
    confidence_rules[[m]] <- conf_vals
  }

  rule_difficulty <- lapply(seq_len(M), function(m) {
    rules_m <- components[[m]]
    vals <- fit$EAP$eta[rules_m, m]
    names(vals) <- rownames(fit$EAP$eta)[rules_m]
    vals
  })
  names(rule_difficulty) <- comp_names

  item_difficulty <- lapply(seq_len(M), function(m) {
    items_m <- which(C[, m] == 1)
    vals <- fit$EAP$beta[items_m, m]
    names(vals) <- rownames(fit$EAP$beta)[items_m]
    vals
  })
  names(item_difficulty) <- comp_names

  plots <- NULL
  if (include_plots) {
    plots <- lapply(seq_len(M), function(m) {
      .gmltm_plot_component(
        comp_name = comp_names[m],
        rule_diff = rule_difficulty[[m]],
        theta_eap = theta_EAP_student[m],
        theta_lo  = theta_lower[m],
        theta_hi  = theta_upper[m],
        gamma     = gamma_m[m],
        mastery   = mastery[m],
        y = y, digits = digits
      )
    })
    names(plots) <- comp_names

    plots$badges <- .gmltm_plot_badges(
      comp_names        = comp_names,
      theta_eap          = theta_EAP_student,
      domain             = mastery,
      tau_list           = tau_m,
      domain_rules_list  = domain_rules,
      digits             = digits
    )
  }

  result <- list(
    student = student_label,
    theta = list(EAP = theta_EAP_student, lower = theta_lower, upper = theta_upper,
                credible = credible),
    cutline = list(y = y, gamma = gamma_m, tau = tau_m, fixed = !is.null(theta_cut)),
    domain = mastery,
    domain_rules = domain_rules,
    confidence = confidence_idx,
    confidence_rules = confidence_rules,
    difficulty = list(rules = rule_difficulty, items = item_difficulty),
    plots = plots
  )
  class(result) <- "GMLTM_student_report"
  return(result)
}

#' Deprecated Alias for \code{student_report}
#'
#' @description
#' \code{informe_estudiante} is a deprecated alias for
#' \code{\link{student_report}}, kept for backward compatibility with GMLTM
#' (< 2.0.0). It will be removed in a future release; use
#' \code{\link{student_report}} instead.
#'
#' @param ... Arguments passed on to \code{\link{student_report}}.
#'
#' @return See \code{\link{student_report}}.
#' @keywords internal
#' @export
informe_estudiante <- function(...) {
  .Deprecated("student_report")
  student_report(...)
}

#' Print Method for GMLTM_student_report Objects
#'
#' @param x An object of class \code{"GMLTM_student_report"}, as
#'   returned by \code{\link{student_report}}.
#' @param digits Number of decimal digits to display. Default is 3.
#' @param ... Currently ignored.
#'
#' @return Invisibly returns \code{x}.
#' @export
print.GMLTM_student_report <- function(x, digits = 3, ...) {
  cat(sprintf("Individual GMLTM-D report - Student: %s\n", x$student))
  cat("--------------------------------------------------------\n")
  cat(sprintf("Mastery cutline: y = %.2f%s\n", x$cutline$y,
             if (isTRUE(x$cutline$fixed)) " (gamma set manually via theta_cut)" else ""))
  cat("\n")

  comp_names <- names(x$domain)
  for (m in seq_along(comp_names)) {
    cm <- comp_names[m]
    state_lab <- if (x$domain[m] == 1) "MASTERY" else "NON-MASTERY"
    gamma_val <- x$cutline$gamma[m]
    gamma_lab <- if (is.finite(gamma_val)) sprintf("%.*f", digits, gamma_val) else as.character(gamma_val)
    cat(sprintf(
      "%s: theta = %.*f [%.*f, %.*f]  |  gamma = %s  |  %s (confidence = %.*f)\n",
      cm, digits, x$theta$EAP[m], digits, x$theta$lower[m], digits, x$theta$upper[m],
      gamma_lab, state_lab, digits, x$confidence[m]
    ))

    rule_names_m <- names(x$domain_rules[[cm]])
    for (k in seq_along(rule_names_m)) {
      rule_state <- if (x$domain_rules[[cm]][k] == 1) "mastery" else "non-mastery"
      cat(sprintf(
        "    - %s: tau = %.*f  |  %s (confidence = %.*f)\n",
        rule_names_m[k], digits, x$cutline$tau[[cm]][k],
        rule_state, digits, x$confidence_rules[[cm]][k]
      ))
    }
  }
  cat("\n")
  invisible(x)
}
