# =====================================================================
# student_report_batch(): paginated batch reporting on top of
# student_report(). Every requested student gets the cheap numeric
# summary (theta EAP, component/rule mastery, confidence -- all already
# available from fit$EAP/fit$posterior$theta), but the expensive part
# (building the ggplot2 objects in student_report()$plots) only runs for
# the students in the current page (student_ids[(offset+1):(offset+limit)]).
# =====================================================================

# One-row data.frame with the lightweight numeric fields of a single
# student_report() result (theta EAP / domain / confidence per component,
# and per rule). Column names are stable across students because they all
# come from the same Q/components, so plain rbind() is safe.
.gmltm_report_to_row <- function(report) {
  comp_names <- names(report$domain)
  row <- data.frame(student = report$student, stringsAsFactors = FALSE)

  for (cm in comp_names) {
    row[[paste0("theta_", cm)]]      <- unname(report$theta$EAP[cm])
    row[[paste0("domain_", cm)]]     <- unname(report$domain[cm])
    row[[paste0("confidence_", cm)]] <- unname(report$confidence[cm])
  }

  for (cm in comp_names) {
    rule_names <- names(report$domain_rules[[cm]])
    for (rk in rule_names) {
      row[[paste0("domain_rule_", rk)]]     <- unname(report$domain_rules[[cm]][rk])
      row[[paste0("confidence_rule_", rk)]] <- unname(report$confidence_rules[[cm]][rk])
    }
  }

  row
}

#' Paginated Batch Diagnostic Reports for GMLTM-D Models
#'
#' @description
#' Runs \code{\link{student_report}} over many examinees at once. A
#' lightweight numeric summary (posterior EAP of \eqn{\theta}, mastery and
#' decision confidence per component and per rule) is always computed for
#' every requested student, since it only reads quantities already stored in
#' \code{fit$EAP}/\code{fit$posterior$theta}. The comparatively expensive
#' step -- building the \code{ggplot2} objects returned by
#' \code{\link{student_report}} -- is only performed for the students in the
#' current page, so large batches stay cheap to page through.
#'
#' @param fit An object of class \code{"GMLTM"} (as returned by
#'   \code{\link{GMLTM}} or \code{\link{GMLTM_corr}}).
#' @param Q A \eqn{p \times K} Q-matrix, identical to the one used to fit
#'   \code{fit}, specifying which rules each item requires.
#' @param components A named list grouping rules into components, identical
#'   to the one used to fit \code{fit}.
#' @param student_ids Vector of examinee identifiers to include in the
#'   batch, in either form accepted by \code{\link{student_report}}'s
#'   \code{student_id} (row names of \code{fit$data}, or integer row
#'   indices). Default is \code{NULL}, meaning every student in
#'   \code{fit$data} (\code{seq_len(nrow(fit$data))}).
#' @param offset Integer. Number of students in \code{student_ids} to skip
#'   before the current page. Default is \code{0}.
#' @param limit Integer. Maximum number of students in the current page
#'   (i.e. with full reports/plots generated). Default is \code{10}.
#' @param y,theta_cut,credible,digits Passed through to
#'   \code{\link{student_report}}; see its documentation.
#'
#' @details
#' The current page covers
#' \code{student_ids[(offset + 1):(offset + limit)]}. If \code{offset} is at
#' or beyond \code{length(student_ids)}, the page is empty (\code{reports}
#' has length 0) rather than an error; \code{summary} is unaffected, since it
#' always covers every requested student regardless of paging.
#'
#' @return A list of class \code{"GMLTM_batch_report"} with elements:
#' \describe{
#'   \item{\code{summary}}{A \code{data.frame} with one row per requested
#'     student (\code{student}, \code{theta_<component>},
#'     \code{domain_<component>}, \code{confidence_<component>} for each
#'     component, and \code{domain_rule_<rule>}, \code{confidence_rule_<rule>}
#'     for each rule), built without generating any plots.}
#'   \item{\code{reports}}{Named list (by resolved student label) of full
#'     \code{\link{student_report}} objects, including their \code{plots},
#'     for the students in the current page only.}
#'   \item{\code{page}}{A list with \code{offset}, \code{limit}, and
#'     \code{total} (the number of requested students in \code{student_ids}).}
#' }
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
#'
#'   # First page: students 1-10 (summary covers all 149 examinees either way)
#'   page1 <- student_report_batch(fit, Q, components, limit = 10)
#'   page1
#'   page1$reports[["1"]]$plots$badges
#'
#'   # Next page: students 11-20, individual reports for students 1-10 are
#'   # still fully available in page1$reports
#'   page2 <- student_report_batch(fit, Q, components, offset = 10, limit = 10)
#'   page2$reports[["11"]]
#' }
#'
#' @family diagnostic reporting functions
#' @export
student_report_batch <- function(fit, Q, components, student_ids = NULL,
                                  offset = 0, limit = 10,
                                  y = 0.50, theta_cut = NULL,
                                  credible = 0.95, digits = 3) {

  if (!inherits(fit, "GMLTM"))
    stop("'fit' must be an object of class 'GMLTM', as returned by GMLTM() or GMLTM_corr().")
  if (is.null(student_ids))
    student_ids <- seq_len(nrow(fit$data))
  if (length(student_ids) == 0)
    stop("'student_ids' must have at least one element.")
  if (!is.numeric(offset) || length(offset) != 1 || offset < 0 || offset != as.integer(offset))
    stop("'offset' must be a single non-negative integer.")
  if (!is.numeric(limit) || length(limit) != 1 || limit <= 0 || limit != as.integer(limit))
    stop("'limit' must be a single positive integer.")
  offset <- as.integer(offset)
  limit  <- as.integer(limit)

  n_total <- length(student_ids)
  window_idx <- seq.int(offset + 1L, offset + limit)
  window_idx <- window_idx[window_idx >= 1L & window_idx <= n_total]
  window_ids <- student_ids[window_idx]

  reports <- vector("list", length(window_ids))
  for (i in seq_along(window_ids)) {
    reports[[i]] <- student_report(
      fit = fit, Q = Q, components = components, student_id = window_ids[i],
      y = y, theta_cut = theta_cut, credible = credible, digits = digits,
      include_plots = TRUE
    )
  }
  names(reports) <- vapply(reports, function(r) r$student, character(1))

  summary_rows <- vector("list", n_total)
  for (i in seq_len(n_total)) {
    match_pos <- match(student_ids[i], window_ids)
    report_i <- if (!is.na(match_pos)) {
      reports[[match_pos]]
    } else {
      student_report(
        fit = fit, Q = Q, components = components, student_id = student_ids[i],
        y = y, theta_cut = theta_cut, credible = credible, digits = digits,
        include_plots = FALSE
      )
    }
    summary_rows[[i]] <- .gmltm_report_to_row(report_i)
  }
  summary_df <- do.call(rbind, summary_rows)
  rownames(summary_df) <- NULL

  result <- list(
    summary = summary_df,
    reports = reports,
    page = list(offset = offset, limit = limit, total = n_total)
  )
  class(result) <- "GMLTM_batch_report"
  result
}

#' Print Method for GMLTM_batch_report Objects
#'
#' @param x An object of class \code{"GMLTM_batch_report"}, as returned by
#'   \code{\link{student_report_batch}}.
#' @param digits Number of decimal digits to display. Default is 3.
#' @param ... Currently ignored.
#'
#' @return Invisibly returns \code{x}.
#' @export
print.GMLTM_batch_report <- function(x, digits = 3, ...) {
  cat(sprintf("GMLTM-D batch report -- %d student(s)\n", x$page$total))
  cat("--------------------------------------------------------\n")
  print(x$summary, digits = digits)
  cat("\n")

  n_page <- length(x$reports)
  if (n_page == 0) {
    cat(sprintf(
      "No students in this page (offset = %d is at or beyond the last student, total = %d).\n",
      x$page$offset, x$page$total))
  } else {
    first <- x$page$offset + 1L
    last  <- x$page$offset + n_page
    cat(sprintf(
      "Full reports (with plots) generated for students %d-%d of %d.\n",
      first, last, x$page$total))
    if (last < x$page$total) {
      next_last <- min(last + x$page$limit, x$page$total)
      cat(sprintf("Use offset = %d to see students %d-%d.\n", last, last + 1L, next_last))
    }
  }
  invisible(x)
}
