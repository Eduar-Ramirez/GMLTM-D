test_that("summary has one row per student for a large batch (n = 1000)", {
  fit <- make_mock_fit(Q_normal, components_normal, n = 1000, iters = 30)

  batch <- student_report_batch(fit, Q_normal, components_normal, limit = 10)

  expect_s3_class(batch, "GMLTM_batch_report")
  expect_equal(nrow(batch$summary), 1000)
  expect_equal(batch$page$total, 1000)
  expect_true(all(c("student", "theta_global", "domain_global",
                     "confidence_global", "domain_rule_rule1") %in% names(batch$summary)))
})

test_that("reports only contains the students in the current page (limit)", {
  fit <- make_mock_fit(Q_normal, components_normal, n = 1000, iters = 30)

  batch <- student_report_batch(fit, Q_normal, components_normal, offset = 0, limit = 10)
  expect_equal(length(batch$reports), 10)
  expect_setequal(names(batch$reports), paste0("S", 1:10))

  # Every report in the page is a full student_report(), plots included
  expect_s3_class(batch$reports[["S1"]], "GMLTM_student_report")
  expect_true(is.list(batch$reports[["S1"]]$plots))
  expect_s3_class(batch$reports[["S1"]]$plots$badges, "ggplot")

  batch2 <- student_report_batch(fit, Q_normal, components_normal, offset = 10, limit = 15)
  expect_equal(length(batch2$reports), 15)
  expect_setequal(names(batch2$reports), paste0("S", 11:25))
  # But the summary still covers every student regardless of the page
  expect_equal(nrow(batch2$summary), 1000)
})

test_that("an out-of-range offset yields an empty page, not an error", {
  fit <- make_mock_fit(Q_normal, components_normal, n = 20, iters = 30)

  batch <- student_report_batch(fit, Q_normal, components_normal, offset = 20, limit = 10)
  expect_equal(length(batch$reports), 0)
  expect_equal(nrow(batch$summary), 20)
  expect_output(print(batch), "No students in this page")

  # An offset that only partially overflows returns a partial (shorter) page
  batch_partial <- student_report_batch(fit, Q_normal, components_normal, offset = 15, limit = 10)
  expect_equal(length(batch_partial$reports), 5)
  expect_setequal(names(batch_partial$reports), paste0("S", 16:20))
})

test_that("print.GMLTM_batch_report reports the next-page offset", {
  fit <- make_mock_fit(Q_normal, components_normal, n = 30, iters = 30)

  batch <- student_report_batch(fit, Q_normal, components_normal, offset = 0, limit = 10)
  expect_output(print(batch), "Use offset = 10 to see students 11-20")
})

test_that("explicit student_ids are respected and preserved in order", {
  fit <- make_mock_fit(Q_normal, components_normal, n = 20, iters = 30)

  ids <- c(5, 3, 5, 8)
  batch <- student_report_batch(fit, Q_normal, components_normal,
                                 student_ids = ids, offset = 0, limit = 2)
  expect_equal(nrow(batch$summary), length(ids))
  expect_equal(batch$summary$student, paste0("S", ids))
  expect_equal(length(batch$reports), 2)
})

test_that("invalid offset/limit are rejected", {
  fit <- make_mock_fit(Q_normal, components_normal, n = 10, iters = 30)

  expect_error(student_report_batch(fit, Q_normal, components_normal, offset = -1),
               "non-negative integer")
  expect_error(student_report_batch(fit, Q_normal, components_normal, limit = 0),
               "positive integer")
  expect_error(student_report_batch(fit, Q_normal, components_normal, limit = 2.5),
               "positive integer")
})
