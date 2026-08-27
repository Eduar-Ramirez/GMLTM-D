test_that("domain_rules has one entry per model rule when no rule is omitted", {
  fit <- make_mock_fit(Q_normal, components_normal)
  report <- student_report(fit, Q_normal, components_normal, student_id = 1)

  expect_equal(length(unlist(report$domain_rules)), nrow(fit$EAP$eta))
  expect_equal(length(unlist(report$confidence_rules)), nrow(fit$EAP$eta))
  expect_setequal(names(report$domain_rules), names(components_normal))
})

test_that("domain_rules and confidence_rules hold valid values", {
  fit <- make_mock_fit(Q_normal, components_normal)
  report <- student_report(fit, Q_normal, components_normal, student_id = 2)

  domain_vals <- unlist(report$domain_rules)
  expect_true(all(domain_vals %in% c(0L, 1L)))

  conf_vals <- unlist(report$confidence_rules)
  expect_true(all(conf_vals >= 0 & conf_vals <= 1))

  expect_true(all(is.finite(unlist(report$cutline$tau))))
})

test_that("a rule with no items of its own in its component is omitted, not an error", {
  fit <- make_mock_fit(Q_edge, components_edge)

  expect_warning(
    report <- student_report(fit, Q_edge, components_edge, student_id = 1),
    "no items in this component require this rule"
  )

  # rule2 (component "global") has no items of its own and must be dropped,
  # not merely set to NA.
  expect_equal(length(report$domain_rules$global), 1)
  expect_false("rule2" %in% names(report$domain_rules$global))
  expect_true("rule1" %in% names(report$domain_rules$global))

  # rule3 (component "local") is unaffected by the omission of rule2.
  expect_equal(length(report$domain_rules$local), 1)
  expect_true("rule3" %in% names(report$domain_rules$local))

  expect_true(all(unlist(report$domain_rules) %in% c(0L, 1L)))
})

test_that("plots$badges is produced alongside the per-component continuum plots", {
  fit <- make_mock_fit(Q_normal, components_normal)
  report <- student_report(fit, Q_normal, components_normal, student_id = 1)

  expect_true("badges" %in% names(report$plots))
  expect_s3_class(report$plots$badges, "ggplot")
  expect_s3_class(report$plots$global, "ggplot")
})

test_that("print.GMLTM_student_report runs without error on the new fields", {
  fit <- make_mock_fit(Q_normal, components_normal)
  report <- student_report(fit, Q_normal, components_normal, student_id = 1)
  expect_output(print(report), "MASTERY|NON-MASTERY")
})
