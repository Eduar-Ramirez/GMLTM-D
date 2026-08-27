ppc_Q <- matrix(
  c(1, 0, 0,
    0, 1, 0,
    1, 1, 0,
    0, 0, 1,
    1, 0, 1,
    0, 1, 1),
  nrow = 6, ncol = 3, byrow = TRUE
)
ppc_components <- list(global = c(1, 2), local = c(3))

test_that("simulated proportions are always between 0 and 1", {
  ppc <- prior_predictive_check(ppc_Q, ppc_components, N = 40, S = 50)

  expect_true(all(ppc$global >= 0 & ppc$global <= 1))
  expect_true(all(ppc$by_item >= 0 & ppc$by_item <= 1))
  expect_true(all(ppc$by_component >= 0 & ppc$by_component <= 1, na.rm = TRUE))
})

test_that("output dimensions match S, N and the Q/components structure", {
  S <- 50
  ppc <- prior_predictive_check(ppc_Q, ppc_components, N = 40, S = S)

  expect_length(ppc$global, S)
  expect_equal(dim(ppc$by_item), c(S, nrow(ppc_Q)))
  expect_equal(dim(ppc$by_component), c(S, length(ppc_components)))
  expect_equal(colnames(ppc$by_component), names(ppc_components))
  expect_equal(ppc$N, 40L)
  expect_equal(ppc$S, S)
})

test_that("by_component is NULL when by_component = FALSE", {
  ppc <- prior_predictive_check(ppc_Q, ppc_components, N = 40, S = 50, by_component = FALSE)
  expect_null(ppc$by_component)
})

test_that("changing priors$alpha$family perceptibly changes the simulated distribution", {
  set.seed(2026)
  ppc_normal <- prior_predictive_check(
    ppc_Q, ppc_components, N = 200, S = 500,
    priors = list(alpha = list(mu = 0, sigma = 1, family = "normal"))
  )
  ppc_lognormal <- prior_predictive_check(
    ppc_Q, ppc_components, N = 200, S = 500,
    priors = list(alpha = list(mu = 0, sigma = 1, family = "lognormal"))
  )

  # Same sigma, only the alpha family differs: the Log-Normal(0, 1) prior on
  # alpha has a visibly heavier right tail than the half-Normal(0, 1), which
  # should inflate the spread of simulated global proportions perceptibly
  # (checked across several seeds during development: consistently >40%
  # higher variance, never close to 1).
  var_normal    <- var(ppc_normal$global)
  var_lognormal <- var(ppc_lognormal$global)
  expect_gt(var_lognormal, 1.2 * var_normal)
})

test_that("prior_predictive_check validates its own priors like GMLTM()", {
  expect_error(
    prior_predictive_check(ppc_Q, ppc_components, N = 40,
                            priors = list(alpha = list(family = "uniform"))),
    'must be either "normal" or "lognormal"'
  )
  expect_error(
    prior_predictive_check(ppc_Q, ppc_components, N = -1),
    "positive integer"
  )
  expect_error(
    prior_predictive_check(ppc_Q, ppc_components, N = 40, S = 0),
    "positive integer"
  )
})

test_that("plot_prior_predictive_check returns a ggplot for each type", {
  ppc <- prior_predictive_check(ppc_Q, ppc_components, N = 40, S = 50)

  expect_s3_class(plot_prior_predictive_check(ppc), "ggplot")
  expect_s3_class(plot_prior_predictive_check(ppc, type = "by_item"), "ggplot")
  expect_s3_class(plot_prior_predictive_check(ppc, type = "by_component"), "ggplot")

  ppc_no_comp <- prior_predictive_check(ppc_Q, ppc_components, N = 40, S = 50,
                                         by_component = FALSE)
  expect_error(plot_prior_predictive_check(ppc_no_comp, type = "by_component"),
               "by_component = FALSE")
})
