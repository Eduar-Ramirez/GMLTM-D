test_that("LLTM rejects invalid priors (negative sigma)", {
  expect_error(
    LLTM(matrix(1, 5, 3), matrix(1, 3, 2),
         priors = list(theta = list(mu = 0, sigma = -1))),
    "sigma"
  )
})

test_that("LLTM rejects invalid priors (zero sigma)", {
  expect_error(
    LLTM(matrix(1, 5, 3), matrix(1, 3, 2),
         priors = list(eta = list(mu = 0, sigma = 0))),
    "sigma"
  )
})

test_that("modifyList preserves unspecified prior defaults", {
  default <- list(theta = list(mu = 0, sigma = 1), eta = list(mu = 0, sigma = 1))
  custom  <- list(theta = list(sigma = 2))
  result  <- utils::modifyList(default, custom)
  expect_equal(result$eta$sigma, 1)
  expect_equal(result$theta$sigma, 2)
  expect_equal(result$theta$mu, 0)
})

test_that("GMLTM rejects invalid Beta prior shapes", {
  Q <- matrix(1, 3, 2)
  components <- list(c1 = 1, c2 = 2)
  expect_error(
    GMLTM(matrix(1, 5, 3), Q, components,
          priors = list(c = list(shape1 = -1, shape2 = 20))),
    "shape1"
  )
  expect_error(
    GMLTM(matrix(1, 5, 3), Q, components,
          priors = list(c = list(shape1 = 3, shape2 = 0))),
    "shape2"
  )
})

test_that("GMLTM fits successfully with priors$alpha$family = 'lognormal'", {
  Q <- matrix(c(1, 0, 0, 1, 1, 1, 0, 1), nrow = 4, ncol = 2, byrow = TRUE)
  components <- list(a = 1, b = 2)
  data <- as.data.frame(matrix(stats::rbinom(6 * 4, 1, 0.5), 6, 4))

  fit <- suppressWarnings(GMLTM(
    data, Q, components, iters = 20, iter_warmup = 20,
    chains = 1, cores = 1,
    priors = list(alpha = list(mu = 0, sigma = 1, family = "lognormal"))
  ))

  expect_s3_class(fit, "GMLTM")
  expect_equal(fit$priors$alpha$family, "lognormal")
})

test_that("GMLTM rejects an invalid priors$alpha$family", {
  Q <- matrix(1, 3, 2)
  components <- list(c1 = 1, c2 = 2)
  expect_error(
    GMLTM(matrix(1, 5, 3), Q, components,
          priors = list(alpha = list(family = "uniform"))),
    'must be either "normal" or "lognormal"'
  )
})

test_that("the data_list passed to rstan carries the correct alpha_lognormal flag", {
  captured <- NULL
  testthat::local_mocked_bindings(
    stan_model = function(...) structure(list(), class = "stanmodel"),
    sampling = function(object, data, ...) {
      captured <<- data$alpha_lognormal
      stop("mock: captured data_list")
    },
    .package = "rstan"
  )

  Q <- matrix(c(1, 0, 0, 1, 1, 1, 0, 1), nrow = 4, ncol = 2, byrow = TRUE)
  components <- list(a = 1, b = 2)
  data <- matrix(0, 6, 4)

  expect_error(
    GMLTM(data, Q, components, priors = list(alpha = list(family = "lognormal"))),
    "mock: captured data_list"
  )
  expect_equal(captured, 1L)

  expect_error(
    GMLTM(data, Q, components, priors = list(alpha = list(family = "normal"))),
    "mock: captured data_list"
  )
  expect_equal(captured, 0L)

  expect_error(GMLTM(data, Q, components), "mock: captured data_list")
  expect_equal(captured, 0L)
})
