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
