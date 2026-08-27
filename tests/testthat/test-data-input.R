# Regression test: GMLTM()/GMLTM_corr()/MLTM()/LLTM() used to build the
# Stan `y` vector via `unname(unlist(data))`, which only flattens correctly
# for a data.frame -- unlist() on a plain matrix returns it unchanged
# (dim preserved), so Stan received a matrix where it expected a flat
# vector and failed. Fixed via as.vector(as.matrix(data)), which flattens
# both input types identically (column-major, matching `item`/`ID`).

.capture_y <- function(fit_call) {
  captured <- NULL
  testthat::local_mocked_bindings(
    stan_model = function(...) structure(list(), class = "stanmodel"),
    sampling = function(object, data, ...) {
      captured <<- data$y
      stop("mock: captured data_list")
    },
    .package = "rstan"
  )
  testthat::expect_error(fit_call(), "mock: captured data_list")
  captured
}

test_that("a plain matrix and an equivalent data.frame produce the same flattened y", {
  Q <- matrix(c(1, 0, 0, 1, 1, 1, 0, 1), nrow = 4, ncol = 2, byrow = TRUE)
  components <- list(a = 1, b = 2)
  raw <- matrix(c(1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1),
                nrow = 6, ncol = 4)
  data_matrix <- raw
  data_frame  <- as.data.frame(raw)

  y_matrix <- .capture_y(function() GMLTM(data_matrix, Q, components))
  y_frame  <- .capture_y(function() GMLTM(data_frame, Q, components))

  expect_null(dim(y_matrix))
  expect_equal(y_matrix, y_frame)
  expect_equal(y_matrix, as.vector(raw))
})

test_that("GMLTM_corr/MLTM/LLTM also flatten a plain matrix correctly", {
  Q <- matrix(c(1, 0, 0, 1, 1, 1, 0, 1), nrow = 4, ncol = 2, byrow = TRUE)
  components <- list(a = 1, b = 2)
  raw <- matrix(c(1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1),
                nrow = 6, ncol = 4)

  y_corr  <- .capture_y(function() GMLTM_corr(raw, Q, components))
  y_mltm  <- .capture_y(function() MLTM(raw, Q, components))
  y_lltm  <- .capture_y(function() LLTM(raw, Q))

  expect_equal(y_corr, as.vector(raw))
  expect_equal(y_mltm, as.vector(raw))
  expect_equal(y_lltm, as.vector(raw))
})

test_that("fitting with a plain matrix runs end to end without error (all four models)", {
  set.seed(1)
  p <- 4; K <- 2; n <- 6
  Q <- matrix(c(1, 0, 0, 1, 1, 1, 0, 1), nrow = p, ncol = K, byrow = TRUE)
  components <- list(a = 1, b = 2)
  data_matrix <- matrix(stats::rbinom(n * p, 1, 0.6), n, p)

  expect_s3_class(
    suppressWarnings(GMLTM(data_matrix, Q, components, iters = 20, iter_warmup = 20,
                           chains = 1, cores = 1)),
    "GMLTM"
  )
  expect_s3_class(
    suppressWarnings(GMLTM_corr(data_matrix, Q, components, iters = 20, iter_warmup = 20,
                                chains = 1, cores = 1)),
    "GMLTM_corr"
  )
  expect_s3_class(
    suppressWarnings(MLTM(data_matrix, Q, components, iters = 20, iter_warmup = 20,
                          chains = 1, cores = 1)),
    "MLTM"
  )
  expect_s3_class(
    suppressWarnings(LLTM(data_matrix, Q, iters = 20, iter_warmup = 20,
                          chains = 1, cores = 1)),
    "LLTM"
  )
})
