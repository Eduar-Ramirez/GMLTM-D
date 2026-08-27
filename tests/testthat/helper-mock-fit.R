# Mock "GMLTM"-classed fit object with the minimal fields student_report()
# and student_report_batch() read (fit$data, fit$EAP$theta/eta/alpha/beta/
# guessing, fit$posterior$theta), so tests run instantly instead of
# requiring a real Stan fit. Shared across test-student_report*.R files
# (testthat auto-sources helper-*.R before running any tests).
make_mock_fit <- function(Q, components, n = 5, iters = 200, seed = 1) {
  set.seed(seed)
  p <- nrow(Q)
  K <- ncol(Q)
  M <- length(components)
  comp_names <- names(components)

  alpha <- matrix(stats::runif(p * M, 0.8, 1.2), nrow = p, ncol = M,
                   dimnames = list(paste0("item", seq_len(p)), comp_names))
  beta <- matrix(stats::rnorm(p * M, 0, 0.5), nrow = p, ncol = M,
                 dimnames = list(paste0("item", seq_len(p)), comp_names))
  guessing <- stats::setNames(rep(0.1, p), paste0("item", seq_len(p)))

  eta <- matrix(stats::rnorm(K * M, 0, 0.5), nrow = K, ncol = M,
               dimnames = list(paste0("rule", seq_len(K)), comp_names))

  theta_eap <- matrix(stats::rnorm(n * M, 0, 1), nrow = n, ncol = M,
                      dimnames = list(paste0("S", seq_len(n)), comp_names))

  theta_draws <- array(0, dim = c(iters, n, M))
  for (i in seq_len(n)) {
    for (m in seq_len(M)) {
      theta_draws[, i, m] <- stats::rnorm(iters, mean = theta_eap[i, m], sd = 0.3)
    }
  }

  data <- matrix(NA_integer_, nrow = n, ncol = p,
                 dimnames = list(paste0("S", seq_len(n)), paste0("item", seq_len(p))))

  fit <- list(
    data = data,
    EAP = list(theta = theta_eap, eta = eta, alpha = alpha, beta = beta, guessing = guessing),
    posterior = list(theta = theta_draws)
  )
  class(fit) <- "GMLTM"
  fit
}

# Normal case: every rule has at least one item of its own within its
# component (no omission expected).
Q_normal <- matrix(
  c(1, 0, 0,
    0, 1, 0,
    1, 1, 0,
    0, 0, 1),
  nrow = 4, ncol = 3, byrow = TRUE
)
components_normal <- list(global = c(1, 2), local = c(3))

# Edge case: rule 2 (assigned to "global") never appears in Q at all, so it
# has no items of its own within the "global" component.
Q_edge <- matrix(
  c(1, 0, 0,
    1, 0, 0,
    0, 0, 1,
    0, 0, 1),
  nrow = 4, ncol = 3, byrow = TRUE
)
components_edge <- list(global = c(1, 2), local = c(3))
