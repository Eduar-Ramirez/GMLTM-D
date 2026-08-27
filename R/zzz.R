.onAttach <- function(libname, pkgname) {
  if (requireNamespace("rstan", quietly = TRUE)) {
    rstan::rstan_options(auto_write = TRUE)
  }
}
