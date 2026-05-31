## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.
* All Stan models are compiled at runtime via rstan::stan_model().
  Compilation occurs on first use per session; subsequent calls reuse
  the cached compiled model due to rstan::rstan_options(auto_write = TRUE).
* Examples with Stan sampling are wrapped in \donttest{} due to
  compilation and sampling time (typically 30-120 seconds).
* The SystemRequirements field lists C++17 and GNU make, required by rstan.
