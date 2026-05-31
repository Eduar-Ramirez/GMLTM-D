# GMLTM 0.1.0

* `GMLTM1()` and `GMLTM2()` have been removed. Their functionality is
  fully covered by the `priors` argument of `GMLTM()`. See the vignette
  for equivalent prior specifications.
* Initial release to CRAN.
* Implements LLTM (Fischer, 1973), MLTM-D (Embretson & Yang, 2013),
  and GMLTM-D (Ramirez et al., 2024) via Bayesian HMC with Stan.
* All models support user-defined prior distributions via the `priors` argument.
* Includes reliability estimation, posterior predictive checks, and
  ICC visualization tools.
* Stan backend uses `rstan` for full CRAN compatibility.
