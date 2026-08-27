# GMLTM 2.0.0

* New `GMLTM_corr()`: a variant of `GMLTM()` that replaces the independent
  prior on theta (Sigma = sigma^2 * I) with a multivariate Normal prior with
  a general correlation matrix Sigma among cognitive components, estimated
  via a Cholesky factor with an LKJ prior (argument `lkj_eta`). Implements
  innovation (iii) of the GMLTM-D as formulated in Ramirez et al. (2024):
  admitting correlations between components. Returns `EAP$Sigma` and
  `quantiles$Sigma` in addition to the elements already returned by `GMLTM()`.
* New `extract_correlation()`: extracts the posterior correlation matrix
  Sigma from a `GMLTM_corr` fit, reports pairwise credible intervals for each
  pair of components (flagging whether the interval excludes 0), and draws a
  correlation heatmap consistent with the style used in
  `generate_Q_with_interactions()`.
* `compute_model_validation()` now labels models by class in `$Summary` and
  `$Comparison` (e.g. `"GMLTM-D (Sigma libre)"` for `GMLTM_corr` fits,
  `"GMLTM-D (Sigma=I)"` for `GMLTM` fits, `"MLTM-D"`, `"LLTM"`), so `GMLTM`,
  `GMLTM_corr`, `MLTM`, and `LLTM` fits can be freely mixed in the same
  comparison list.
* New `student_report()`: builds a per-component and per-rule mastery report
  for a single examinee (theta EAP and credible interval, mastery cutlines
  `gamma_m`/`tau_km`, decision-confidence indices, and rule/item difficulty),
  following the mastery-diagnosis procedure for the MLTM-D in Embretson
  (2019). Includes a scannable badge-style summary plot
  (`report$plots$badges`) alongside the existing per-component continuum
  plots. This function was previously named `informe_estudiante()`; that
  name is kept as a deprecated alias and will be removed in a future
  release.
* New `student_report_batch()`: pages `student_report()` over many
  examinees. Always computes the lightweight numeric summary (theta EAP,
  component/rule mastery and confidence) for every requested student, but
  only builds the `ggplot2` plots for the students in the current
  `offset`/`limit` page, so large batches stay cheap to page through.
* `priors$alpha` in `GMLTM()`/`GMLTM_corr()` now accepts a `family` element,
  either `"normal"` (default; half-Normal, unchanged) or `"lognormal"`
  (Log-Normal), letting users match the alpha prior family actually used
  when fitting. New `prior_predictive_check()` simulates data directly from
  a declared `priors` specification (including `alpha$family`) under the
  GMLTM-D's generative structure, before any data is fitted, following the
  prior predictive checking workflow of Gelman et al. (2020, Sect. 2.4);
  `plot_prior_predictive_check()` visualizes the resulting proportion-correct
  distribution against a configurable substantively plausible range.
* Fixed a bug in `GMLTM()`, `GMLTM_corr()`, `MLTM()`, and `LLTM()`: when
  `data` was a plain `matrix` (as opposed to a `data.frame`), the response
  vector was built via `unlist(data)`, which does not flatten a matrix (its
  `dim` attribute is left untouched) and made model fitting fail with a
  dimension-mismatch error from Stan. Fixed by flattening with
  `as.vector(as.matrix(data))`, which produces the same column-major order
  for both `data.frame` and `matrix` inputs.
* All package documentation and source comments are now in English.
* Internal reorganization for maintainability, with no change in behavior:
  `reliability.R` and `conditional_reliability.R` were split into smaller,
  topic-focused files (`reliability-enhanced.R`, `reliability-diagnostics.R`,
  `conditional-reliability-compare.R`); related exported functions are now
  grouped via `@family` tags in their documentation; and a package-level
  overview page (`?GMLTM`) was added.
* Removed two orphaned, unused precompiled Stan model objects that had been
  accidentally left in `inst/`.

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
