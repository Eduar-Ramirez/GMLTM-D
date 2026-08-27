## Submission

This is a major update (0.1.0 -> 2.0.0). Highlights (full detail in NEWS.md):

* New `GMLTM_corr()` and `extract_correlation()`: a variant of `GMLTM()` with
  correlated latent components, and a function to summarize the estimated
  correlation matrix.
* New `student_report()` (component- and rule-level mastery report for a
  single examinee) and `student_report_batch()` (paginated batch version).
  `student_report()` replaces `informe_estudiante()`, kept as a deprecated
  alias.
* `priors$alpha` in `GMLTM()`/`GMLTM_corr()` now accepts `family = "normal"`
  (default, unchanged) or `"lognormal"`. New `prior_predictive_check()` and
  `plot_prior_predictive_check()` simulate and visualize data from a declared
  prior before fitting, following Gelman et al. (2020, Bayesian Workflow).
* Fixed a bug where `data` supplied as a plain `matrix` (rather than a
  `data.frame`) made `GMLTM()`/`GMLTM_corr()`/`MLTM()`/`LLTM()` fail.
* All source comments and documentation are now in English; internal
  reorganization of the largest files for maintainability, with no change in
  exported behavior.

No previously exported function was removed or had its signature changed in
a backward-incompatible way; `informe_estudiante()` is deprecated (with
`.Deprecated()`) rather than removed.

## R CMD check results

0 errors | 0 warnings | 2 notes

* "unable to verify current time" is a network/environment limitation of the
  checking machine, not a package issue (also present, and noted as such, in
  the 0.1.0 submission).
* "checking HTML version of manual ... NOTE" (HTML validation problems with
  `<main>`) is produced by an outdated local `tidy` (Apple's 2006 build of
  HTML Tidy, which predates HTML5); it is not reproducible on CRAN's own
  check machines, which use a current `tidy`/`libtidy`.

`R CMD check --as-cran --run-donttest` was also run locally (all examples,
including every `\donttest` example that fits a real model): all passed.
