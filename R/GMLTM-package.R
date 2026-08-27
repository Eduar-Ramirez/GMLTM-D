#' GMLTM: Generalized Multicomponent Latent Trait Model for Diagnosis
#'
#' @description
#' Bayesian estimation of Item Response Theory models that decompose item
#' difficulty into cognitive operations or rules: the Linear Logistic Test
#' Model (LLTM), the Multicomponent Latent Trait Model for Diagnosis
#' (MLTM-D), and the Generalized Multicomponent Latent Trait Model for
#' Diagnosis (GMLTM-D). All models are fit via Hamiltonian Monte Carlo using
#' \pkg{rstan}.
#'
#' @section Model fitting:
#' \code{\link{LLTM}}, \code{\link{MLTM}}, \code{\link{GMLTM}}, and
#' \code{\link{GMLTM_corr}} fit the four supported models and share the
#' same \code{data}/\code{Q}/\code{priors} interface. See
#' \code{vignette("GMLTM-intro", package = "GMLTM")} for a worked example.
#'
#' @section Reliability:
#' \code{\link{reliability}} and \code{\link{enhanced_mltm_reliability}}
#' estimate marginal and component-level reliability from a fitted model.
#' \code{\link{conditional_reliability_tif}} and
#' \code{\link{compare_conditional_reliability}} estimate reliability as a
#' function of the ability scale. \code{\link{quick_reliability_check}},
#' \code{\link{check_reliability_data_quality}}, and
#' \code{\link{demo_reliability_analysis}} support quick diagnostics.
#'
#' @section Model checking and diagnostics:
#' \code{\link{ppchecks}} and \code{\link{marginal_Pchecks}} perform
#' posterior predictive checks; \code{\link{compute_model_validation}}
#' computes WAIC/LOO-based validation statistics;
#' \code{\link{extract_correlation}} inspects correlations among latent
#' components; \code{\link{student_report}} builds a per-examinee mastery
#' report.
#'
#' @section Plotting and item structure:
#' \code{\link{plot_ICC_grouped}} and \code{\link{plot_ICC_individual}}
#' draw item characteristic curves; \code{\link{generate_Q_with_interactions}}
#' expands a Q-matrix with rule-interaction columns.
#'
#' @references
#' Fischer, G. H. (1973). The linear logistic test model as an instrument
#' in educational research. \emph{Acta Psychologica}, \bold{37}(6), 359--374.
#' \doi{10.1016/0001-6918(73)90003-6}
#'
#' Embretson, S. E., & Yang, X. (2013). A multicomponent latent trait model
#' for diagnosis. \emph{Psychometrika}, \bold{78}, 14--36.
#' \doi{10.1007/s11336-012-9296-y}
#'
#' Embretson, S. E. (2019). Diagnostic modeling of skill hierarchies and
#' cognitive processes with MLTM-D. In M. von Davier & Y.-S. Lee (Eds.),
#' \emph{Handbook of Diagnostic Classification Models} (pp. 185--208).
#' Springer. \doi{10.1007/978-3-030-05584-4_9}
#'
#' Ramirez, E. S., Jimenez, M., Franco, V. R., & Alvarado, J. M. (2024).
#' Delving into the complexity of analogical reasoning: A detailed exploration
#' with the Generalized Multicomponent Latent Trait Model for Diagnosis.
#' \emph{Journal of Intelligence}, \bold{12}, 67.
#' \doi{10.3390/jintelligence12070067}
#'
#' @keywords internal
"_PACKAGE"
