# GMLTM

**Generalized Multicomponent Latent Trait Model for Diagnosis**

[![CRAN status](https://www.r-pkg.org/badges/version/GMLTM)](https://CRAN.R-project.org/package=GMLTM)
[![License: GPL-3](https://img.shields.io/badge/License-GPL%203-blue.svg)](https://opensource.org/licenses/GPL-3.0)

The `GMLTM` package provides Bayesian estimation of Item Response Theory models that decompose item difficulty into cognitive operations or rules. It implements the Linear Logistic Test Model (LLTM; Fischer, 1973), the Multicomponent Latent Trait Model for Diagnosis (MLTM-D; Embretson & Yang, 2013), and the Generalized Multicomponent Latent Trait Model for Diagnosis (GMLTM-D; Ramírez et al., 2024). All models are estimated via Hamiltonian Monte Carlo using Stan through the `rstan` interface.

## Installation

### From CRAN (recommended)

```r
install.packages("GMLTM")
```

### From GitHub (development version)

```r
# Install devtools if needed
install.packages("devtools")

# Install from GitHub
devtools::install_github("Eduar-Ramirez/GMLTM-D", force = TRUE)
```

## Requirements

The package requires `rstan` as the Stan backend. Install it from CRAN before using `GMLTM`:

```r
install.packages("rstan")

# Recommended configuration
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
```

## Basic Usage

```r
library(GMLTM)

# Load example data
data(analogy)

# Define Q-matrix (items x cognitive rules)
Q <- matrix(...)  # your Q-matrix here

# Define component structure
components <- list(
  transformation = c(1, 2, 3),
  relational     = c(4, 5)
)

# Fit GMLTM-D model
fit <- GMLTM(analogy, Q, components,
             iters = 2000, iter_warmup = 1000,
             chains = 2, cores = 2)

# Extract EAP estimates
fit$EAP$eta      # rule difficulty
fit$EAP$beta     # item difficulty
fit$EAP$alpha    # discrimination
fit$EAP$guessing # guessing parameters

# Marginal reliability
reliability(fit)

# Model fit
compute_model_validation(fit)
```

## Custom Prior Distributions

A key feature of `GMLTM` is support for user-defined prior distributions via the `priors` argument. This enables prior sensitivity analysis — refitting models with different priors to verify that conclusions are robust.

```r
# Conservative priors (default) — Beta(3,20) for guessing (mean ~0.13)
fit_conservative <- GMLTM(analogy, Q, components,
  iters = 2000, iter_warmup = 1000, chains = 2,
  priors = list(
    theta = list(mu = 0, sigma = 1),
    eta   = list(mu = 0, sigma = 1),
    alpha = list(mu = 0, sigma = 1),
    c     = list(shape1 = 3, shape2 = 20)
  ))

# Moderate priors — Beta(2,5) for guessing (mean ~0.29)
fit_moderate <- GMLTM(analogy, Q, components,
  iters = 2000, iter_warmup = 1000, chains = 2,
  priors = list(
    theta = list(mu = 0, sigma = 2),
    eta   = list(mu = 0, sigma = 2),
    c     = list(shape1 = 2, shape2 = 5)
  ))

# Diffuse priors — Beta(1,1) uniform for guessing
fit_diffuse <- GMLTM(analogy, Q, components,
  iters = 2000, iter_warmup = 1000, chains = 2,
  priors = list(
    theta = list(mu = 0, sigma = 5),
    eta   = list(mu = 0, sigma = 5),
    c     = list(shape1 = 1, shape2 = 1)
  ))

# Compare models using LOO-CV
loo::loo_compare(
  loo::loo(as.matrix(fit_conservative$fit, pars = "log_lik")),
  loo::loo(as.matrix(fit_moderate$fit,    pars = "log_lik")),
  loo::loo(as.matrix(fit_diffuse$fit,     pars = "log_lik"))
)
```

### Prior parameters by model

| Parameter | Distribution | Default | Models |
|-----------|-------------|---------|--------|
| `theta` (ability) | Normal(`mu`, `sigma`) | N(0, 1) | LLTM, MLTM, GMLTM |
| `eta` (rule difficulty) | Normal(`mu`, `sigma`) | N(0, 1) | LLTM, MLTM, GMLTM |
| `alpha` (discrimination) | Half-Normal(`sigma`) | HN(1) | MLTM, GMLTM |
| `c` (guessing) | Beta(`shape1`, `shape2`) | Beta(3, 20) | GMLTM only |

## Main Functions

| Function | Description |
|----------|-------------|
| `GMLTM()` | Fit the GMLTM-D model |
| `MLTM()` | Fit the MLTM-D model |
| `LLTM()` | Fit the LLTM model |
| `reliability()` | Marginal reliability estimation |
| `ppchecks()` | Posterior predictive checks (histogram) |
| `marginal_Pchecks()` | Marginal proportion checks with credible intervals |
| `compute_model_validation()` | LOO-CV and WAIC model fit indices |
| `plot_ICC_grouped()` | Item characteristic curves (grouped, 3×3 layout) |
| `plot_ICC_individual()` | Item characteristic curves (individual) |
| `conditional_reliability_tif()` | Conditional reliability via Test Information Function |
| `generate_Q_with_interactions()` | Extend Q-matrix with rule interactions |

## References

Fischer, G. H. (1973). The linear logistic test model as an instrument in educational research. *Acta Psychologica*, 37(6), 359–374.

Embretson, S. E., & Yang, X. (2013). A multicomponent latent trait model for diagnosis. *Psychometrika*, 78, 14–36.

Ramírez, E. S., Jiménez, M., Franco, V. R., & Alvarado, J. M. (2024). Delving into the complexity of analogical reasoning: A detailed exploration with the Generalized Multicomponent Latent Trait Model for Diagnosis. *Journal of Intelligence*, 12, 67. https://doi.org/10.3390/jintelligence12070067
