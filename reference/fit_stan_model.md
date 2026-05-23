# Fit a Bayesian two-stage model using Stan

Fits a mixed effects model with imputation using Stan via the
instantiate package.

## Usage

``` r
fit_stan_model(
  data,
  main_model_formula,
  imputation_model_formula,
  imputation_distribution = c("normal", "bernoulli", "beta_binomial",
    "negative_binomial"),
  mixture_components = NULL,
  inference_method = c("mcmc", "pathfinder", "laplace", "optimize"),
  pathfinder_draws = 1000L,
  pathfinder_num_paths = 4L,
  laplace_draws = 1000L,
  optimize_algorithm = c("lbfgs", "bfgs", "newton"),
  optimize_iter = 2000L,
  use_pathfinder_init = FALSE,
  n_chains = 4L,
  iter_warmup = 1000L,
  iter_sampling = 1000L,
  adapt_delta = 0.8,
  seed = 777L,
  parallel_chains = 1L
)
```

## Arguments

- data:

  A data frame containing the outcome and covariates

- main_model_formula:

  One-sided formula or string for covariates in the main model (e.g.,
  `~ age + splines::ns(bmi, 3)`). Intercept is automatically removed.

- imputation_model_formula:

  One-sided formula or string for covariates in the imputation model
  (e.g., `~ age + factor(site)`). Intercept is automatically removed.

- imputation_distribution:

  Distribution for the imputation model: "normal" for continuous x,
  "bernoulli" for binary x, "beta_binomial" for bounded count data, or
  "negative_binomial" for unbounded count data (default: "normal")

- mixture_components:

  Number of mixture components for the imputation model, or NULL
  (default) for a single-component model. When specified (must be \>=
  2), uses a mixture of normals for the imputation distribution, which
  can handle multimodal covariate distributions. Only available when
  `imputation_distribution = "normal"`.

- inference_method:

  Inference method to use:

  - "mcmc" (default): Full MCMC sampling via NUTS. Most accurate but
    slowest.

  - "pathfinder": Variational inference via Pathfinder algorithm. Fast
    approximate inference with posterior draws. Good for initial
    exploration.

  - "laplace": Laplace approximation. Finds posterior mode and
    approximates with Gaussian. Very fast with good accuracy for
    well-behaved posteriors.

  - "optimize": Maximum a posteriori (MAP) point estimate only. Fastest
    but no uncertainty quantification. Returns CmdStanMLE object (no
    posterior draws).

- pathfinder_draws:

  Number of approximate posterior draws from Pathfinder (default: 1000L)

- pathfinder_num_paths:

  Number of Pathfinder paths to run (default: 4L)

- laplace_draws:

  Number of approximate posterior draws from Laplace (default: 1000L)

- optimize_algorithm:

  Optimization algorithm for MAP estimation: "lbfgs" (default), "bfgs",
  or "newton"

- optimize_iter:

  Maximum iterations for optimization (default: 2000L)

- use_pathfinder_init:

  Logical; if TRUE, use Pathfinder variational inference to initialize
  MCMC chains. Only applies when `inference_method = "mcmc"`. This can
  dramatically improve sampling efficiency for complex models with many
  latent parameters (default: FALSE)

- n_chains:

  Number of MCMC chains (default: 4)

- iter_warmup:

  Number of warmup iterations per chain (default: 1000)

- iter_sampling:

  Number of sampling iterations per chain (default: 1000)

- adapt_delta:

  Target acceptance rate for HMC (default: 0.8)

- seed:

  Random seed for reproducibility (default: 777L)

- parallel_chains:

  Number of chains to run in parallel (default: 1L)

## Value

A CmdStan fit object. Type depends on `inference_method`: CmdStanMCMC
(mcmc), CmdStanPathfinder (pathfinder), CmdStanLaplace (laplace), or
CmdStanMLE (optimize). All support `$summary()`, but only
mcmc/pathfinder/laplace support `$draws()` for posterior samples.
