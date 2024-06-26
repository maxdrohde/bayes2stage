% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_mixed_effects_data.R
\name{generate_mixed_effects_data}
\alias{generate_mixed_effects_data}
\title{Generate continuous mixed-effects data for the purpose of testing the Bayesian
Two-Stage design methods. We generate longitudinal data with random intercepts
and slopes. There are two continuous covariate, the expensive covariate
x_e and the inexpensive covariate x_z.}
\usage{
generate_mixed_effects_data(
  N = 1000,
  Ms = 5,
  alpha_main = 1,
  beta_x_e = 2,
  beta_x_z = 0.5,
  beta_t = 2,
  beta_t_xe_interaction = 1,
  error_sd = 1,
  rand_intercept_sd = 3,
  rand_slope_sd = 3,
  rand_eff_corr = 0.4,
  x_cov = 0.6,
  xe_var = 1,
  xz_var = 1
)
}
\arguments{
\item{N}{Number of subjects to generate}

\item{Ms}{A vector containing the number of timepoints to generate}

\item{alpha_main}{Intercept}

\item{beta_x_e}{x_e effect}

\item{beta_x_z}{x_z effect}

\item{beta_t}{time effect}

\item{beta_t_xe_interaction}{time by x_e interaction}

\item{error_sd}{SD of error term}

\item{rand_intercept_sd}{SD of random intercepts}

\item{rand_slope_sd}{SD of random slopes}

\item{rand_eff_corr}{Correlation between random effects}

\item{x_cov}{Covariance of the covariates}

\item{xe_var}{Variance of x_e}

\item{xz_var}{Variance of x_z}
}
\value{
A simulated mixed-effects dataset
}
\description{
Generate continuous mixed-effects data for the purpose of testing the Bayesian
Two-Stage design methods. We generate longitudinal data with random intercepts
and slopes. There are two continuous covariate, the expensive covariate
x_e and the inexpensive covariate x_z.
}
