library(bayes2stage)
library(nimble)
library(nimbleEcology)
library(nimbleMacros)
library(tidyverse)

# source("custom_nimble_functions.R")

# Test data
dataset <- bayes2stage::generate_data(N=500,
                                      x_dist = "beta_binomial",
                                      x_disp_param = 2.4,
                                      x_size = 50L)

# Plot data
bayes2stage:::data_plot(dataset)

bds_df <-
  bayes2stage::bds_design(
    dataset,
    fixed_effects_formula = y ~ t + z,
    sampling_type = "intercept",
    cutoff_high = 0.9,
    cutoff_low = 0.1,
    sampling_N = 50,
    prop_high = 0.40,
    prop_middle = 0.20,
    prop_low = 0.40)

samples <-
  fit_model(bds_df,
            main_model_covariates = c("z"),
            imputation_model_covariates = c("z"),
            imputation_model_distribution = "beta_binomial",
            x_size = 50L,
            nchain = 4,
            niter = 10000,
            nburnin = 5000)
