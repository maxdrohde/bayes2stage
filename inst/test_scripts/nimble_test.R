library(bayes2stage)
library(nimble)
library(tidyverse)

# Test data
dataset <- bayes2stage::generate_data(N=200)

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

ods_df <-
  bayes2stage::ods_design(
    dataset,
    sampling_type = "intercept",
    cutoff_high = 0.9,
    cutoff_low = 0.1,
    sampling_N = 50,
    prop_high = 0.40,
    prop_middle = 0.20,
    prop_low = 0.40)

#bds_df <- bds_df |> tidyr::drop_na()
#bds_df$id <- as.integer(factor(bds_df$id))

samples <- fit_nimble_model(bds_df, print_model_details = TRUE)

check_mcmc(samples)

samples <- fit_nimble_model(bds_df, niter = 20000, nchains = 4)

check_mcmc(samples)
