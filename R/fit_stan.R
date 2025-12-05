mod1 <- cmdstanr::cmdstan_model("~/r_package_development/bayes2stage/inst/stan_models/mixed_effects_imputation.stan")
mod2 <- cmdstanr::cmdstan_model("~/r_package_development/bayes2stage/inst/stan_models/mixed_effects_imputation_centered.stan")

df <- bayes2stage::generate_data(N = 2000,
                                 x_dist = "normal") |>
  bayes2stage::srs_design(N = 500)

data_list <- format_data_mcmc(df,
                              main_vars = c("z"),
                              imputation_vars = c("z"))

fit1 <- mod1$sample(
  data = data_list,
  seed = 1234,
  chains = 10,
  parallel_chains = 10,
  iter_warmup = 2000,
  iter_sampling = 2000,
  adapt_delta = 0.99
)

fit2 <- mod2$sample(
  data = data_list,
  seed = 1234,
  chains = 10,
  parallel_chains = 10,
  iter_warmup = 2000,
  iter_sampling = 2000,
  adapt_delta = 0.99
)

  posterior::as_draws(fit1) |>
  posterior::subset_draws(variable = c("beta_x", "beta[1]", "beta_t")) |>
  posterior::summarise_draws()

posterior::summarise_draws(fit1, "beta")

library(shinystan)
my_sso <- launch_shinystan(fit1)
