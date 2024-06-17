## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(cmdstanr)
library(arrow)

get_draws <- function(stan_fit){
  variables <-
    c("sigma_rand_effects[1]",
      "sigma_rand_effects[2]",
      "corr_rand_effects[2,1]",
      "corr_rand_effects[1,2]",
      "alpha_main",
      "beta_t",
      "beta_x_e",
      "beta[1]",
      "sigma_main",
      "alpha_imp",
      "gamma[1]",
      "sigma_imp")

  stan_fit$draws(variables = variables,
                 format = "df")
}

## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
model <- cmdstanr::cmdstan_model(stan_file = "stan/mixed_effects_imputation.stan")

## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
run_sim <- function(){

  iter_warmup <- 250
  iter_sampling <- 250

  df <-
    bayes2stage::generate_mixed_effects_data(N = 100,
                                             Ms = 6,
                                             alpha_main = 1,
                                             beta_x_e = 1,
                                             beta_x_z = 5,
                                             beta_t = 1,
                                             beta_t_xe_interaction = 1,
                                             error_sd = 4,
                                             rand_intercept_sd = 4,
                                             rand_slope_sd = 1,
                                             rand_eff_corr = 0.5,
                                             x_cov = 0.7,
                                             xe_var = 1,
                                             xz_var = 1)

  stage2_df_ods <- bayes2stage::ods_design(df, type = "slope")
  stage2_df_blup <- bayes2stage::bds_design(df, type = "slope")
  stage2_df_srs <- bayes2stage::srs_design(df, proportion = 0.25)

################################################################################

  fit_full <- model$sample(data = bayes2stage::create_stan_data(df),
                         chains = 8,
                         parallel_chains = 8,
                         show_messages = FALSE,
                         show_exceptions = FALSE,
                         iter_warmup = iter_warmup,
                         iter_sampling = iter_sampling)

################################################################################
  cc_srs_df <-
  stage2_df_srs |> drop_na()
  ids <- data.frame(id = unique(cc_srs_df$id),
                  id_new = 1:length(unique(cc_srs_df$id)))
  cc_srs_df <- left_join(cc_srs_df, ids, by = "id")
  cc_srs_df$id <- cc_srs_df$id_new

  fit_cc_srs <- model$sample(data = bayes2stage::create_stan_data(cc_srs_df),
                         chains = 8,
                         parallel_chains = 8,
                         show_messages = FALSE,
                         show_exceptions = FALSE,
                         iter_warmup = iter_warmup,
                         iter_sampling = iter_sampling)

################################################################################
  cc_ods_df <-
  stage2_df_ods |> drop_na()
  ids <- data.frame(id = unique(cc_ods_df$id),
                  id_new = 1:length(unique(cc_ods_df$id)))
  cc_ods_df <- left_join(cc_ods_df, ids, by = "id")
  cc_ods_df$id <- cc_ods_df$id_new

  fit_cc_ods <- model$sample(data = bayes2stage::create_stan_data(cc_ods_df),
                   chains = 8,
                   parallel_chains = 8,
                   show_messages = FALSE,
                   show_exceptions = FALSE,
                   iter_warmup = iter_warmup,
                   iter_sampling = iter_sampling)

################################################################################
  cc_blup_df <-
  stage2_df_blup |> drop_na()
  ids <- data.frame(id = unique(cc_blup_df$id),
                  id_new = 1:length(unique(cc_blup_df$id)))
  cc_blup_df <- left_join(cc_blup_df, ids, by = "id")
  cc_blup_df$id <- cc_blup_df$id_new

  fit_cc_blup <- model$sample(data = bayes2stage::create_stan_data(cc_blup_df),
                   chains = 8,
                   parallel_chains = 8,
                   show_messages = FALSE,
                  show_exceptions = FALSE,
                   iter_warmup = iter_warmup,
                   iter_sampling = iter_sampling)


################################################################################
  fit_impute_srs <- model$sample(data = bayes2stage::create_stan_data(stage2_df_srs),
                 chains = 8,
                 parallel_chains = 8,
                  show_messages = FALSE,
                 show_exceptions = FALSE,
                 iter_warmup = iter_warmup,
                 iter_sampling = iter_sampling)

################################################################################
  fit_impute_ods <- model$sample(data = bayes2stage::create_stan_data(stage2_df_ods),
                 chains = 8,
                 parallel_chains = 8,
                 show_messages = FALSE,
                 show_exceptions = FALSE,
                 iter_warmup = iter_warmup,
                 iter_sampling = iter_sampling)

################################################################################
  fit_impute_blup <- model$sample(data = bayes2stage::create_stan_data(stage2_df_blup),
                 chains = 8,
                 parallel_chains = 8,
                 show_messages = FALSE,
                 show_exceptions = FALSE,
                 iter_warmup = iter_warmup,
                 iter_sampling = iter_sampling)

  d1 <- get_draws(fit_full)
  d2 <- get_draws(fit_cc_srs)
  d3 <- get_draws(fit_cc_ods)
  d4 <- get_draws(fit_cc_blup)
  d5 <- get_draws(fit_impute_srs)
  d6 <- get_draws(fit_impute_ods)
  d7 <- get_draws(fit_impute_blup)

  d1$model <- "Full Dataset"
  d2$model <- "Complete Case: SRS"
  d3$model <- "Complete Case: ODS"
  d4$model <- "Complete Case: BLUP"
  d5$model <- "Bayesian Imputation: SRS"
  d6$model <- "Bayesian Imputation: ODS"
  d7$model <- "Bayesian Imputation: BLUP"

  d <- purrr::bind_rows(d1, d2, d3, d4, d5, d6, d7)

  d |> arrow::write_parquet(
    sink = glue::glue("data/data_{Sys.time()}.parquet")
    )
}


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
purrr::walk(1:1, ~run_sim(), .progress = TRUE)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
files <- fs::dir_ls("data")

df <-
  purrr::imap(files, ~arrow::read_parquet(.x) |> mutate(dataset = .y)) |>
  purrr::bind_rows()

arrow::write_parquet(df, sink = "combined.parquet")
