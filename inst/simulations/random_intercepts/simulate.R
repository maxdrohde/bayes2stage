## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(cmdstanr)
library(arrow)

get_draws <- function(stan_fit){
  variables <-
    c("beta_x_e")

  stan_fit$draws(variables = variables,
                 format = "df")
}

## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
model <- cmdstanr::cmdstan_model(stan_file = "mixed_effects_imputation.stan")

## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
run_sim <- function(){

  N <- 200
  x_cov <- 0.8

  iter_warmup <- 500
  iter_sampling <- 500

  df <-
    bayes2stage::generate_mixed_effects_data(N = N,
                                             Ms = 6,
                                             alpha_main = 1,
                                             beta_x_e = 1,
                                             beta_x_z = 1,
                                             beta_t = 1.5,
                                             beta_t_xe_interaction = 0,
                                             error_sd = 4,
                                             rand_intercept_sd = 4,
                                             rand_slope_sd = 1,
                                             rand_eff_corr = 0.3,
                                             x_cov = x_cov,
                                             xe_var = 1,
                                             xz_var = 1)

  stage2_df_bds <- bayes2stage::bds_design(df, type = "intercept")

################################################################################
  fit_impute_bds <-
    model$sample(data = bayes2stage::create_stan_data(stage2_df_bds),
                 chains = 8,
                 parallel_chains = 8,
                 show_messages = FALSE,
                 show_exceptions = FALSE,
                 iter_warmup = iter_warmup,
                 iter_sampling = iter_sampling)

  d1 <- get_draws(fit_impute_bds)
  d1$model <- "BDS: Imputation"

################################################################################
  cc_bds_df <-
    stage2_df_bds |> drop_na()
  ids <- data.frame(id = unique(cc_bds_df$id),
                    id_new = 1:length(unique(cc_bds_df$id)))
  cc_bds_df <- left_join(cc_bds_df, ids, by = "id")
  cc_bds_df$id <- cc_bds_df$id_new

  fit_cc_bds <-
    model$sample(data = bayes2stage::create_stan_data(cc_bds_df),
                 chains = 8,
                 parallel_chains = 8,
                 show_messages = FALSE,
                 show_exceptions = FALSE,
                 iter_warmup = iter_warmup,
                 iter_sampling = iter_sampling)

  d2 <- get_draws(fit_cc_bds)
  d2$model <- "BDS: Complete Case"

################################################################################

  d <- dplyr::bind_rows(d1, d2)

  d$x_cov <- x_cov
  d$N <- N

  d |> arrow::write_parquet(
    sink = glue::glue("data/data_{Sys.time()}.parquet")
    )
}


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
purrr::walk(1:5, ~run_sim(), .progress = TRUE)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
files <- fs::dir_ls("data")

df <-
  purrr::imap(files, ~arrow::read_parquet(.x) |> mutate(dataset = .y)) |>
  dplyr::bind_rows()

arrow::write_parquet(df, sink = "combined.parquet")
