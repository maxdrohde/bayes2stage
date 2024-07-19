run_sim <- function(N,
                    sampling_type,
                    sampling_design,
                    impute,
                    gamma2,
                    beta_x_e,
                    beta_x_z,
                    gamma_sd,
                    spline){

  iter_warmup <- 1000
  iter_sampling <- 2000

  stopifnot(sampling_type %in% c("intercept", "slope"))

  df <-
    bayes2stage:::generate_mixed_effects_data_quadratic(
      N = N,
      Ms = 6,
      alpha_main = 1,
      beta_x_e = beta_x_e,
      beta_x_z = beta_x_z,
      beta_t = 1.5,
      beta_t_xe_interaction = 1,
      error_sd = 4,
      rand_intercept_sd = 4,
      rand_slope_sd = 1,
      rand_eff_corr = 0.3,
      gamma0 = 1,
      gamma1 = 1,
      gamma2 = gamma2,
      gamma_sd = gamma_sd
    )

  model <- cmdstanr::cmdstan_model(stan_file = "mixed_effects_imputation.stan")

  if (sampling_design == "ODS") {
    stage2_df <- bayes2stage::ods_design(df, type = sampling_type)
  } else if (sampling_design == "BDS") {
    stage2_df <- bayes2stage::bds_design(df, type = sampling_type)
  } else if (sampling_design == "SRS") {
    stage2_df <- bayes2stage::srs_design(df, proportion = 0.25)
  } else if (sampling_design == "FULL"){
    stage2_df <- df
  } else {
    stop("Invalid argument for `sampling_design`")
  }

  if (!impute) {
    stage2_df <- tidyr::drop_na(stage2_df)

    ids <- data.frame(id = unique(stage2_df$id),
                      id_new = 1:length(unique(stage2_df$id)))

    stage2_df <- dplyr::left_join(stage2_df, ids, by = "id")
    stage2_df$id <- stage2_df$id_new
  }

  fit <- model$sample(data = bayes2stage::create_stan_data(stage2_df, spline = spline),
                      chains = 2,
                      parallel_chains = 1,
                      max_treedepth = 14L,
                      show_messages = FALSE,
                      show_exceptions = FALSE,
                      iter_warmup = iter_warmup,
                      iter_sampling = iter_sampling)

  res <-
    fit$draws(variables = c("alpha_main",
                            "beta_t",
                            "beta_x_e",
                            "beta_xe_t_interaction",
                            "beta[1]"),
                 format = "matrix") |>
    colMeans() |>
    as.data.frame.list()

  res$N <- N
  res$sampling_type <- sampling_type
  res$sampling_design <- sampling_design
  res$impute <-impute
  res$gamma2 <- gamma2
  res$true_beta_x_e <- beta_x_e
  res$true_beta_x_z <- beta_x_z
  res$gamma_sd <- gamma_sd
  res$spline <- spline
  res$iter_warmup <- iter_warmup
  res$iter_sampling <- iter_sampling

  return(res)
}

simulate <- function(i,
                     ITER){

  GRID_N <- c(200L)
  GRID_sampling_type <- c("intercept", "slope")
  GRID_sampling_design <- c("ODS", "BDS", "SRS", "FULL")
  GRID_impute <- c(TRUE, FALSE)
  GRID_gamma2 <- c(0L)
  GRID_beta_x_e <- c(1L)
  GRID_beta_x_z <- c(1L, 4L)
  GRID_gamma_sd <- c(1, 3L)
  GRID_spline <- c(FALSE)

  param_df <-
    tidyr::crossing(
      GRID_N,
      GRID_sampling_type,
      GRID_sampling_design,
      GRID_impute,
      GRID_gamma2,
      GRID_beta_x_e,
      GRID_beta_x_z,
      GRID_gamma_sd,
      GRID_spline
    )

  # Remove full and !impute
  param_df <-
    param_df |>
    dplyr::filter(!((GRID_sampling_design == "FULL") & (GRID_impute == FALSE)))

  row <- param_df[i,]

  out <- purrr::map(1:ITER, ~run_sim(N = row$GRID_N[[1]],
                                   sampling_type = row$GRID_sampling_type[[1]],
                                   sampling_design = row$GRID_sampling_design[[1]],
                                   impute = row$GRID_impute[[1]],
                                   gamma2 = row$GRID_gamma2[[1]],
                                   beta_x_e = row$GRID_beta_x_e[[1]],
                                   beta_x_z = row$GRID_beta_x_z[[1]],
                                   gamma_sd = row$GRID_gamma_sd[[1]],
                                   spline = row$GRID_spline[[1]]),
                    .progress = TRUE) |>
    purrr::list_rbind()

  readr::write_csv(out, glue::glue("results_{i}.csv"))
}

