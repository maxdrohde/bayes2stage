test_that("model fitting", {
  df <- generate_data(x_dist = "normal", N = 100)

  bayes2stage:::fit_model(
    df,
    main_model_covariates = c("z"),
    imputation_model_covariates = c("z"),
    imputation_model_distribution = "normal",
    n_chains = 1,
    niter = 200,
    nburnin = 100,
    print_code = TRUE)
})
