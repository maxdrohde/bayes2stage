test_that("stage2 sampling", {

  # Data parameters
  n_subjects <- 200L
  M <- 5L
  n_sampled <- as.integer(0.2 * n_subjects)

  # Sampling parameters
  sampling_type <- "intercept"
  cutoff_high <- 0.9
  cutoff_low <- 0.1
  prop_high <- 0.40
  prop_middle <- 0.20
  prop_low <- 0.40


  df <- generate_data(x_dist = "normal",
                      N = n_subjects)

  srs_df <- bayes2stage::srs_design(df, n_sampled = n_sampled)

  ods_df <- bayes2stage::ods_design(df,
                                    sampling_type = sampling_type,
                                    cutoff_high = cutoff_high,
                                    cutoff_low = cutoff_low,
                                    n_sampled = n_sampled,
                                    prop_high = prop_high,
                                    prop_middle = prop_middle,
                                    prop_low = prop_low)

  bds_df <- bayes2stage::bds_design(df,
                                    fixed_effects_formula = y ~ t + z,
                                    sampling_type = sampling_type,
                                    cutoff_high = cutoff_high,
                                    cutoff_low = cutoff_low,
                                    n_sampled = n_sampled,
                                    prop_high = prop_high,
                                    prop_middle = prop_middle,
                                    prop_low = prop_low)
})
