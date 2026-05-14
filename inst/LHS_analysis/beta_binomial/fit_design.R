################################################################################
# Fit one design for the credible regions analysis (beta-binomial imputation)
#
# Usage: Rscript fit_design.R <design_name>
#
# Design names:
#   full, srs, ods_intercept, ods_slope, ods_ellipse,
#   bds_intercept, bds_slope, bds_ellipse
#
# Output: cr_results/<design_name>.parquet
# Skips if output file already exists.
################################################################################

Sys.setenv(RENV_PROJECT = normalizePath("../../.."))
source("../../../renv/activate.R")

args <- commandArgs(trailingOnly = TRUE)
design_name <- args[[1]]

SCRIPT_DIR <- normalizePath("..")
source("../lhs_config.R")

# Use schildcrout spec but override to beta-binomial imputation
SPEC <- MODEL_SPECS$schildcrout
SPEC$imp_distribution <- "beta_binomial"

output_file <- glue::glue("cr_results/{design_name}.parquet")
if (fs::file_exists(output_file)) {
    cli::cli_alert_info("Already exists: {output_file}")
    quit(save = "no")
}

fs::dir_create("cr_results")

################################################################################
# Load data and apply design
################################################################################

df <- load_lhs_data()
N_SUBJECTS <- length(unique(df$id))
N_SAMPLED <- as.integer(SAMPLING_FRACTION * N_SUBJECTS)

set.seed(777L)

cli::cli_h1("Fitting design: {design_name} (beta-binomial)")

DESIGN_ARGS <- list(
    cutoff_high = CUTOFF_HIGH, cutoff_low = CUTOFF_LOW,
    n_sampled = N_SAMPLED,
    prop_high = PROP_HIGH, prop_middle = PROP_MIDDLE, prop_low = PROP_LOW
)

design_df <- switch(design_name,
    full = df,
    srs = bayes2stage::srs_design(df, N_SAMPLED),
    ods_intercept = do.call(bayes2stage::ods_design,
        c(list(data = df, sampling_type = "intercept"), DESIGN_ARGS)),
    ods_slope = do.call(bayes2stage::ods_design,
        c(list(data = df, sampling_type = "slope"), DESIGN_ARGS)),
    ods_ellipse = do.call(bayes2stage::ellipse_design,
        c(list(data = df, estimation_method = "ods"), DESIGN_ARGS)),
    bds_intercept = do.call(bayes2stage::bds_design,
        c(list(data = df, fixed_effects_formula = SPEC$bds_formula,
               sampling_type = "intercept"), DESIGN_ARGS)),
    bds_slope = do.call(bayes2stage::bds_design,
        c(list(data = df, fixed_effects_formula = SPEC$bds_formula,
               sampling_type = "slope"), DESIGN_ARGS)),
    bds_ellipse = do.call(bayes2stage::ellipse_design,
        c(list(data = df, estimation_method = "bds",
               fixed_effects_formula = SPEC$bds_formula), DESIGN_ARGS)),
    cli::cli_abort("Unknown design: {design_name}")
)

################################################################################
# Fit model
################################################################################

fit <- bayes2stage::fit_stan_model(
    data = design_df,
    main_model_formula = SPEC$main_formula,
    imputation_model_formula = SPEC$imp_formula,
    imputation_distribution = SPEC$imp_distribution,
    inference_method = "mcmc",
    n_chains = 4L,
    iter_warmup = 1000L,
    iter_sampling = 2000L,
    parallel_chains = 2L,
    adapt_delta = INFERENCE_ARGS$adapt_delta,
    seed = INFERENCE_ARGS$seed
)

################################################################################
# Extract draws + diagnostics and save
################################################################################

draws <- posterior::as_draws_df(
    fit$draws(variables = c("beta_x", "beta_x_t_interaction"))
)

diag <- fit$diagnostic_summary(quiet = TRUE)

out <- data.frame(
    beta_x = draws$beta_x,
    beta_xt = draws$beta_x_t_interaction,
    design = design_name,
    divergent = sum(diag$num_divergent),
    max_treedepth = sum(diag$num_max_treedepth),
    ebfmi_min = min(diag$ebfmi[is.finite(diag$ebfmi)])
)

arrow::write_parquet(out, output_file)
cli::cli_alert_success("Saved {nrow(out)} draws to {output_file}")
