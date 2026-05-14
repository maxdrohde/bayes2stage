################################################################################
# Fit one design (no covariates spec) for the credible regions analysis
#
# Usage: Rscript fit_design.R <design_name> <distribution>
#
# Design names:
#   full, srs, ods_intercept, ods_slope, ods_ellipse,
#   bds_intercept, bds_slope, bds_ellipse
#
# Distributions: normal, beta_binomial, negative_binomial
#
# Output: cr_results/<design_name>_<distribution>.parquet
# Skips if output file already exists.
################################################################################

Sys.setenv(RENV_PROJECT = normalizePath("../../.."))
source("../../../renv/activate.R")

args <- commandArgs(trailingOnly = TRUE)
design_name <- args[[1]]
distribution <- args[[2]]

SCRIPT_DIR <- normalizePath("..")
source("../lhs_config.R")

SPEC <- MODEL_SPECS$no_covariates
SPEC$imp_distribution <- distribution

output_file <- glue::glue("cr_results/{design_name}_{distribution}.parquet")
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

cli::cli_h1("Fitting design: {design_name} ({distribution}, no covariates)")

DESIGN_ARGS <- list(
    cutoff_high = CUTOFF_HIGH, cutoff_low = CUTOFF_LOW,
    n_sampled = N_SAMPLED,
    prop_high = PROP_HIGH, prop_middle = PROP_MIDDLE, prop_low = PROP_LOW
)

design_df <- switch(design_name,
    full = df,
    srs = bayes2stage::srs_design(df, N_SAMPLED),
    srs_no_imp = {
        srs_df <- bayes2stage::srs_design(df, N_SAMPLED)
        dplyr::filter(srs_df, !is.na(x)) |>
            dplyr::mutate(id = as.integer(factor(id)))
    },
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
    n_chains = 2L,
    iter_warmup = 500L,
    iter_sampling = 500L,
    parallel_chains = 1L,
    adapt_delta = INFERENCE_ARGS$adapt_delta,
    seed = INFERENCE_ARGS$seed
)

################################################################################
# Extract draws + diagnostics and save
################################################################################

available <- fit$metadata()$stan_variables
key_params <- bayes2stage::get_key_parameters(available)
# get_key_parameters returns indexed names like sigma_re[1], but Stan metadata
# uses the base name sigma_re. Map indexed names back to base names.
base_params <- unique(gsub("\\[.*\\]$", "", key_params))
draw_vars <- intersect(base_params, available)

draws <- posterior::as_draws_df(fit$draws(variables = draw_vars))

diag <- fit$diagnostic_summary(quiet = TRUE)

# Drop posterior metadata columns (.chain, .iteration, .draw)
drop_cols <- c(".chain", ".iteration", ".draw")
out <- as.data.frame(draws)[, setdiff(names(draws), drop_cols)]
out$design <- design_name
out$imputation <- distribution
out$divergent <- sum(diag$num_divergent)
out$max_treedepth <- sum(diag$num_max_treedepth)
out$ebfmi_min <- min(diag$ebfmi[is.finite(diag$ebfmi)])

arrow::write_parquet(out, output_file)
cli::cli_alert_success("Saved {nrow(out)} draws to {output_file}")
