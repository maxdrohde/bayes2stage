################################################################################
# Run one LHS replicate
#
# Usage: Rscript run_rep.R <rep_id>
#
#   rep_id = 0:  Fit full-data models (all subjects measured)
#   rep_id >= 1: Draw subsamples and fit all (design x spec x method) combos
#
# Output: results/rep_<rep_id>.parquet
# Skips if output file already exists.
################################################################################

################################################################################
# Activate renv from package root
################################################################################

Sys.setenv(RENV_PROJECT = normalizePath("../.."))
source("../../renv/activate.R")

################################################################################
# Command Line Arguments
################################################################################

args <- commandArgs(trailingOnly = TRUE)
rep_id <- as.integer(args[[1]])

SCRIPT_DIR <- normalizePath(".")
source("lhs_config.R")

output_file <- glue::glue("results/rep_{rep_id}.parquet")
if (fs::file_exists(output_file)) {
    cli::cli_alert_info("Results already exist for rep {rep_id}: {output_file}")
    quit(save = "no")
}

fs::dir_create("results")

################################################################################
# Load Data
################################################################################

df <- load_lhs_data()
N_SUBJECTS <- length(unique(df$id))
N_SAMPLED <- as.integer(SAMPLING_FRACTION * N_SUBJECTS)

cli::cli_h1("LHS Replicate: rep_id={rep_id}")
cli::cli_alert_info("N subjects: {N_SUBJECTS}, N sampled: {N_SAMPLED}")

################################################################################
# Fitting helper
################################################################################

fit_and_summarize <- function(data, spec_name, design, method = "stan") {
    spec <- MODEL_SPECS[[spec_name]]
    tryCatch({
        fit <- if (method == "stan") {
            bayes2stage::fit_stan_model(
                data = data,
                main_model_formula = spec$main_formula,
                imputation_model_formula = spec$imp_formula,
                imputation_distribution = spec$imp_distribution,
                inference_method = "mcmc",
                n_chains = INFERENCE_ARGS$n_chains,
                iter_warmup = INFERENCE_ARGS$iter_warmup,
                iter_sampling = INFERENCE_ARGS$iter_sampling,
                parallel_chains = INFERENCE_ARGS$parallel_chains,
                adapt_delta = INFERENCE_ARGS$adapt_delta,
                seed = INFERENCE_ARGS$seed
            )
        } else {
            bayes2stage::fit_acml_ods(
                ods_df = data,
                cutoff_low = CUTOFF_LOW,
                cutoff_high = CUTOFF_HIGH,
                main_model_formula = spec$main_formula
            )
        }
        summ <- bayes2stage::model_summary(fit)
        summ$model_spec <- spec_name
        summ$design <- design
        summ$replicate <- rep_id
        summ$method <- method
        summ
    }, error = function(e) {
        data.frame(
            parameter = NA_character_,
            model_spec = spec_name,
            design = design,
            replicate = rep_id,
            method = method,
            error = conditionMessage(e)
        )
    })
}

################################################################################
# Run
################################################################################

if (rep_id == 0L) {

    # Full-data fits: all subjects measured, no subsampling
    cli::cli_h2("Full-data fits")
    fits <- purrr::map(names(MODEL_SPECS), \(spec_name) {
        cli::cli_alert("Fitting {spec_name} (full data)")
        fit_and_summarize(df, spec_name, design = "full")
    })

} else {

    set.seed(777L + rep_id)
    fits <- list()

    # -- SRS --
    cli::cli_h2("SRS")
    srs_df <- bayes2stage::srs_design(df, N_SAMPLED)
    srs_complete <- dplyr::filter(srs_df, !is.na(x)) |>
        dplyr::mutate(id = as.integer(factor(id)))

    for (spec_name in names(MODEL_SPECS)) {
        cli::cli_alert("Fitting {spec_name} (srs)")
        fits <- c(fits, list(
            fit_and_summarize(srs_df, spec_name, "srs"),
            fit_and_summarize(srs_complete, spec_name, "srs_no_imp")
        ))
    }

    # -- ODS and BDS --
    for (sampling_type in c("intercept", "slope")) {
        cli::cli_h2("ODS/BDS ({sampling_type})")

        ods_df <- bayes2stage::ods_design(
            df, sampling_type = sampling_type,
            cutoff_high = CUTOFF_HIGH, cutoff_low = CUTOFF_LOW,
            n_sampled = N_SAMPLED,
            prop_high = PROP_HIGH, prop_middle = PROP_MIDDLE,
            prop_low = PROP_LOW
        )

        ods_label <- paste0("ods_", sampling_type)
        bds_label <- paste0("bds_", sampling_type)

        for (spec_name in names(MODEL_SPECS)) {
            bds_df <- bayes2stage::bds_design(
                df,
                fixed_effects_formula = MODEL_SPECS[[spec_name]]$bds_formula,
                sampling_type = sampling_type,
                cutoff_high = CUTOFF_HIGH, cutoff_low = CUTOFF_LOW,
                n_sampled = N_SAMPLED,
                prop_high = PROP_HIGH, prop_middle = PROP_MIDDLE,
                prop_low = PROP_LOW
            )

            cli::cli_alert("Fitting {spec_name} ({ods_label}, {bds_label})")
            new_fits <- list(
                fit_and_summarize(ods_df, spec_name, ods_label, "stan"),
                fit_and_summarize(bds_df, spec_name, bds_label, "stan")
            )

            if (MODEL_SPECS[[spec_name]]$imp_distribution != "beta_binomial") {
                new_fits <- c(new_fits, list(
                    fit_and_summarize(ods_df, spec_name, ods_label, "acml")
                ))
            }
            fits <- c(fits, new_fits)
        }
    }
}

################################################################################
# Save
################################################################################

res <- dplyr::bind_rows(fits)

if (!"error" %in% names(res)) {
    res$error <- NA_character_
}

arrow::write_parquet(res, output_file)
cli::cli_alert_success("Saved {nrow(res)} rows to {output_file}")
