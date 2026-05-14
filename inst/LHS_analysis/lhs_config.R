################################################################################
# LHS Analysis Configuration
#
# Shared configuration for LHS case study scripts. Source this file to access:
#   - All sampling/inference constants
#   - MODEL_SPECS list
#   - load_lhs_data(): Returns cleaned LHS data frame
################################################################################

################################################################################
#                    SAMPLING & INFERENCE SETTINGS
################################################################################

N_REPS <- 100L
N_CORES <- 14L

SAMPLING_FRACTION <- 0.25
CUTOFF_HIGH <- 0.9
CUTOFF_LOW <- 0.1
PROP_HIGH <- 0.40
PROP_MIDDLE <- 0.20
PROP_LOW <- 0.40

INFERENCE_ARGS <- list(
    n_chains = 2L,
    iter_warmup = 500L,
    iter_sampling = 500L,
    parallel_chains = 1L,
    adapt_delta = 0.8,
    seed = 777L
)

################################################################################
#                    MODEL SPECIFICATIONS
################################################################################

MODEL_SPECS <- list(
    no_covariates = list(
        label = "No covariates (normal)",
        main_formula = "~ 1",
        imp_formula = "~ 1",
        imp_distribution = "normal",
        bds_formula = y ~ t
    ),
    schildcrout = list(
        label = "Schildcrout 2019 (normal)",
        main_formula = paste(
            "~ packyear + f31cigs + bmi0 + female + age + site +",
            "bmi_change + t:age + t:female + t:site"
        ),
        imp_formula = "~ packyear + f31cigs + bmi0 + female + age + site",
        imp_distribution = "normal",
        bds_formula = y ~ t + packyear + f31cigs + bmi0 + female + age + site +
            bmi_change + t:age + t:female + t:site
    )
)

################################################################################
#                    DATA LOADING
################################################################################

N_SNPS <- 43L

#' Load and prepare LHS data
#'
#' @param subset Logical; if TRUE, use a random subset of subjects
#' @param n_subset Integer; number of subjects in subset
#' @return A data frame ready for analysis
load_lhs_data <- function(subset = FALSE, n_subset = 200L) {
    # Resolve path relative to this script's location
    script_dir <- if (exists("SCRIPT_DIR")) {
        SCRIPT_DIR
    } else {
        "."
    }
    data_path <- file.path(script_dir, "clean_data", "lhs_cleaned.parquet")

    df_raw <- arrow::read_parquet(data_path)

    df <- df_raw |>
        dplyr::rename(x = x_count) |>
        dplyr::select(-x_mean, -rs177852) |>
        dplyr::mutate(site = as.factor(site), n_trials = N_SNPS)

    df$y <- df$fvc

    if (subset) {
        set.seed(777)
        subset_ids <- sample(unique(df$id), size = n_subset)
        df <- dplyr::filter(df, id %in% subset_ids)
        df$id <- as.integer(factor(df$id))
        cli::cli_alert_warning("Using subset: {n_subset} subjects")
    }

    df
}

################################################################################
#                    DISPLAY LABELS
################################################################################

DESIGN_LABELS <- c(
    full = "Full data",
    srs = "SRS",
    srs_no_imp = "SRS (no imp)",
    ods_intercept = "ODS (intercept)",
    ods_slope = "ODS (slope)",
    bds_intercept = "BDS (intercept)",
    bds_slope = "BDS (slope)"
)

PARAM_LABELS <- c(
    beta_x = "X effect",
    beta_x_t_interaction = "X \u00d7 Time",
    beta_t = "Time slope",
    alpha_main = "Intercept"
)

SPEC_LABELS <- c(
    no_covariates = "No covariates | Imputation distribution: Normal",
    schildcrout = "Schildcrout 2019 covariates | Imputation distribution: Normal"
)

KEY_PARAMS <- c("beta_x", "beta_x_t_interaction", "beta_t", "alpha_main")
