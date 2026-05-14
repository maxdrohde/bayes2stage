#' Fit ACML model for outcome-dependent sampling
#'
#' Fits an ascertainment-corrected maximum likelihood model for ODS designs.
#'
#' @param ods_df Data frame from ODS sampling. Must contain columns: `id`, `target`,
#'   `selected`, `category`, `x`, `y`, `t`, `sampling_type`.
#' @param cutoff_low Lower quantile cutoff (between 0 and 1) for sampling regions.
#' @param cutoff_high Upper quantile cutoff (between 0 and 1) for sampling regions.
#' @param main_model_formula One-sided formula or string for additional covariates
#'   in the main model (e.g., `~ z`). Terms `t`, `x`, and `x:t` are always included.
#'   Default: `~ 1` (no additional covariates).
#' @param random_effects_formula Formula for random effects. Default: `~ 1 + t`
#'   (random intercept and slope on time).
#' @return A data frame with class `"acml_fit"` containing parameter estimates
#'   and confidence intervals. Columns: `variable`, `mean`, `sd`, `2.5%`, `97.5%`.
#' @export
fit_acml_ods <- function(ods_df,
                         cutoff_low,
                         cutoff_high,
                         main_model_formula = ~ 1,
                         random_effects_formula = ~ 1 + t) {

    # Validate inputs
    check_cols(ods_df, c("id", "target", "selected", "category", "x", "y", "t", "sampling_type"))

    valid_types <- c("intercept", "slope")
    actual_types <- unique(ods_df$sampling_type)
    unsupported <- setdiff(actual_types, valid_types)
    if (length(unsupported) > 0L) {
        cli::cli_abort(c(
            "{.fun fit_acml_ods} only supports univariate ODS designs.",
            "x" = "Found unsupported sampling type{?s}: {.val {unsupported}}.",
            "i" = "Ellipse designs require Bayesian estimation via {.fun fit_stan_model}."
        ))
    }

    if (!is.numeric(cutoff_low) || length(cutoff_low) != 1L ||
        cutoff_low < 0 || cutoff_low > 1) {
        cli::cli_abort("{.arg cutoff_low} must be a number between 0 and 1.")
    }
    if (!is.numeric(cutoff_high) || length(cutoff_high) != 1L ||
        cutoff_high < 0 || cutoff_high > 1) {
        cli::cli_abort("{.arg cutoff_high} must be a number between 0 and 1.")
    }
    if (cutoff_low >= cutoff_high) {
        cli::cli_abort("{.arg cutoff_low} must be less than {.arg cutoff_high}.")
    }

    # Convert string to formula if needed
    if (is.character(main_model_formula)) {
        main_model_formula <- stats::as.formula(main_model_formula)
    }
    if (is.character(random_effects_formula)) {
        random_effects_formula <- stats::as.formula(random_effects_formula)
    }

    # Build full fixed effects formula: y ~ t + x + x:t + <covariates>
    covariate_terms <- labels(stats::terms(main_model_formula))
    if (length(covariate_terms) > 0L) {
        full_formula_str <- paste0("y ~ t + x + x:t + ",
                                   paste(covariate_terms, collapse = " + "))
    } else {
        full_formula_str <- "y ~ t + x + x:t"
    }
    fixed_effects_formula <- stats::as.formula(full_formula_str)

    cli::cli_alert_info("Fitting ACML model...")

    # Create a data frame of ID and target (one row per ID)
    subj <-
        ods_df |>
        dplyr::select(id, target) |>
        dplyr::distinct(id, .keep_all = TRUE)

    # Compute the quantiles for the target intercept or slope
    q_low <- as.numeric(stats::quantile(subj$target,
                                        probs = cutoff_low,
                                        na.rm = TRUE))

    q_high <- as.numeric(stats::quantile(subj$target,
                                         probs = cutoff_high,
                                         na.rm = TRUE))

    # Filter to selected subjects
    samp <- dplyr::filter(ods_df, selected)

    # acml.lmem2 uses data[, col] subsetting which requires a plain data.frame
    # (tibbles return 1-column tibbles instead of vectors)
    samp <- as.data.frame(samp)

    # Don't use weights
    samp$Weights <- 1

    # Expand to matrices (same values repeated for each sampled row)
    cutM <- matrix(rep(c(q_low, q_high), nrow(samp)),
                   ncol = 2,
                   byrow = TRUE)

    # Compute ODS sampling probabilities
    lv <- c("Low", "Middle", "High")
    by_id <- dplyr::group_split(ods_df, id)
    cat_by_id <- purrr::map_chr(by_id, \(x) as.character(x$category[[1]]))
    sel_by_id <- purrr::map_lgl(by_id, \(x) any(!is.na(x$x)))
    N_by <- table(factor(cat_by_id, levels = lv))
    nsel_by <- table(factor(cat_by_id[sel_by_id], levels = lv))
    SampProb_vec <- ifelse(N_by > 0, as.numeric(nsel_by) / as.numeric(N_by), 0)
    names(SampProb_vec) <- lv
    probM <- matrix(SampProb_vec, nrow(samp), length(lv), byrow = TRUE)

    # Build lmer formula with random effects
    re_terms <- paste0("(", deparse(random_effects_formula[[2]]), " | id)")
    lmer_formula <- stats::as.formula(paste(deparse1(fixed_effects_formula), "+", re_terms))

    # Get sensible starting values from a naive lmer model
    m0 <- tryCatch(
        lme4::lmer(lmer_formula, data = samp, REML = FALSE),
        error = function(e) {
            cli::cli_abort(c(
                "Initial lmer model failed:",
                "x" = e$message
            ), call = NULL)
        }
    )

    beta0 <- lme4::fixef(m0)
    vc <- lme4::VarCorr(m0)$id
    sd_b0 <- attr(vc, "stddev")[1]
    sd_b1 <- attr(vc, "stddev")[2]
    rho01 <- attr(vc, "correlation")[1, 2]
    sd_e <- stats::sigma(m0)

    z_rho <- log((1 + rho01) / (1 - rho01))
    InitVals <- c(beta0, log(sd_b0), log(sd_b1), z_rho, log(sd_e))

    # Fit ACML
    fit <- tryCatch(
        acml.lmem2(
            formula.fixed = fixed_effects_formula,
            formula.random = random_effects_formula,
            data = samp,
            id = id,
            w.function = as.character(samp$sampling_type),
            InitVals = InitVals,
            cutpoints = cutM,
            SampProb = probM,
            Weights = Weights
        ),
        error = function(e) {
            cli::cli_abort(c(
                "ACML fitting failed:",
                "x" = e$message
            ), call = NULL)
        }
    )

    ## 95% Wald CIs (robust)
    est <- fit$coefficients
    V <- fit$robcov
    se <- sqrt(diag(V))
    z_crit <- stats::qnorm(0.975)

    # Build variable names from lme4 model matrix (accounts for factor expansion)
    re_terms_parsed <- all.vars(random_effects_formula)
    n_re <- length(re_terms_parsed) + 1L  # +1 for intercept
    n_re_corr <- choose(n_re, 2)

    fixed_names <- c("alpha_main", paste0("beta_", names(beta0)[-1L]))
    re_sd_names <- paste0("log_sigma_re[", seq_len(n_re), "]")
    re_corr_names <- if (n_re_corr > 0L) "z_corr_re" else character(0)
    error_names <- "log_sigma_main"

    var_names <- c(fixed_names, re_sd_names, re_corr_names, error_names)

    ci_est <- data.frame(
        variable = var_names,
        mean = est,
        sd = se,
        `2.5%` = est - z_crit * se,
        `97.5%` = est + z_crit * se,
        check.names = FALSE
    )

    class(ci_est) <- c("acml_fit", class(ci_est))

    cli::cli_alert_success("ACML fitting complete.")

    return(ci_est)
}
