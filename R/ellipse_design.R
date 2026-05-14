get_mahalanobis_targets <- function(data, estimation_method, fixed_effects_formula = NULL) {
    if (estimation_method == "ods") {
        check_cols(data, c("y", "t", "id"))
        out <- data |>
            dplyr::summarize(
                coefs = list(stats::lm(y ~ t)$coefficients),
                .by = id
            )
        out$intercept <- purrr::map_dbl(out$coefs, "(Intercept)")
        out$slope <- purrr::map_dbl(out$coefs, "t")
        out$coefs <- NULL
    } else {
        if (is.null(fixed_effects_formula)) {
            cli::cli_abort("{.arg fixed_effects_formula} is required when {.arg estimation_method} is {.val bds}.")
        }
        check_cols(data, c("t", "id"))
        if (!("t" %in% all.vars(fixed_effects_formula))) {
            cli::cli_abort("{.var t} must be in {.arg fixed_effects_formula}.")
        }
        lmer_formula <- stats::update(fixed_effects_formula, . ~ . + (1 + t | id))
        mod <- lme4::lmer(formula = lmer_formula, data = data, REML = FALSE)
        re_df <- lme4::ranef(mod) |> as.data.frame()
        out <- dplyr::inner_join(
            re_df |>
                dplyr::filter(term == "(Intercept)") |>
                dplyr::select(id = grp, intercept = condval),
            re_df |>
                dplyr::filter(term == "t") |>
                dplyr::select(id = grp, slope = condval),
            by = "id"
        )
        out$id <- as.integer(as.character(out$id))
    }

    mat <- as.matrix(out[, c("intercept", "slope")])
    out$target <- stats::mahalanobis(mat, center = colMeans(mat), cov = stats::cov(mat))

    return(out)
}

#' Set x to missing based on an ellipse design
#'
#' Stratifies subjects by Mahalanobis distance in the bivariate
#' (intercept, slope) space. Subjects far from the center (high
#' Mahalanobis distance) are oversampled.
#'
#' @param data Dataset to use
#' @param estimation_method How to estimate per-subject intercepts and slopes:
#'   `"ods"` uses per-subject OLS, `"bds"` uses BLUPs from a mixed model.
#' @param fixed_effects_formula Formula for the fixed-effects when fitting the
#'   model to estimate BLUPs. Required when `estimation_method = "bds"`.
#' @param cutoff_low Which quantile to use as the cutoff for the Low category
#' @param cutoff_high Which quantile to use as the cutoff for the High category
#' @param n_sampled How many subjects should be sampled?
#' @param prop_high What proportion to sample from the High category?
#' @param prop_middle What proportion to sample from the Middle category?
#' @param prop_low What proportion to sample from the Low category?
#' @return A dataset where the x values are selected based on an ellipse design.
#'   Includes columns `target` (Mahalanobis distance), `category`,
#'   `intercept`, and `slope`.
#' @export
ellipse_design <- function(
    data,
    estimation_method = c("ods", "bds"),
    fixed_effects_formula = NULL,
    cutoff_low,
    cutoff_high,
    n_sampled,
    prop_high,
    prop_middle,
    prop_low
) {
    check_cols(data, "x")
    estimation_method <- match.arg(estimation_method)

    targets <- get_mahalanobis_targets(data, estimation_method, fixed_effects_formula)

    sampling_result <- perform_stratified_sampling(
        id_target_df = targets,
        n_sampled = n_sampled,
        cutoff_low = cutoff_low,
        cutoff_high = cutoff_high,
        prop_low = prop_low,
        prop_middle = prop_middle,
        prop_high = prop_high
    )

    stage2_df <- set_missing(data, sampling_result$selected_ids)
    stage2_df <- dplyr::left_join(stage2_df, sampling_result$id_target_df, by = "id")
    stage2_df$sampling_type <- "ellipse"

    return(stage2_df)
}
