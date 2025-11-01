#' fit_acml_ods
#'
#' @param ods_df ods_df
#' @return params with conf ints
#' @export
fit_acml_ods <- function(ods_df,
                         cutoff_low,
                         cutoff_high){

  # Create a data frame of ID and target
  # One row per ID
  subj <-
    ods_df |>
    dplyr::select(id, target) |>
    dplyr::distinct(id, .keep_all = TRUE)

  # Compute the quantiles for the target intercept or slope
  q_low  <- as.numeric(quantile(subj$target,
                                probs = cutoff_low,
                                na.rm = TRUE))

  q_high <- as.numeric(quantile(subj$target,
                                probs = cutoff_high,
                                na.rm = TRUE))

  # Filter to selected subjects
  samp <- filter(ods_df, selected)

  # Don't use weights
  samp$Weights <- 1

  # Expand to matrices (same values repeated for each sampled row)
  cutM  <- matrix(rep(c(q_low, q_high), nrow(samp)),
                  ncol = 2,
                  byrow = TRUE)

  # Compute ODS sampling probabilities
  lv <- c("Low","Middle","High")
  by_id <- dplyr::group_split(ods_df, id)
  cat_by_id <- map_chr(by_id, ~ as.character(.x$category[[1]]))
  sel_by_id <- map_lgl(by_id, ~ any(!is.na(.x$x)))
  N_by    <- table(factor(cat_by_id, levels = lv))
  nsel_by <- table(factor(cat_by_id[sel_by_id], levels = lv))
  SampProb_vec <- ifelse(N_by > 0, as.numeric(nsel_by) / as.numeric(N_by), 0)
  names(SampProb_vec) <- lv
  probM <- matrix(SampProb_vec, nrow(samp), length(lv), byrow = TRUE)

  # Get sensible starting values from a naive lmer model
  m0 <- lme4::lmer(y ~ t + z + x + x:t + (t | id),
             data = samp,
             REML = FALSE)

  beta0 <- lme4::fixef(m0)
  vc    <- lme4::VarCorr(m0)$id
  sd_b0 <- attr(vc, "stddev")[1]
  sd_b1 <- attr(vc, "stddev")[2]
  rho01 <- attr(vc, "correlation")[1, 2]
  sd_e  <- sigma(m0)

  z_rho   <- log((1 + rho01) / (1 - rho01))
  InitVals <- c(beta0, log(sd_b0), log(sd_b1), z_rho, log(sd_e))

  # Fit ACML
  fit <- acml.lmem2(
    formula.fixed   = y ~ t + z + x + x:t,
    formula.random  = ~ 1 + t,
    data      = samp,
    id        = id,
    w.function = rep(sampling_type, nrow(samp)),
    InitVals   = InitVals,
    cutpoints  = cutM,
    SampProb   = probM,
    Weights    = Weights
  )

  ## 95% Wald CIs (robust)
  est <- fit$coefficients
  V   <- fit$robcov
  se  <- sqrt(diag(V))
  z   <- qnorm(0.975)

  ci_est <- data.frame(
    parameter = c(
      "(Intercept)",
      "t",
      "z",
      "x",
      "x:t",
      "log(sd_(Intercept))",
      "log(sd_t)",
      "z_rho((Intercept),t)",
      "log(sd_eps)"),
    estimate = est,
    se = se,
    lower = est - z*se,
    upper = est + z*se
  )

  return(ci_est)
}
