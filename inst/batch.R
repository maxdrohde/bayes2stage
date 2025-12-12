suppressPackageStartupMessages({
  library(mirai)
  library(fs)
  library(arrow)
  library(data.table)
  library(janitor)
  library(dplyr)
  library(purrr)
})

source("simulation_config.R")

# Get true parameter values for ALL simulation settings
all_params <- get_all_true_parameter_values()

# ==============================================================================
# MERGE PARQUET FILES FIRST (before parallel processing)
# ==============================================================================

processed_dir <- "processed_data"
fs::dir_create(processed_dir)
cached_file <- file.path(processed_dir, "combined_data.parquet")

if (fs::file_exists(cached_file)) {
  message(glue::glue("Cached data already exists at {cached_file}"))
} else {
  message("Merging all parquet files before parallel processing...")

  if (!fs::dir_exists("./results")) {
    stop("Directory './results' does not exist.")
  }

  # Load main results
  message("Reading parquet files from ./results...")
  res <-
    arrow::open_dataset("./results") |>
    dplyr::collect() |>
    janitor::clean_names() |>
    dplyr::select(
      sim_setting, parameter, estimate = mean, se = sd, type,
      rhat, ess_bulk = n_eff_bulk, ess_tail = n_eff_tail,
      divergent_transitions, max_treedepth_exceeded, ebfmi_min
    ) |>
    dplyr::mutate(
      type = clean_type_names(type),
      type = factor(type, levels = TYPE_LEVELS)
    )
  message(glue::glue("Loaded {nrow(res)} rows from ./results"))

  # Load ACML results (if available)
  acml_path <- "./results_acml"
  has_acml <- fs::dir_exists(acml_path)
  if (has_acml) {
    has_acml <- length(fs::dir_ls(acml_path, glob = "*.parquet", recurse = TRUE)) > 0
  }

  if (has_acml) {
    message("Reading parquet files from ./results_acml...")
    res_acml_raw <-
      arrow::open_dataset(acml_path) |>
      dplyr::collect() |>
      janitor::clean_names()

    # Filter out error rows (those with 'message' column populated or missing 'variable')
    n_total <- nrow(res_acml_raw)

    # Build error mask safely (columns may not exist if all success/all error)
    has_error <- rep(FALSE, n_total)
    if ("message" %in% names(res_acml_raw)) {
      has_error <- has_error | !is.na(res_acml_raw$message)
    }
    if ("variable" %in% names(res_acml_raw)) {
      has_error <- has_error | is.na(res_acml_raw$variable)
    } else {
      # No 'variable' column means all rows are errors
      has_error <- rep(TRUE, n_total)
    }

    n_errors <- sum(has_error)

    if (n_errors > 0) {
      message(glue::glue(
        "Filtered out {n_errors} ACML error rows ({round(100 * n_errors / n_total, 2)}% of total)"
      ))
    }

    # Check if we have valid data to process
    n_valid <- n_total - n_errors
    if (n_valid == 0 || !"variable" %in% names(res_acml_raw)) {
      message("No valid ACML results to include.")
      res_acml <- dplyr::tibble()
    } else {
      res_acml <-
        res_acml_raw |>
        dplyr::filter(!has_error) |>
        dplyr::mutate(type = "ACML ODS") |>
        dplyr::select(sim_setting = sim_setting, parameter = variable, estimate = mean, se = sd, type) |>
        dplyr::mutate(
          # Map ACML parameter names to match main results
          parameter = dplyr::case_when(
            parameter == "beta_z" ~ "beta[1]",
            TRUE ~ parameter
          ),
          type = factor(type, levels = TYPE_LEVELS),
          # ACML doesn't have MCMC diagnostics
          rhat = NA_real_,
          ess_bulk = NA_real_,
          ess_tail = NA_real_,
          divergent_transitions = NA_integer_,
          max_treedepth_exceeded = NA_integer_,
          ebfmi_min = NA_real_
        )
    }
  } else {
    res_acml <- dplyr::tibble()
  }

  # Combine and cache
  full_data <- dplyr::bind_rows(res, res_acml)
  message(glue::glue("Saving merged data to {cached_file}..."))
  arrow::write_parquet(full_data, cached_file)
  message(glue::glue("Done merging. Total rows: {nrow(full_data)}"))
}

# ==============================================================================
# PARALLEL PROCESSING
# ==============================================================================

# Use available cores minus one, at least 1
n_cores <- tryCatch(
  {
    cores <- parallel::detectCores(logical = FALSE)
    max(1L, cores - 1L)
  },
  error = function(e) 1L
)

message(glue::glue("Starting parallel analysis with {n_cores} cores..."))

# Filter all_params to only include (sim_setting, parameter) combos that exist in data
cached_data <- arrow::read_parquet(cached_file)
available_combos <- cached_data |>
  dplyr::distinct(sim_setting, parameter)
rm(cached_data) # Free memory

original_count <- nrow(all_params)
all_params <- all_params |>
  dplyr::semi_join(available_combos, by = c("sim_setting", "parameter"))

if (nrow(all_params) < original_count) {
  message(glue::glue(
    "Note: Filtered from {original_count} to {nrow(all_params)} parameter combinations based on available data"
  ))
}

if (nrow(all_params) == 0) {
  stop("No matching (sim_setting, parameter) combinations found in data!")
}

mirai::daemons(n_cores)

message(glue::glue("Generating {nrow(all_params)} analysis plots..."))
tasks <- mirai::mirai_map(
  .x = seq_len(nrow(all_params)),
  .f = function(idx, all_params) {
    row <- all_params[idx, ]
    sim_setting <- row$sim_setting
    p_name <- row$parameter
    p_val <- row$true_value
    command_str <- paste(
      "Rscript analyze-simulation-results.R",
      shQuote(p_name), p_val, sim_setting
    )
    message(glue::glue("Running: {command_str}"))
    exit_code <- system(command_str)
    list(exit_code = exit_code, command = command_str)
  },
  .args = list(all_params = all_params)
)[.progress, .stop]

results <- tasks[]

# Shut down daemons explicitly (on.exit doesn't work when sourced interactively)
mirai::daemons(0)

failed <- vapply(results, function(res) is.list(res) && res$exit_code != 0, logical(1))
if (any(failed)) {
  bad <- results[failed][[1]]
  stop("Command failed: ", bad$command, " (exit code ", bad$exit_code, ")")
}

# ==============================================================================
# MERGE ALL FIGURES
# ==============================================================================

cat("\n")
message("Merging all figures into a single PDF...")

suppressPackageStartupMessages(library(magick))

plot_files <- glue::glue(
  "plots/setting_{all_params$sim_setting}/simulation_plots_{all_params$parameter}.png"
)
plot_files <- plot_files[fs::file_exists(plot_files)]

if (length(plot_files) > 0) {
  images <- image_read(plot_files)
  image_write(images, "plots/all_parameters_combined.pdf", format = "pdf")
  message(glue::glue("Combined {length(plot_files)} figures into plots/all_parameters_combined.pdf"))
} else {
  warning("No plot files found to merge.")
}

cat("\n")
message("Creating per-parameter PDFs across settings...")

unique_params <- unique(all_params$parameter)
fs::dir_create("plots/by_parameter")

for (param in unique_params) {
  param_files <- glue::glue(
    "plots/setting_{all_params$sim_setting[all_params$parameter == param]}/simulation_plots_{param}.png"
  )
  param_files <- param_files[fs::file_exists(param_files)]

  if (length(param_files) > 0) {
    images <- image_read(param_files)
    # Clean parameter name for filename (replace brackets)
    clean_param <- gsub("[][]", "", param)
    output_file <- glue::glue("plots/by_parameter/{clean_param}_across_settings.pdf")
    image_write(images, output_file, format = "pdf")
  }
}

message(glue::glue("Created {length(unique_params)} per-parameter PDFs in plots/by_parameter/"))

cat("\n")
message("Generating data sample plots...")
source("plot_data_samples.R")
