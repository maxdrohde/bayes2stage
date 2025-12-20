# Using explicit namespacing instead of library() calls

source("simulation_config.R")

# Get true parameter values for ALL simulation settings
all_params <- get_all_true_parameter_values()

cli::cli_h1("Merge Parquet Files")

processed_dir <- "processed_data"
fs::dir_create(processed_dir)
cached_file <- file.path(processed_dir, "combined_data.parquet")

if (fs::file_exists(cached_file)) {
  cli::cli_alert_info("Cached data already exists at {.path {cached_file}}")
} else {
  cli::cli_alert_info("Merging all parquet files before parallel processing...")

  # Load main results (if available)
  has_main <- fs::dir_exists("./results")
  if (has_main) {
    has_main <- length(fs::dir_ls(
      "./results",
      glob = "*.parquet",
      recurse = TRUE
    )) >
      0
  }

  if (has_main) {
    cli::cli_alert("Reading parquet files from {.path ./results}...")
    res_raw <-
      arrow::open_dataset("./results") |>
      dplyr::collect() |>
      janitor::clean_names()

    # Filter out ACML error rows (those with 'message' column populated)
    n_total <- nrow(res_raw)
    has_error <- rep(FALSE, n_total)
    if ("message" %in% names(res_raw)) {
      has_error <- !is.na(res_raw$message)
    }
    n_errors <- sum(has_error)

    if (n_errors > 0) {
      cli::cli_alert_warning(
        "Filtered out {.val {n_errors}} error rows ({round(100 * n_errors / n_total, 2)}% of total)"
      )
    }

    # Handle potential column name differences (posterior::rhat vs rhat)
    # janitor::clean_names() converts posterior::rhat -> posterior_rhat
    if ("posterior_rhat" %in% names(res_raw)) {
      res_raw <- res_raw |>
        dplyr::mutate(
          rhat = dplyr::coalesce(posterior_rhat, rhat),
          ess_bulk = dplyr::coalesce(posterior_ess_bulk, ess_bulk),
          ess_tail = dplyr::coalesce(posterior_ess_tail, ess_tail)
        )
    }

    full_data <-
      res_raw |>
      dplyr::filter(!has_error) |>
      dplyr::select(
        sim_setting,
        sim_iter,
        parameter,
        estimate = mean,
        se = sd,
        type,
        q2_5,
        q97_5,
        rhat,
        ess_bulk,
        ess_tail,
        divergent_transitions,
        max_treedepth_exceeded,
        ebfmi_min
      ) |>
      dplyr::mutate(
        type = clean_type_names(type),
        type = factor(type, levels = TYPE_LEVELS)
      )
    cli::cli_alert_success(
      "Loaded {.val {nrow(full_data)}} rows from {.path ./results}"
    )
  } else {
    cli::cli_abort("No data found in {.path ./results}.")
  }
  cli::cli_alert("Saving merged data to {.path {cached_file}}...")
  arrow::write_parquet(full_data, cached_file)
  cli::cli_alert_success("Done merging. Total rows: {.val {nrow(full_data)}}")
}

cli::cli_h1("Parallel Processing")

# Use available cores minus one, at least 1
n_cores <- tryCatch(
  {
    cores <- parallel::detectCores(logical = FALSE)
    max(1L, cores - 1L)
  },
  error = \(e) 1L
)

cli::cli_alert_info("Starting parallel analysis with {.val {n_cores}} cores...")

# Filter all_params to only include (sim_setting, parameter) combos that exist in data
cached_data <- arrow::read_parquet(cached_file)
available_combos <- cached_data |>
  dplyr::distinct(sim_setting, parameter)
rm(cached_data) # Free memory

original_count <- nrow(all_params)
all_params <- all_params |>
  dplyr::semi_join(available_combos, by = c("sim_setting", "parameter"))

if (nrow(all_params) < original_count) {
  cli::cli_alert_info(
    "Filtered from {.val {original_count}} to {.val {nrow(all_params)}} parameter combinations based on available data"
  )
}

if (nrow(all_params) == 0) {
  cli::cli_abort(
    "No matching (sim_setting, parameter) combinations found in data!"
  )
}

mirai::daemons(n_cores)

cli::cli_alert("Generating {.val {nrow(all_params)}} analysis plots...")
tasks <- mirai::mirai_map(
  .x = seq_len(nrow(all_params)),
  .f = \(idx, all_params) {
    row <- all_params[idx, ]
    sim_setting <- row$sim_setting
    p_name <- row$parameter
    p_val <- row$true_value
    command_str <- paste(
      "Rscript analyze-simulation-results.R",
      shQuote(p_name),
      p_val,
      sim_setting
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

failed <- purrr::map_lgl(results, \(res) is.list(res) && res$exit_code != 0)
if (any(failed)) {
  bad <- results[failed][[1]]
  cli::cli_abort("Command failed: {bad$command} (exit code {bad$exit_code})")
}

cli::cli_h1("Merge Figures")

cli::cli_alert("Merging all figures into a single PDF...")

plot_files <- glue::glue(
  "plots/setting_{all_params$sim_setting}/simulation_plots_{all_params$parameter}.png"
)
plot_files <- plot_files[fs::file_exists(plot_files)]

if (length(plot_files) > 0) {
  images <- magick::image_read(plot_files)
  magick::image_write(
    images,
    "plots/all_parameters_combined.pdf",
    format = "pdf"
  )
  cli::cli_alert_success(
    "Combined {.val {length(plot_files)}} figures into {.path plots/all_parameters_combined.pdf}"
  )
} else {
  cli::cli_alert_warning("No plot files found to merge.")
}

cli::cli_alert("Creating per-parameter PDFs across settings...")

unique_params <- unique(all_params$parameter)
fs::dir_create("plots/by_parameter")

for (param in unique_params) {
  param_files <- glue::glue(
    "plots/setting_{all_params$sim_setting[all_params$parameter == param]}/simulation_plots_{param}.png"
  )
  param_files <- param_files[fs::file_exists(param_files)]

  if (length(param_files) > 0) {
    images <- magick::image_read(param_files)
    clean_param <- gsub("[][]", "", param)
    output_file <- glue::glue(
      "plots/by_parameter/{clean_param}_across_settings.pdf"
    )
    magick::image_write(images, output_file, format = "pdf")
  }
}

cli::cli_alert_success(
  "Created {.val {length(unique_params)}} per-parameter PDFs in {.path plots/by_parameter/}"
)

cli::cli_h1("Data Sample Plots")

cli::cli_alert("Generating data sample plots...")
source("plot_data_samples.R")

cli::cli_alert_success("Batch processing complete!")
