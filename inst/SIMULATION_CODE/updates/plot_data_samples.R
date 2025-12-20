################################################################################
# Plot Data Samples
#
# This script generates 50 datasets using the same data generating mechanism
# as run_sim.R and creates a combined PDF of all plots.
################################################################################

################################################################################
# Packages
################################################################################

suppressWarnings(suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
}))

################################################################################
# Shared Configuration
################################################################################

source("simulation_config.R")

################################################################################
# Setup
################################################################################

# Number of datasets to generate per setting
n_samples <- DEFAULT_N_DATA_SAMPLES

# Output directory for individual plots
base_plots_dir <- "plots/data_samples"
fs::dir_create(base_plots_dir)

################################################################################
# Simulation scenario grid
################################################################################

grid <- get_simulation_grid()
n_settings <- nrow(grid)

################################################################################
# Generate plots for all settings
################################################################################

# Set up parallel workers (once, reuse for all settings)
n_cores <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
mirai::daemons(n_cores)

all_plot_files <- list()

for (i in seq_len(n_settings)) {
  params <- dplyr::slice(grid, i)
  current_N <- params[["N"]]
  current_frac <- params[["sampling_fraction"]]

  # Create setting-specific directory
  setting_dir <- glue::glue("{base_plots_dir}/setting_{i}")
  fs::dir_create(setting_dir)

  message(glue::glue("Generating {n_samples} data plots for setting {i} (N={current_N})..."))

  # Generate plots in parallel using mirai_map
  tasks <- mirai::mirai_map(
    .x = seq_len(n_samples),
    .f = \(j, i, params, setting_dir) {
      # Source config to get generate_simulation_data function
      source("simulation_config.R", local = TRUE)

      # Set seed using Cantor mapping (same as run_sim.R)
      seed <- bayes2stage::cantor_seed(i, j)
      set.seed(seed)

      # Generate data
      df <- generate_simulation_data(params)

      # Get N and subset_size for title
      N <- params[["N"]]
      sampling_fraction <- params[["sampling_fraction"]]
      subset_size <- min(DEFAULT_PLOT_SUBSET_SIZE, nrow(df))

      # Create plot
      p <- bayes2stage::plot_data(df, subset_size = subset_size) +
        patchwork::plot_annotation(
          title = glue::glue("Setting {i}, Dataset {j} (N={N}, sampled={as.integer(N * sampling_fraction)})"),
          theme = ggplot2::theme(
            plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5)
          )
        )

      # Save plot
      plot_file <- glue::glue("{setting_dir}/data_sample_{stringr::str_pad(j, 2, pad = '0')}.png")

      ggplot2::ggsave(
        plot = p,
        filename = plot_file,
        width = 12,
        height = 8,
        units = "in",
        dpi = 150
      )

      as.character(plot_file)
    },
    .args = list(i = i, params = params, setting_dir = setting_dir)
  )[.progress]

  # Collect results for this setting
  plot_files <- unlist(tasks[])
  all_plot_files[[i]] <- plot_files

  # Combine into per-setting PDF
  message(glue::glue("Combining setting {i} plots into PDF..."))
  images <- magick::image_read(plot_files)
  output_pdf <- glue::glue("{base_plots_dir}/data_samples_setting_{i}.pdf")
  magick::image_write(images, output_pdf, format = "pdf")
  message(glue::glue("Created {output_pdf}"))
}

# Shut down daemons
mirai::daemons(0)

################################################################################
# Combine all settings into one master PDF
################################################################################

message("Combining all settings into master PDF...")

all_files <- unlist(all_plot_files)
images <- magick::image_read(all_files)
master_pdf <- "plots/data_samples_all_settings.pdf"
magick::image_write(images, master_pdf, format = "pdf")

message(glue::glue("Done! Created {length(all_files)} plots across {n_settings} settings."))
message(glue::glue("Master PDF: {master_pdf}"))
