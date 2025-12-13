libs <- file.path(R_PACKAGE_DIR, "libs", R_ARCH)
dir.create(libs, recursive = TRUE, showWarnings = FALSE)
for (file in c("symbols.rds", Sys.glob(paste0("*", SHLIB_EXT)))) {
  if (file.exists(file)) {
    file.copy(file, file.path(libs, file))
  }
}
inst_stan <- file.path("..", "inst", "stan")
if (dir.exists(inst_stan)) {
  warning(
    "Stan models in inst/stan/ are deprecated in {instantiate} ",
    ">= 0.0.4.9001 (2024-01-03). Please put them in src/stan/ instead."
  )
  if (file.exists("stan")) {
    warning("src/stan/ already exists. Not copying models from inst/stan/.")
  } else {
    message("Copying inst/stan/ to src/stan/.")
    fs::dir_copy(path = inst_stan, new_path = "stan")
  }
}
bin <- file.path(R_PACKAGE_DIR, "bin")
if (!file.exists(bin)) {
  dir.create(bin, recursive = TRUE, showWarnings = FALSE)
}
bin_stan <- file.path(bin, "stan")
fs::dir_copy(path = "stan", new_path = bin_stan)
callr::r(
  func = function(bin_stan) {
    models <- instantiate::stan_package_model_files(path = bin_stan)
    cmdstan_path <- cmdstanr::cmdstan_path()

    # Compile first model sequentially to build precompiled header
    message(sprintf("Compiling %d Stan models...", length(models)))
    message("Building precompiled header with first model...")
    cmdstanr::cmdstan_model(models[[1]], stanc_options = list("O1"), quiet = TRUE)

    if (length(models) > 1L) {
      remaining <- models[-1]
      n_workers <- min(parallel::detectCores() - 1L, length(remaining))
      message(sprintf("Compiling remaining %d models in parallel using %d daemons...",
                      length(remaining), n_workers))
      mirai::daemons(n_workers)
      purrr::walk(
        remaining,
        purrr::in_parallel(
          \(m) {
            options(cmdstanr_path = cmdstan_path)
            cmdstanr::cmdstan_model(m, stanc_options = list("O1"), quiet = TRUE)
          },
          cmdstan_path = cmdstan_path
        )
      )
      mirai::daemons(0)
    }
    message("Done compiling Stan models.")
  },
  args = list(bin_stan = bin_stan),
  show = TRUE,
  stderr = "2>&1"
)
