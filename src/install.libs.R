libs <- file.path(R_PACKAGE_DIR, "libs", R_ARCH)
dir.create(libs, recursive = TRUE, showWarnings = FALSE)
for (file in c("symbols.rds", Sys.glob(paste0("*", SHLIB_EXT)))) {
    if (file.exists(file)) {
        file.copy(file, file.path(libs, file))
    }
}

copy_dir <- if (requireNamespace("fs", quietly = TRUE)) {
    fs::dir_copy
} else {
    function(path, new_path) {
        dir.create(new_path, recursive = TRUE, showWarnings = FALSE)
        file.copy(list.files(path, full.names = TRUE), new_path, recursive = TRUE)
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
        copy_dir(path = inst_stan, new_path = "stan")
    }
}

bin <- file.path(R_PACKAGE_DIR, "bin")
if (!file.exists(bin)) {
    dir.create(bin, recursive = TRUE, showWarnings = FALSE)
}
bin_stan <- file.path(bin, "stan")
copy_dir(path = "stan", new_path = bin_stan)

if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    message(
        "NOTE: cmdstanr not found. Stan models will not be compiled ",
        "during installation.\n",
        "Install cmdstanr and reinstall this package to compile Stan models."
    )
} else {
    compile_models <- function(bin_stan) {
        models <- instantiate::stan_package_model_files(path = bin_stan)
        cmdstan_path <- cmdstanr::cmdstan_path()

        message(sprintf("Compiling %d Stan models...", length(models)))
        message("Building precompiled header with first model...")
        cmdstanr::cmdstan_model(
            models[[1]], stanc_options = list("O1"), quiet = TRUE
        )

        if (length(models) > 1L) {
            remaining <- models[-1]
            has_mirai <- requireNamespace("mirai", quietly = TRUE)
            has_purrr_parallel <- requireNamespace("purrr", quietly = TRUE) &&
                "in_parallel" %in% getNamespaceExports("purrr")

            if (has_mirai && has_purrr_parallel) {
                n_workers <- min(
                    parallel::detectCores() - 1L, length(remaining)
                )
                message(sprintf(
                    "Compiling remaining %d models in parallel (%d workers)...",
                    length(remaining), n_workers
                ))
                mirai::daemons(n_workers)
                on.exit(mirai::daemons(0), add = TRUE)
                purrr::walk(
                    remaining,
                    purrr::in_parallel(
                        \(m) {
                            options(cmdstanr_path = cmdstan_path)
                            cmdstanr::cmdstan_model(
                                m, stanc_options = list("O1"), quiet = TRUE
                            )
                        },
                        cmdstan_path = cmdstan_path
                    )
                )
            } else {
                message(sprintf(
                    "Compiling remaining %d models sequentially...",
                    length(remaining)
                ))
                for (m in remaining) {
                    cmdstanr::cmdstan_model(
                        m, stanc_options = list("O1"), quiet = TRUE
                    )
                }
            }
        }
        message("Done compiling Stan models.")
    }

    tryCatch(
        if (requireNamespace("callr", quietly = TRUE)) {
            callr::r(
                func = compile_models,
                args = list(bin_stan = bin_stan),
                show = TRUE,
                stderr = "2>&1"
            )
        } else {
            compile_models(bin_stan)
        },
        error = function(e) {
            message(
                "NOTE: Stan model compilation failed. The package will ",
                "still be installed,\nbut models will need to be compiled ",
                "before use.\nError: ", conditionMessage(e)
            )
        }
    )
}
