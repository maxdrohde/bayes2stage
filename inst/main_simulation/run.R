library(furrr)

future::plan(multisession, workers = 12)


# Run simulation
furrr::future_walk(1:56,
                   ~ bayes2stage:::simulate(.x, ITER = 2000),
                   .options = furrr_options(seed=TRUE, chunk_size = 2L),
                   .progress = FALSE)
