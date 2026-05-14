test_that("generate_data returns correct structure", {
    df <- generate_data(N = 100L, M = 5L, x_dist = "normal")

    expect_equal(nrow(df), 100L * 5L)
    expect_equal(length(unique(df$id)), 100L)

    expect_true(all(c("y", "t", "x", "z", "id") %in% names(df)))

    expect_equal(min(df$t), 0)
    expect_equal(max(df$t), 1)

    x_by_id <- tapply(df$x, df$id, function(v) length(unique(v)))
    expect_true(all(x_by_id == 1L))
})

test_that("all distribution types work", {
    expect_no_error(generate_data(x_dist = "normal", N = 50L))
    expect_no_error(generate_data(x_dist = "bernoulli", N = 50L))
    expect_no_error(generate_data(x_dist = "poisson", N = 50L))
    expect_no_error(generate_data(x_dist = "binomial", N = 50L, x_size = 10L))
    expect_no_error(generate_data(x_dist = "negative_binomial", N = 50L, x_disp_param = 1L))
    expect_no_error(generate_data(x_dist = "beta_binomial", N = 50L, x_disp_param = 1L, x_size = 10L))
})

test_that("binomial distributions require size parameter", {
    expect_error(generate_data(x_dist = "binomial", N = 50L))
    expect_error(generate_data(x_dist = "beta_binomial", N = 50L, x_disp_param = 1L))
})

test_that("dispersion parameter required for count distributions", {
    expect_error(generate_data(x_dist = "negative_binomial", N = 50L))
    expect_error(generate_data(x_dist = "beta_binomial", N = 50L, x_size = 10L))
})

test_that("invalid distribution type errors", {
    expect_error(generate_data(x_dist = "invalid", N = 50L))
})

test_that("bernoulli distribution works", {
    df <- generate_data(N = 100L, x_dist = "bernoulli")
    expect_true(all(df$x %in% c(0L, 1L)))
})

test_that("direction = 'z_given_x' generates data correctly", {
    df <- generate_data(
        N = 100L,
        M = 5L,
        direction = "z_given_x",
        x_dist = "bernoulli",
        x_prevalence = 0.3
    )

    expect_equal(nrow(df), 100L * 5L)
    expect_true(all(df$x %in% c(0L, 1L)))

    x_by_id <- tapply(df$x, df$id, function(v) length(unique(v)))
    expect_true(all(x_by_id == 1L))
})

test_that("direction = 'z_given_x' requires bernoulli x_dist", {
    expect_error(
        generate_data(direction = "z_given_x", x_dist = "normal", N = 50L),
        "must be.*bernoulli"
    )
})

test_that("x_prevalence validation works", {
    expect_error(
        generate_data(direction = "z_given_x", x_dist = "bernoulli", x_prevalence = 0, N = 50L),
        "x_prevalence"
    )
    expect_error(
        generate_data(direction = "z_given_x", x_dist = "bernoulli", x_prevalence = 1, N = 50L),
        "x_prevalence"
    )
})

test_that("variable M works", {
    df <- generate_data(N = 100L, M = c(3L, 4L, 5L, 6L))

    obs_per_subject <- table(df$id)
    expect_true(all(obs_per_subject %in% c(3L, 4L, 5L, 6L)))
    expect_true(length(unique(obs_per_subject)) > 1L)
})

test_that("time is scaled 0-1 with variable M", {
    df <- generate_data(N = 50L, M = c(3L, 5L, 7L))

    time_range_by_id <- tapply(df$t, df$id, function(v) c(min(v), max(v)))
    expect_true(all(sapply(time_range_by_id, function(r) r[1] == 0 && r[2] == 1)))
})

test_that("M must be >= 2", {
    expect_error(generate_data(N = 50L, M = 1L), "M.*>= 2")
    expect_error(generate_data(N = 50L, M = c(1L, 3L, 5L)), "M.*>= 2")
})
