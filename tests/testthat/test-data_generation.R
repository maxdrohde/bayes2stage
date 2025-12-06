test_that("generate_data returns correct structure", {
  df <- generate_data(N = 100, M = 5, x_dist = "normal")

  # Check dimensions

  expect_equal(nrow(df), 100 * 5)
  expect_equal(length(unique(df$id)), 100)

  # Check required columns exist
  expect_true(all(c("y", "t", "x", "z", "id") %in% names(df)))

  # Check time is scaled 0-1
  expect_equal(min(df$t), 0)
  expect_equal(max(df$t), 1)

  # Check x is constant within subject
  x_by_id <- tapply(df$x, df$id, function(v) length(unique(v)))
  expect_true(all(x_by_id == 1))
})

test_that("all distribution types work", {
  expect_no_error(generate_data(x_dist = "normal", N = 50))
  expect_no_error(generate_data(x_dist = "poisson", N = 50))
  expect_no_error(generate_data(x_dist = "binomial", N = 50, x_size = 10))
  expect_no_error(generate_data(x_dist = "negative_binomial", N = 50, x_disp_param = 1))
  expect_no_error(generate_data(x_dist = "beta_binomial", N = 50, x_disp_param = 1, x_size = 10))
})

test_that("binomial distributions require size parameter", {
  expect_error(generate_data(x_dist = "binomial", N = 50))
  expect_error(generate_data(x_dist = "beta_binomial", N = 50, x_disp_param = 1))
})

test_that("dispersion parameter required for count distributions", {
  expect_error(generate_data(x_dist = "negative_binomial", N = 50))
  expect_error(generate_data(x_dist = "beta_binomial", N = 50, x_size = 10))
})

test_that("invalid distribution type errors",
{
  expect_error(generate_data(x_dist = "invalid", N = 50))
})
