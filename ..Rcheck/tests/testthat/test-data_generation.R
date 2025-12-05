test_that("all the distribution types work", {
  generate_data(x_dist = "normal", N = 500)
  generate_data(x_dist = "poisson", N = 500)
  generate_data(x_dist = "binomial", N = 500, x_size = 1)
  generate_data(x_dist = "negative_binomial", N = 500, x_disp_param = 1)
  generate_data(x_dist = "beta_binomial", N = 500, x_disp_param = 1, x_size = 1)
})

