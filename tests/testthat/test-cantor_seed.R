test_that("cantor_seed returns correct known values", {
    expect_equal(cantor_seed(1L, 1L), 4L)
    expect_equal(cantor_seed(1L, 2L), 8L)
    expect_equal(cantor_seed(2L, 1L), 7L)
    expect_equal(cantor_seed(5L, 10L), 130L)
    expect_equal(cantor_seed(3L, 3L), 24L)
})

test_that("cantor_seed returns integer type", {
    result <- cantor_seed(1L, 1L)
    expect_type(result, "integer")
})

test_that("cantor_seed produces unique values for all pairs", {
    seeds <- integer(100L)
    idx <- 1L
    for (i in 1:10) {
        for (j in 1:10) {
            seeds[idx] <- cantor_seed(i, j)
            idx <- idx + 1L
        }
    }
    expect_equal(length(unique(seeds)), 100L)
})

test_that("cantor_seed is asymmetric: (i,j) != (j,i) for i != j", {
    expect_false(cantor_seed(1L, 2L) == cantor_seed(2L, 1L))
    expect_false(cantor_seed(3L, 7L) == cantor_seed(7L, 3L))
    expect_false(cantor_seed(10L, 1L) == cantor_seed(1L, 10L))
})

test_that("cantor_seed errors when i < 1", {
    expect_error(cantor_seed(0L, 1L), "must be >= 1")
    expect_error(cantor_seed(-5L, 1L), "must be >= 1")
})

test_that("cantor_seed errors when j < 1", {
    expect_error(cantor_seed(1L, 0L), "must be >= 1")
    expect_error(cantor_seed(1L, -3L), "must be >= 1")
})

test_that("cantor_seed errors on overflow", {
    # cantor_seed(65535, 1) = 2147516417 > .Machine$integer.max
    expect_error(cantor_seed(65535L, 1L), "overflow")
    expect_error(cantor_seed(32768L, 32768L), "overflow")
    # Extreme case where i + j would overflow integer range
    expect_error(cantor_seed(2000000000L, 2000000000L), "overflow")
})

test_that("cantor_seed succeeds just below overflow boundary", {
    # cantor_seed(65534, 1) = 2147450881, under .Machine$integer.max
    expect_no_error(cantor_seed(65534L, 1L))
    expect_no_error(cantor_seed(32767L, 32767L))
})
