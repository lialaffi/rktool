test_that(".rk_bin_widths returns the width of the rk regions", {
  w_on  <- rktool:::.rk_bin_widths("1_1_on")
  w_off <- rktool:::.rk_bin_widths("1_1_off")
  expect_equal(w_on, 1/9)
  expect_equal(w_on + w_off, 0.2)
  expect_true(is.na(rktool:::.rk_bin_widths("foo")))
})

test_that("normalize_rk_counts works with method = 'width' (wide)", {
  b <- data.frame(`1_1_on` = c(2, 4), `1_1_off` = c(1, 2), check.names = FALSE)
  out <- normalize_rk_counts(b, method = "width")
  expect_equal(out[["1_1_on_norm"]], c(2, 4) / (1/9))
})

test_that("normalize_rk_counts works with a null table (long)", {
  b <- data.frame(bin = c("1_1_on", "1_1_off"), n = c(10, 5))
  norm_df <- data.frame(type = NA_character_, bin = c("1_1_on", "1_1_off"),
                        normalization_constant = c(0.5, 0.25))
  out <- normalize_rk_counts(b, norm_df, format = "long", bin_col = bin)
  expect_equal(out$n_norm, c(20, 20))
  expect_equal(out$norm, c(0.5, 0.25))
})

test_that("normalize_rk_counts gives NA for bins missing from norm_df (wide)", {
  b <- data.frame(type = c("a", "b"), `1_1_on` = c(2, 4), `1_25_on` = c(1, 1),
                  check.names = FALSE)
  norm_df <- data.frame(type = c("a", "b"), bin = "1_1_on",
                        normalization_constant = c(0.5, 0.25))
  bins <- c("1_1_on", "1_25_on")

  # without type_col (one constant per bin)
  out <- normalize_rk_counts(b, norm_df[1, c("bin", "normalization_constant")],
                             bin_cols = bins)
  expect_equal(out[["1_1_on_norm"]], c(4, 8))
  expect_equal(out[["1_25_on_norm"]], c(NA_real_, NA_real_))

  # with type_col
  out_t <- normalize_rk_counts(b, norm_df, type_col = type, bin_cols = bins)
  expect_equal(out_t[["1_1_on_norm"]], c(4, 16))
  expect_equal(out_t[["1_25_on_norm"]], c(NA_real_, NA_real_))
})

test_that("normalize_rk_counts requires norm_df with method = 'null'", {
  b <- data.frame(`1_1_on` = 1, check.names = FALSE)
  expect_error(normalize_rk_counts(b))
})
