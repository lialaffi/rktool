test_that("calculate_tk works from onsets", {
  d <- data.frame(onset = c(3, 0, 4, 1))
  out <- calculate_tk(d, onset_col = onset)
  expect_equal(out$onset, c(0, 1, 3, 4))
  expect_equal(out$tk, c(NA, 1, 2, 1))
})

test_that("calculate_tk + calculate_rk work from onsets", {
  d <- data.frame(onset = c(0, 1, 3, 4))
  out <- calculate_rk(calculate_tk(d, onset_col = onset), tk_col = tk)
  expect_equal(out$tk, c(NA, 1, 2, 1))
  expect_equal(out$rk, c(NA, 1/3, 2/3, NA))
})

test_that("calculate_rk works by group and with tk_max", {
  d <- data.frame(g = c("a", "a", "a", "b", "b"), tk = c(1, 2, 1, 1, 1))
  out <- calculate_rk(d, tk_col = tk, by = g)
  expect_equal(out$rk, c(1/3, 2/3, NA, 0.5, NA))

  d2 <- data.frame(tk = c(1, 2, 1, 10))
  out2 <- calculate_rk(d2, tk_col = tk, tk_max = 5)
  expect_equal(out2$rk, c(1/3, 2/3, NA, NA))
})

test_that("calculate_rk skips negative tk with a warning", {
  d <- data.frame(tk = c(1, 2, -1, 3, 1))
  expect_warning(out <- calculate_rk(d, tk_col = tk), "1 negative value")
  expect_equal(out$rk, c(1/3, NA, NA, 3/4, NA))
})
