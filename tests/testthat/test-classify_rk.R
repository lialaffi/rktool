test_that("classify_rk assigns on/off labels and returns a factor", {
  out <- classify_rk(c(0.5, 0.3, 1/3, 0.1), n_max = 20)
  expect_s3_class(out, "factor")
  expect_equal(as.character(out), c("1_1_on", "1_2_off", "1_2_on", "1_9_on"))
  expect_true(all(c("1_20_off", "1_1_on", "20_1_on") %in% levels(out)))
})

test_that("classify_rk uses n_max = 5 by default", {
  out <- classify_rk(c(0.5, 0.1))
  expect_equal(levels(out)[c(1, length(levels(out)))], c("1_5_off", "5_1_on"))
  expect_equal(as.character(out), c("1_1_on", NA))   # 0.1 is 1_9: beyond n_max
})

test_that("classify_rk respects peaks and gives NA outside them", {
  out <- classify_rk(c(0.5, 0.3), peaks = c("1_1", "2_1"))
  expect_equal(as.character(out), c("1_1_on", NA))
  expect_equal(levels(out), c("1_1_off", "1_1_on", "2_1_off", "2_1_on"))
})

test_that("classify_rk handles n_max = 1, NA and values outside (0, 1)", {
  out <- classify_rk(c(0.5, 0.41, NA, -1, 2), n_max = 1)
  expect_equal(levels(out), c("1_1_off", "1_1_on"))
  expect_equal(as.character(out), c("1_1_on", "1_1_off", NA, NA, NA))
})

test_that("classify_rk intervals are left-open (a, b]", {
  out <- classify_rk(c(0.4, 0.5, 0.6), peaks = "1_1")
  expect_equal(as.character(out), c(NA, "1_1_on", "1_1_off"))
})

test_that("classify_rk validates input", {
  expect_error(classify_rk(0.5, peaks = "2_3"))
  expect_error(classify_rk(0.5, peaks = c("1_1", "1_1")))
  expect_error(classify_rk(0.5, n_max = 0))
})
