test_that("count_rk returns wide and long counts", {
  d <- data.frame(
    id = c("a", "a", "b"),
    bins = factor(c("1_1_on", "1_1_off", "1_1_on"))
  )
  w <- count_rk(d, by = "id", bin_col = "bins")
  expect_equal(w[["1_1_on"]], c(1, 1))
  expect_equal(w[["1_1_off"]], c(1, 0))

  l <- count_rk(d, by = "id", bin_col = "bins", format = "long")
  expect_equal(nrow(l), 4)
  expect_equal(sum(l$n), 3)
})

test_that("count_rk adds the missing on/off bin of an observed category", {
  d <- data.frame(id = c("a", "a", "b"), rk = c(0.5, 0.49, 0.3))
  d$bins <- classify_rk(d$rk, peaks = c("1_1", "1_2", "1_3"))
  w <- count_rk(d, by = id, bin_col = bins)
  # 1_1: only "on" observed; 1_2: only "off" observed; 1_3: never observed
  # columns follow the factor levels (order of 'peaks')
  expect_equal(names(w), c("id", "1_1_off", "1_1_on", "1_2_off", "1_2_on"))
  expect_equal(w[["1_1_off"]], c(0, 0))
  expect_equal(w[["1_2_on"]], c(0, 0))

  l <- count_rk(d, by = id, bin_col = bins, format = "long")
  expect_equal(nrow(l), 8)
  expect_equal(sum(l$n), 3)

  # character bins, no factor levels: same completion, off before on
  d$bins <- as.character(d$bins)
  expect_setequal(names(count_rk(d, by = id, bin_col = bins)),
                  c("id", "1_2_off", "1_2_on", "1_1_off", "1_1_on"))
})

test_that("count_rk accepts bare column names", {
  d <- data.frame(
    id = c("a", "a", "b"),
    sex = c("f", "f", "m"),
    bins = factor(c("1_1_on", "1_1_off", "1_1_on"))
  )
  expect_equal(count_rk(d, by = id, bin_col = bins, keep_cols = sex),
               count_rk(d, by = "id", bin_col = "bins", keep_cols = "sex"))
})
