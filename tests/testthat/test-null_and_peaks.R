test_that("compute_null_distributions returns proportions summing to ~1", {
  set.seed(1)
  d <- data.frame(tk = runif(200, 0.1, 1))
  out <- compute_null_distributions(d, tk_col = tk, distribution = "uniform",
                                    n_sim = 1e4, n_max = 20, seed = 1)
  expect_equal(sum(out$normalization_constant), 1, tolerance = 1e-3)
  expect_true(is.data.frame(attr(out, "params")))
})

test_that("find_rk_peaks finds the mode of the 1_1 category", {
  set.seed(2)
  rk <- stats::rnorm(500, 0.5, 0.03)
  out <- find_rk_peaks(rk, peaks = "1_1", plot = FALSE)
  expect_equal(out$peak_x, 0.5, tolerance = 0.03)
})

test_that("test_peaks runs on simulated normalizations", {
  set.seed(3)
  counts <- data.frame(`1_1_on` = rpois(20, 10) + 1, `1_1_off` = rpois(20, 3),
                       check.names = FALSE)
  mk <- function() data.frame(type = NA_character_,
                              bin = c("1_1_on", "1_1_off"),
                              normalization_constant = c(0.1, 0.1))
  res <- suppressWarnings(
    test_peaks(counts, list(mk(), mk()), verbose = FALSE)
  )
  expect_equal(nrow(res$p_values), 2)
  expect_equal(res$summary$n_run, 2)
})

test_that("test_peaks does not test categories without normalization constants", {
  set.seed(4)
  counts <- data.frame(type = rep(c("a", "b"), each = 10),
                       `1_1_on` = rpois(20, 10) + 1, `1_1_off` = rpois(20, 3),
                       `1_25_on` = rpois(20, 2), `1_25_off` = rpois(20, 2),
                       light_on = 1, light_off = 0, check.names = FALSE)
  nd <- data.frame(bin = c("1_1_on", "1_1_off"), normalization_constant = c(0.1, 0.1))

  # without type: 1_25 not tested, light_on/off not taken as a category
  expect_warning(res <- test_peaks(counts, list(nd), verbose = FALSE), "1_25")
  expect_equal(names(res$p_values), c("run", "p_1_1", "p_1_25"))
  expect_false(is.na(res$p_values$p_1_1))
  expect_true(is.na(res$p_values$p_1_25))
  expect_equal(res$summary$n_run, c(1, 0))
  expect_equal(res$summary$test, c("1_1_on > 1_1_off", "1_25_on > 1_25_off"))

  # the summary label follows 'alternative'
  res_less <- suppressWarnings(
    test_peaks(counts, list(nd), peaks = "1_1", alternative = "less", verbose = FALSE)
  )
  expect_equal(res_less$summary$test, "1_1_on < 1_1_off")

  # with type: type "b" missing from norm_df -> no test for group b
  nd_t <- cbind(type = "a", nd)
  expect_warning(res_t <- test_peaks(counts, list(nd_t), type_col = type, verbose = FALSE),
                 "1_1 \\[b\\]")
  expect_false(is.na(res_t$p_values$p_1_1[res_t$p_values$type == "a"]))
  expect_true(is.na(res_t$p_values$p_1_1[res_t$p_values$type == "b"]))
})
