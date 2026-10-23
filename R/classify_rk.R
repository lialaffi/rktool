#' Parse a rhythmic peak label (format "1_n", "n_1" or "1_1")
#'
#' Internal helper for \code{\link{classify_rk}}: validates and interprets a
#' single peak label supplied via \code{peaks}.
#'
#' @param peak A string in the format \code{"1_n"}, \code{"n_1"} or
#'   \code{"1_1"} (\code{n} an integer >= 1).
#' @return A list with elements \code{n} (integer) and \code{side} (one of
#'   \code{"left"}, \code{"right"}, \code{"center"}).
#' @keywords internal
.parse_rk_peak <- function(peak) {
  parts <- strsplit(peak, "_", fixed = TRUE)[[1]]
  nums <- suppressWarnings(as.integer(parts))
  if (length(parts) != 2 || anyNA(nums) || any(nums < 1)) {
    stop("Malformed peak: '", peak, "'. Expected format: '1_n', 'n_1' or '1_1' (n an integer >= 1).")
  }
  a <- nums[1]; b <- nums[2]
  if (a == 1 && b == 1) {
    list(n = 1, side = "center")
  } else if (a == 1 && b > 1) {
    list(n = b, side = "left")
  } else if (b == 1 && a > 1) {
    list(n = a, side = "right")
  } else {
    stop("Malformed peak: '", peak, "'. Only '1_n' or 'n_1' categories exist (one of the two parts must be 1).")
  }
}

#' Compute breakpoints and labels for a single rhythmic category
#'
#' Internal helper for \code{\link{classify_rk}}, \code{\link{find_rk_peaks}}
#' and \code{\link{normalize_rk_counts}}. Breakpoints are computed with
#' \code{f(n, c) = 1 / (n + c)} (categories \code{"1_n"}) and its mirror
#' \code{g(n, c) = 1 - f(n, c)} (categories \code{"n_1"}).
#'
#' @param n Integer, the category index (e.g. 3 for "1_3" or "3_1").
#' @param side One of \code{"left"} (category "1_n"), \code{"right"}
#'   (category "n_1"), \code{"center"} (category "1_1").
#' @return A list with elements \code{breaks} (5 increasing numeric values)
#'   and \code{labels} (4 strings: off1, on1, on2, off2).
#' @keywords internal
.rk_category_breaks <- function(n, side) {
  f <- function(n, c) 1 / (n + c)
  g <- function(n, c) 1 - f(n, c)
  offs <- c(1.5, 1.25, 1.0, 0.75, 0.5)
  if (side == "center") {
    breaks <- c(f(1, 1.5), f(1, 1.25), f(1, 1.0), g(1, 1.25), g(1, 1.5))
    labels <- paste0("1_1", c("_off1", "_on1", "_on2", "_off2"))
  } else if (side == "left") {
    breaks <- f(n, offs)
    labels <- paste0("1_", n, c("_off1", "_on1", "_on2", "_off2"))
  } else { # right
    breaks <- g(n, rev(offs))
    labels <- paste0(n, "_1", c("_off1", "_on1", "_on2", "_off2"))
  }
  list(breaks = breaks, labels = labels)
}


#' Classify rk values into rhythmic bins
#'
#' Assigns each rhythmic ratio \code{rk} to the "on" or "off" bin of a
#' rhythmic category, from \code{1_n_max} to \code{1_1} up to
#' \code{n_max_1}, or of an explicit (possibly non-contiguous) set of
#' categories given in \code{peaks}.
#'
#' Each category \code{1_n} is centered on \code{1/(n+1)}: its "on" bin is
#' the central part and its "off" bin the two outer parts. The \code{n_1}
#' categories are the mirror image (\code{1 - rk}), and \code{1_1} is
#' centered on 0.5. Bins are returned as labels such as \code{"1_5_off"} /
#' \code{"1_5_on"}.
#'
#' When \code{peaks} is \code{NULL}, the categories are contiguous and cover
#' the whole interval \code{(1/(n_max+1.5), 1 - 1/(n_max+1.5))}; more
#' extreme values become \code{NA}. When \code{peaks} is given, \code{rk}
#' values falling in none of the requested categories become \code{NA}
#' (e.g. \code{c("1_1", "1_2", "2_1", "4_1")} skips \code{"3_1"}).
#'
#' @param rk Numeric vector of rhythmic ratios, values expected in (0, 1).
#' @param peaks Optional character vector of categories to classify, in the
#'   format \code{"1_n"}, \code{"n_1"} or \code{"1_1"} (e.g.
#'   \code{c("1_1", "1_2", "2_1", "4_1")}). If \code{NULL} (default), all
#'   categories from \code{1_n_max} to \code{n_max_1} are used and
#'   \code{n_max} applies; otherwise \code{n_max} is ignored. Malformed or
#'   duplicated labels raise an error.
#' @param n_max Integer, the most extreme category to include when
#'   \code{peaks} is \code{NULL} (default 5). \code{rk} values more extreme
#'   than \code{1_n_max} or \code{n_max_1} become \code{NA}.
#' @return A factor of the same length as \code{rk}, with levels such as
#'   \code{"1_5_off"}, \code{"1_1_on"}, \code{"5_1_off"}. Level order
#'   follows the order of \code{peaks} (or \code{1_n_max -> 1_1 -> n_max_1}
#'   when \code{peaks} is \code{NULL}).
#'
#' @details
#' Breakpoints of category \code{1_n} are \code{1/(n + c)} for \code{c} in
#' \code{c(1.5, 1.25, 1, 0.75, 0.5)}: the "on" bin is
#' \code{(1/(n+1.25), 1/(n+0.75)]}, the "off" bin
#' \code{(1/(n+1.5), 1/(n+1.25)]} plus \code{(1/(n+0.75), 1/(n+0.5)]}.
#' Category \code{n_1} uses \code{1 - 1/(n + c)}; \code{1_1} uses the
#' \code{1_n} formula below 0.5 and the \code{n_1} formula above.
#'
#' @examples
#' # on the example dataset
#' y <- calculate_tk(rk_example, onset_col = onset, by = c(file_id, ind_id))
#' y <- calculate_rk(y, tk_col = tk, by = c(file_id, ind_id))
#' y$bins <- classify_rk(y$rk)
#' table(y$bins)
#'
#' @export
classify_rk <- function(rk, peaks = NULL, n_max = 5) {
  peaks <- .rk_peaks(peaks, n_max)

  defs <- lapply(peaks, function(pk) {
    p <- .parse_rk_peak(pk)
    .rk_category_breaks(p$n, p$side)
  })

  # each sub-bin is the interval (lo, hi]
  lo  <- unlist(lapply(defs, function(d) d$breaks[-5]))
  hi  <- unlist(lapply(defs, function(d) d$breaks[-1]))
  lab <- unlist(lapply(defs, `[[`, "labels"))

  # Collapse sub-bins: off1/off2 -> off, on1/on2 -> on
  collapsed <- sub("[12]$", "", lab)
  lev <- unique(collapsed)
  code <- match(collapsed, lev)

  # elementary intervals between all the (sorted, unique) breakpoints: each
  # one is either exactly a sub-bin or a gap between non-contiguous peaks
  brk <- sort(unique(c(lo, hi)))
  m <- match(brk[-length(brk)], lo)
  elem_code <- code[m]
  elem_code[is.na(m) | hi[m] != brk[-1]] <- NA_integer_

  # findInterval(left.open = TRUE): i means rk in (brk[i], brk[i+1]];
  # 0 and length(brk) (outside all bins) index past elem_code and give NA
  i <- findInterval(rk, brk, left.open = TRUE)
  i[i == 0L] <- NA_integer_

  # build the factor directly from the integer codes (no string handling)
  structure(elem_code[i], levels = lev, class = "factor")
}
