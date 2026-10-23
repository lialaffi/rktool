#' Find rk density peaks for one or more rhythmic categories
#'
#' For each requested rhythmic category (the same categories used by
#' \code{\link{classify_rk}}), finds the location of the mode of the kernel
#' density of \code{rk} within the value range of that category (its "on"
#' and "off" bins together).
#'
#' The density is estimated once on the whole \code{rk} vector (smoother
#' than re-fitting it on each small subset) and its maximum is then taken
#' within each category's range. Ranges use the same breakpoints as
#' \code{\link{classify_rk}}, so results are directly comparable with a
#' classification of the same \code{rk} vector.
#'
#' @param rk Numeric vector of rhythmic ratios (same input expected by
#'   \code{\link{classify_rk}}). \code{NA} values are dropped before
#'   estimating the density.
#' @param peaks Optional character vector of category labels, in the same
#'   format as \code{\link{classify_rk}}'s \code{peaks} argument (e.g.
#'   \code{c("1_1", "1_2", "2_1", "4_1")}). If \code{NULL} (default), all
#'   categories from \code{1_n_max} to \code{n_max_1} are used. Malformed or
#'   duplicated labels raise an error.
#' @param n_max Most extreme category to include when \code{peaks} is
#'   \code{NULL} (default 5). Ignored if \code{peaks} is specified.
#' @param n_points Number of equally spaced points used to evaluate the
#'   kernel density estimate (passed to \code{density()}'s \code{n}
#'   argument, default 2048).
#' @param bw Bandwidth passed to \code{density()} (default \code{"nrd0"}, R's
#'   default rule).
#' @param min_obs Minimum number of observed \code{rk} values required
#'   inside a category's range for its peak to be reported; categories with
#'   fewer observations get \code{NA} peak values, with a warning (default
#'   2).
#' @param plot Logical, if \code{TRUE} (default) draws the overall density
#'   of \code{rk} with a dashed vertical line at each estimated peak, plus a
#'   legend.
#' @param plot_title Optional plot title (default \code{"rk density"}).
#'
#' @return A data frame with one row per requested category and columns:
#'   \item{peak}{Category label (e.g. \code{"1_1"}, \code{"1_2"}).}
#'   \item{range_min, range_max}{Value range of the category, from
#'     \code{\link{classify_rk}}'s breakpoints.}
#'   \item{n_obs}{Number of observed \code{rk} values falling in the range.}
#'   \item{peak_x}{Estimated \code{rk} value at the density peak (mode)
#'     within the range; \code{NA} if \code{n_obs < min_obs}.}
#'   \item{peak_y}{Density height at \code{peak_x}; \code{NA} together with
#'     \code{peak_x}.}
#'
#' @examples
#' y <- calculate_tk(data_rktool, onset_col = onset, by = c(file_id, ind_id))
#' y <- calculate_rk(y, tk_col = tk, by = c(file_id, ind_id))
#'
#' # peak of a single category
#' find_rk_peaks(y$rk, peaks = "1_1", plot = FALSE)
#'
#' # peaks of all categories up to n_max = 5 (default), with plot
#' find_rk_peaks(y$rk)
#'
#' # peaks of a specific, non-contiguous subset of categories
#' find_rk_peaks(y$rk, peaks = c("1_1", "1_2", "2_1", "4_1"), plot = FALSE)
#'
#' @export
find_rk_peaks <- function(rk, peaks = NULL, n_max = 5, n_points = 2048,
                          bw = "nrd0", min_obs = 2, plot = TRUE, plot_title = NULL) {
  peaks <- .rk_peaks(peaks, n_max)

  rk <- rk[!is.na(rk)]
  if (length(rk) < 2) {
    stop("Need at least 2 non-missing rk values to estimate a density.")
  }

  dens <- stats::density(rk, n = n_points, bw = bw)

  rows <- lapply(peaks, function(pk) {
    p <- .parse_rk_peak(pk)
    brk <- .rk_category_breaks(p$n, p$side)$breaks
    range_min <- brk[1]
    range_max <- brk[5]

    n_obs <- sum(rk > range_min & rk <= range_max)
    in_range <- dens$x >= range_min & dens$x <= range_max

    peak_x <- peak_y <- NA_real_
    if (n_obs < min_obs || !any(in_range)) {
      warning("Category '", pk, "' has fewer than ", min_obs,
              " observed rk values in its range; peak set to NA.")
    } else {
      j <- which(in_range)[which.max(dens$y[in_range])]
      peak_x <- dens$x[j]
      peak_y <- dens$y[j]
    }

    data.frame(peak = pk, range_min = range_min, range_max = range_max,
               n_obs = n_obs, peak_x = peak_x, peak_y = peak_y,
               stringsAsFactors = FALSE)
  })

  result <- do.call(rbind, rows)
  rownames(result) <- NULL

  if (plot) {
    graphics::plot(dens,
                   main = if (is.null(plot_title)) "rk density" else plot_title,
                   xlab = "rk", ylab = "density")
    valid <- !is.na(result$peak_x)
    if (any(valid)) {
      cols <- grDevices::rainbow(sum(valid))
      graphics::abline(v = result$peak_x[valid], col = cols, lty = 2)
      graphics::legend("topright", legend = result$peak[valid], col = cols,
                       lty = 2, bty = "n", cex = 0.7)
    }
  }

  result
}
