#' Calculate rhythmic ratios (rk) from intervals (tk)
#'
#' Computes, within each temporal sequence (e.g. recording x individual), the
#' rhythmic ratio \code{rk = tk / (tk + tk+1)}, where \code{tk+1} is the
#' interval following \code{tk} in the same sequence (row order of
#' \code{data}, or the order given by \code{order_col}). To start from
#' onsets, compute the intervals first with \code{\link{calculate_tk}}.
#'
#' @param data A data frame.
#' @param tk_col Column with the intervals (bare or string).
#' @param by Optional grouping column(s) (bare, string or \code{c(...)}),
#'   e.g. \code{c(file_id, ind_id)}. If \code{NULL} (default) the whole
#'   dataset is treated as a single group.
#' @param order_col Optional column used to sort rows within each group
#'   (bare or string). If \code{NULL} (default) the existing row order is
#'   kept.
#' @param tk_max Optional threshold on the intervals. If \code{NULL}
#'   (default) no threshold is applied. Otherwise, \code{rk} is not
#'   calculated (\code{NA}) whenever either of the two intervals involved
#'   (\code{tk} or \code{tk+1}) is greater than \code{tk_max}.
#' @param rk_name Name of the new rk column (default \code{"rk"}).
#'
#' @return An ungrouped data frame (tibble) with all the original columns
#'   plus the \code{rk} column. The \code{rk} of the last interval of each
#'   group (no following interval) is \code{NA}. \code{rk} is also \code{NA}
#'   when \code{tk + tk+1} is not positive (e.g. simultaneous onsets), when
#'   \code{tk} or \code{tk+1} is negative (with a warning reporting how many
#'   negative intervals were found), or when \code{tk} or \code{tk+1}
#'   exceeds \code{tk_max}.
#'
#' @details
#' \code{rk} is always computed between consecutive rows of the same group
#' after ordering: rows with \code{NA} intervals in the middle of a
#' sequence are not removed, so they break the local ratio (resulting
#' \code{NA}). Remove them beforehand if you want to compute rk across them.
#'
#' @examples
#' # starting from onsets: compute the intervals first
#' y <- calculate_tk(data_rktool, onset_col = onset, by = c(file_id, ind_id))
#' y <- calculate_rk(y, tk_col = tk, by = c(file_id, ind_id))
#' summary(y$rk)
#'
#' # ignore ratios involving long pauses (intervals > 5 s)
#' y2 <- calculate_rk(y, tk_col = tk, by = c(file_id, ind_id), tk_max = 5)
#'
#' @export
calculate_rk <- function(data, tk_col, by = NULL, order_col = NULL,
                         tk_max = NULL, rk_name = "rk") {
  stopifnot(is.data.frame(data))

  tk_col     <- .rk_select(rlang::enquo(tk_col), data, "tk_col", one = TRUE, numeric = TRUE)
  by_cols    <- .rk_select(rlang::enquo(by), data, "by")
  order_cols <- .rk_select(rlang::enquo(order_col), data, "order_col")

  if (rk_name %in% names(data)) {
    warning("Column '", rk_name, "' already exists and will be overwritten.")
  }

  n_neg <- sum(data[[tk_col]] < 0, na.rm = TRUE)
  if (n_neg > 0) {
    warning(n_neg, " negative value(s) in '", tk_col, "'. rk is NA for the ratios involving them.")
  }

  if (length(c(by_cols, order_cols)) > 0) {
    data <- dplyr::arrange(data, dplyr::across(dplyr::all_of(c(by_cols, order_cols))))
  }

  data <- dplyr::group_by(data, dplyr::across(dplyr::all_of(by_cols)))

  # rk = tk / (tk + tk+1); NA for the last row of a group, if the sum is not > 0
  # or if either interval is negative
  data <- dplyr::mutate(
    data,
    !!rk_name := {
      t1 <- .data[[tk_col]]
      t2 <- dplyr::lead(t1)
      den <- t1 + t2
      ok <- is.finite(den) & den > 0 & t1 >= 0 & t2 >= 0
      if (!is.null(tk_max)) {
        # rk is not computed if either interval exceeds the threshold
        ok <- ok & t1 <= tk_max & t2 <= tk_max
      }
      rk <- t1 / den
      rk[!ok] <- NA_real_
      rk
    }
  )

  dplyr::ungroup(data)
}
