#' Calculate intervals (tk) from onsets
#'
#' Computes, within each temporal sequence (e.g. recording x individual), the
#' intervals between consecutive onsets: rows are sorted by onset in
#' ascending order and \code{tk = onset - previous onset}. The result can
#' then be passed to \code{\link{calculate_rk}}.
#'
#' @param data A data frame.
#' @param onset_col Column with the onset times (bare or string).
#' @param by Optional grouping column(s) (bare, string or \code{c(...)}),
#'   e.g. \code{c(file_id, ind_id)}. If \code{NULL} (default) the whole
#'   dataset is treated as a single group.
#' @param tk_name Name of the new interval column (default \code{"tk"}).
#'
#' @return An ungrouped data frame (tibble) with all the original columns
#'   plus the \code{tk} column, sorted by \code{by} and onset. The first
#'   \code{tk} of each group (no previous onset) is \code{NA}.
#'
#' @details
#' Rows with \code{NA} onsets are not removed: they are sorted last within
#' their group and produce \code{NA} intervals. Remove them beforehand if
#' needed.
#'
#' @examples
#' # intervals within each recording x individual
#' y <- calculate_tk(rk_example, onset_col = onset, by = c(file_id, ind_id))
#' head(y)
#'
#' # then the rhythmic ratios
#' y <- calculate_rk(y, tk_col = tk, by = c(file_id, ind_id))
#'
#' @export
calculate_tk <- function(data, onset_col, by = NULL, tk_name = "tk") {
  stopifnot(is.data.frame(data))

  onset_col <- .rk_select(rlang::enquo(onset_col), data, "onset_col", one = TRUE, numeric = TRUE)
  by_cols   <- .rk_select(rlang::enquo(by), data, "by")

  if (tk_name %in% names(data) && tk_name != onset_col) {
    warning("Column '", tk_name, "' already exists and will be overwritten.")
  }

  data <- dplyr::arrange(data, dplyr::across(dplyr::all_of(c(by_cols, onset_col))))
  data <- dplyr::group_by(data, dplyr::across(dplyr::all_of(by_cols)))

  # tk: onset minus previous onset (NA for the first row of a group)
  data <- dplyr::mutate(
    data,
    !!tk_name := .data[[onset_col]] - dplyr::lag(.data[[onset_col]])
  )

  dplyr::ungroup(data)
}
