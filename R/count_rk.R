#' Count occurrences of each rhythmic bin per group, in wide or long format
#'
#' Counts how many times each rhythmic bin (e.g. the \code{bins} column
#' produced by \code{\link{classify_rk}}) occurs within each group (e.g.
#' \code{file_id}, \code{ind_id}). The result can be returned in two
#' formats:
#' \itemize{
#'   \item \code{"wide"} (default): one row per combination of \code{by} and
#'     one column per bin, containing the corresponding count.
#'   \item \code{"long"}: one row per combination of \code{by} and bin, with
#'     the bin label in a column named as \code{bin_col} and the count in a
#'     column \code{n}.
#' }
#' Both formats contain the same information: group x bin combinations
#' absent from the data are included with the value \code{fill} (default 0)
#' also in the long format. Moreover, for every category with at least one
#' observed bin, both \code{<category>_on} and \code{<category>_off} are
#' always present: if only one of the two occurs in the data, the other is
#' added with the value \code{fill} (so that, e.g., \code{\link{test_peaks}}
#' can test it).
#'
#' @param data A data frame containing both the grouping columns and the
#'   bin column.
#' @param by Grouping column(s) (bare, string or \code{c(...)}), e.g.
#'   \code{c(file_id, ind_id)}: one count per combination of their values.
#' @param bin_col Column with the bins to count (bare or string), e.g.
#'   \code{bins}. In the wide format its values become column names; in the
#'   long format it is the name of the column with the bin labels. If it is
#'   a factor, the order of its levels is kept (as column order in the wide
#'   format, as factor levels in the long format); otherwise bins are in
#'   alphabetical order.
#' @param format Output format: \code{"wide"} (default) or \code{"long"}.
#' @param keep_cols Optional column(s) to keep in the result (bare, string
#'   or \code{c(...)}), e.g. \code{sex}. They must be constant within each
#'   \code{by} group: otherwise the result contains more than one row for
#'   that group, one for each of their values.
#' @param drop_na_bin Logical, if \code{TRUE} (default) rows with
#'   \code{bin_col == NA} are excluded from the count.
#' @param fill Value used to fill group x bin combinations absent from the
#'   original data (default 0).
#'
#' @return A tibble. With \code{format = "wide"}: one row per combination of
#'   \code{by}, the \code{keep_cols} columns, and one column for each bin
#'   containing the count. With \code{format = "long"}: one row per
#'   combination of \code{by} and bin, the \code{keep_cols} columns, a column
#'   named as \code{bin_col} (factor, with bins in the same order as the wide
#'   columns) and a column \code{n} with the count.
#'
#' @examples
#' y <- calculate_tk(data_rktool, onset_col = onset, by = c(file_id, ind_id))
#' y <- calculate_rk(y, tk_col = tk, by = c(file_id, ind_id))
#' y$bins <- classify_rk(y$rk)
#'
#' # wide format (default): one row per recording x individual
#' count_rk(y, by = c(file_id, ind_id), bin_col = bins, keep_cols = sex)
#'
#' # long format (strings work too)
#' count_rk(y, by = c("file_id", "ind_id"), bin_col = "bins", format = "long")
#'
#' @export
count_rk <- function(data, by, bin_col, format = c("wide", "long"),
                     keep_cols = NULL, drop_na_bin = TRUE, fill = 0) {
  format <- match.arg(format)
  stopifnot(is.data.frame(data))

  by        <- .rk_select(rlang::enquo(by), data, "by")
  bin_col   <- .rk_select(rlang::enquo(bin_col), data, "bin_col", one = TRUE)
  keep_cols <- .rk_select(rlang::enquo(keep_cols), data, "keep_cols")

  level_order <- if (is.factor(data[[bin_col]])) levels(data[[bin_col]]) else NULL

  if (drop_na_bin) {
    data <- data[!is.na(data[[bin_col]]), , drop = FALSE]
  }

  counts <- dplyr::summarise(
    dplyr::group_by(data, dplyr::across(dplyr::all_of(c(by, bin_col)))),
    n = dplyr::n(),
    dplyr::across(dplyr::all_of(keep_cols), dplyr::first),
    .groups = "drop"
  )

  wide <- tidyr::pivot_wider(
    counts,
    names_from = dplyr::all_of(bin_col),
    values_from = "n",
    values_fill = fill
  )

  # every category with at least one observed bin gets both "<category>_off"
  # and "<category>_on": the missing one is added, filled with 'fill' (off
  # before on)
  bin_names <- setdiff(names(wide), c(by, keep_cols))
  new_order <- character(0)
  for (b in bin_names) {
    partner <- if (grepl("_on$", b)) sub("_on$", "_off", b) else if (grepl("_off$", b)) sub("_off$", "_on", b)
    add <- !is.null(partner) && !partner %in% bin_names
    if (add) wide[[partner]] <- rep(fill, nrow(wide))
    new_order <- c(new_order,
                   if (add && grepl("_off$", partner)) c(partner, b)
                   else if (add) c(b, partner)
                   else b)
  }
  wide <- wide[c(by, keep_cols, new_order)]

  if (!is.null(level_order)) {
    wide <- dplyr::select(
      wide,
      dplyr::all_of(c(by, keep_cols)), dplyr::any_of(level_order), dplyr::everything()
    )
  }

  if (format == "wide") {
    return(wide)
  }

  # long format: same information as the wide one, one row per group x bin
  bin_names <- setdiff(names(wide), c(by, keep_cols))

  long <- tidyr::pivot_longer(
    wide,
    cols = dplyr::all_of(bin_names),
    names_to = bin_col,
    values_to = "n"
  )
  # keep bins in the same order as the wide columns
  long[[bin_col]] <- factor(long[[bin_col]], levels = bin_names)

  long
}
