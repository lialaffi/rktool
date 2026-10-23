#' Resolve a tidyselect argument to column names
#'
#' Internal helper shared by the exported functions: evaluates a quosure
#' (bare names, strings, \code{c(...)} or any tidyselect expression) against
#' \code{data} and returns the selected column names. A \code{NULL} quosure
#' gives \code{character(0)}.
#'
#' @param quo A quosure, as returned by \code{rlang::enquo()}.
#' @param data A data frame.
#' @param arg Name of the argument, used in error messages.
#' @param one Logical, if \code{TRUE} exactly one column must be selected.
#' @param numeric Logical, if \code{TRUE} the selected column(s) must be
#'   numeric.
#' @return Character vector of column names.
#' @keywords internal
.rk_select <- function(quo, data, arg, one = FALSE, numeric = FALSE) {
  cols <- if (rlang::quo_is_null(quo)) character(0) else
    names(tidyselect::eval_select(quo, data))
  if (one && length(cols) != 1) stop("'", arg, "' must select exactly one column.")
  if (numeric) {
    for (cc in cols) {
      if (!is.numeric(data[[cc]])) stop("Column '", cc, "' must be numeric.")
    }
  }
  cols
}


#' Validate or generate a set of rhythmic peak labels
#'
#' Internal helper for \code{\link{classify_rk}}, \code{\link{find_rk_peaks}}
#' and \code{\link{plot_rk_density}}: if \code{peaks} is \code{NULL}
#' generates the canonical sequence \code{1_n_max ... 1_1 ... n_max_1};
#' otherwise checks that \code{peaks} is a non-empty character vector
#' without duplicates.
#'
#' @param peaks \code{NULL} or a character vector of peak labels.
#' @param n_max Integer >= 1 (default 5), used only when \code{peaks} is
#'   \code{NULL}.
#' @return Character vector of peak labels.
#' @keywords internal
.rk_peaks <- function(peaks = NULL, n_max = 5) {
  if (is.null(peaks)) {
    if (!is.numeric(n_max) || length(n_max) != 1 || is.na(n_max) ||
        n_max < 1 || n_max != as.integer(n_max)) {
      stop("n_max must be an integer >= 1")
    }
    if (n_max == 1) return("1_1")
    return(c(paste0("1_", n_max:2), "1_1", paste0(2:n_max, "_1")))
  }
  if (!is.character(peaks) || length(peaks) == 0) {
    stop("'peaks' must be a non-empty character vector.")
  }
  if (anyDuplicated(peaks)) {
    stop("Duplicated peaks in 'peaks': ", paste(unique(peaks[duplicated(peaks)]), collapse = ", "))
  }
  peaks
}


#' Key identifying the combination of values of one or more columns
#'
#' Internal helper: pastes the values (as character) of the columns
#' \code{cols} of \code{df} into a single string per row, so that
#' combinations of several columns can be matched with \code{match()},
#' \code{==} or named lookups. Used by \code{\link{compute_null_distributions}}
#' and \code{\link{normalize_rk_counts}}.
#'
#' @param df A data frame.
#' @param cols Character vector of column names (at least one).
#' @return Character vector with one key per row of \code{df}.
#' @keywords internal
.rk_make_key <- function(df, cols) {
  if (length(cols) == 1) return(as.character(df[[cols]]))
  do.call(paste, c(lapply(cols, function(cc) as.character(df[[cc]])), sep = "\r"))
}
