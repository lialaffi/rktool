#' rktool: rhythmic categories from generic event sequences
#'
#' Compute rhythmic ratios (rk), classify them into the "on" / "off" bins of
#' rhythmic categories, count and normalize the bins against a null
#' distribution, and test "on" versus "off" bins for each category.
#'
#' Typical workflow: \code{\link{calculate_tk}} -> \code{\link{calculate_rk}} ->
#' \code{\link{classify_rk}} -> \code{\link{count_rk}} ->
#' \code{\link{compute_null_distributions}} ->
#' \code{\link{normalize_rk_counts}} -> \code{\link{test_peaks}}.
#' The example dataset \code{\link{data_rktool}} can be used to run it.
#'
#' @section Terminology:
#' \describe{
#'   \item{category (or peak)}{A rhythmic category, labelled \code{"1_n"},
#'     \code{"1_1"} or \code{"n_1"} (e.g. \code{"1_2"}, ratio 1:2), centered
#'     on \code{rk = 1/(n+1)}, 0.5 and \code{1 - 1/(n+1)} respectively.}
#'   \item{bin}{The "on" (central) or "off" (outer) part of a category,
#'     labelled e.g. \code{"1_2_on"} / \code{"1_2_off"}.}
#' }
#'
#' Column arguments accept a bare name or a string; those marked
#' \code{c(...)} also accept several columns (any tidyselect expression).
#'
#' @examples
#' # load the example dataset
#' data(data_rktool)
#' head(data_rktool)
#'
#' y <- calculate_tk(data_rktool, onset_col = onset, by = c(file_id, ind_id))
#' y <- calculate_rk(y, tk_col = tk, by = c(file_id, ind_id))
#' y$bins <- classify_rk(y$rk)
#' b <- count_rk(y, by = c(file_id, ind_id, sex), bin_col = bins)
#' head(b)
#'
#' @keywords internal
#' @importFrom rlang .data :=
"_PACKAGE"

