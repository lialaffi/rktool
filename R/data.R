#' Example dataset of event onsets
#'
#' Onset times of events annotated in 183 recordings, for 8 individuals
#' (4 females, 4 males). It is a real dataset meant to run the examples of
#' the package: starting from the onsets, the whole workflow can be
#' reproduced (\code{\link{calculate_tk}}
#' -> \code{\link{calculate_rk}} -> \code{\link{classify_rk}} ->
#' \code{\link{count_rk}} -> \code{\link{compute_null_distributions}} ->
#' \code{\link{normalize_rk_counts}} -> \code{\link{test_peaks}}).
#'
#' Rows are not sorted: \code{\link{calculate_tk}} sorts them by onset
#' within each group.
#'
#' Recordings and individuals are combined in 326 file x individual
#' sequences, some of them very short (a few events, or a single one, which
#' gives no interval).
#'
#' @format A data frame with 11983 rows and 4 columns:
#' \describe{
#'   \item{file_id}{Recording identifier (\code{"file_1"} ... \code{"file_183"}).}
#'   \item{ind_id}{Individual identifier (\code{"ind1"} ... \code{"ind8"}).}
#'   \item{sex}{Sex of the individual: \code{"female"} or \code{"male"}.}
#'   \item{onset}{Onset time of the event, in seconds from the start of the
#'     recording.}
#' }
#'
#' @examples
#' data(data_rktool)
#' head(data_rktool)
#' y <- calculate_tk(data_rktool, onset_col = onset, by = c(file_id, ind_id))
#' y <- calculate_rk(y, tk_col = tk, by = c(file_id, ind_id))
#' summary(y$rk)
"data_rktool"
