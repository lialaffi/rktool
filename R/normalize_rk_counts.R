#' Width of the rk region of rhythmic bins
#'
#' Internal helper for \code{\link{normalize_rk_counts}}: for bin labels as
#' produced by \code{\link{classify_rk}} (e.g. \code{"1_3_on"},
#' \code{"1_1_off"}, \code{"4_1_on"}) returns the total width, in rk units,
#' of the region covered by the bin. The breakpoints are the same used by
#' \code{classify_rk()} (via \code{.parse_rk_peak()} and
#' \code{.rk_category_breaks()}). Being \code{b1 < ... < b5} the breakpoints
#' of a category, the "on" bin covers \code{(b2, b4]} (width \code{b4 - b2})
#' and the "off" bin covers \code{(b1, b2]} and \code{(b4, b5]} (width
#' \code{(b2 - b1) + (b5 - b4)}).
#'
#' @param bins Character vector of bin labels.
#' @return Numeric vector of widths, same length as \code{bins}; \code{NA}
#'   for labels that are \code{NA} or do not have the expected format.
#' @keywords internal
.rk_bin_widths <- function(bins) {
  width_one <- function(label) {
    if (is.na(label)) return(NA_real_)
    m <- regmatches(label, regexec("^([0-9]+_[0-9]+)_(on|off)$", label))[[1]]
    if (length(m) != 3) return(NA_real_)
    p <- tryCatch(.parse_rk_peak(m[2]), error = function(e) NULL)
    if (is.null(p)) return(NA_real_)
    b <- .rk_category_breaks(p$n, p$side)$breaks
    if (m[3] == "on") b[4] - b[2] else (b[2] - b[1]) + (b[5] - b[4])
  }

  vapply(as.character(bins), width_one, numeric(1), USE.NAMES = FALSE)
}


#' Normalize rhythmic bin counts
#'
#' Divides the counts of each rhythmic bin (e.g. produced by
#' \code{\link{count_rk}}) by a normalization constant (by bin, and
#' optionally by type). Two kinds of constants can be used, selected with
#' \code{method}:
#' \itemize{
#'   \item \code{"null"} (default): the proportions produced by
#'     \code{\link{compute_null_distributions}}, i.e. the proportion of each
#'     bin expected under a null distribution of the intervals.
#'   \item \code{"width"}: simply the width, in rk units, of the region
#'     covered by each bin (computed from the same breakpoints used by
#'     \code{\link{classify_rk}}), so that counts become comparable across
#'     bins of different size. No \code{norm_df} is needed.
#' }
#' Both the wide format (one column per bin) and the long format (one
#' row per group x bin) of \code{\link{count_rk}} are supported.
#'
#' @param counts Data frame of counts, as produced by
#'   \code{\link{count_rk}}: wide (one column per rhythmic bin) or long (a
#'   column with the bin labels and a column with the counts). If a type is
#'   used for the normalization, it must also contain the column(s)
#'   identifying it (\code{type_col}).
#' @param norm_df Long data frame with the type column(s) (if any), and the
#'   columns \code{bin} and \code{normalization_constant} (produced by
#'   \code{\link{compute_null_distributions}}). Required when
#'   \code{method = "null"}; ignored when \code{method = "width"}.
#' @param type_col Optional column(s) of \code{counts} identifying the type
#'   (bare, string or \code{c(...)}), e.g. \code{sex}: matched against the
#'   columns of \code{norm_df} with the same names (the \code{type_col} used
#'   in \code{\link{compute_null_distributions}}). With several columns, the
#'   constant is chosen by the combination of their values. If \code{NULL}
#'   (default) no matching by type is performed and \code{norm_df} must
#'   have a single row per bin, applied to all rows. Ignored when
#'   \code{method = "width"}.
#' @param bin_cols Wide format only. Names of the bin columns to normalize
#'   (character vector). If \code{NULL} (default), all columns of
#'   \code{counts} whose name appears among the values of \code{norm_df$bin}
#'   are used automatically (with \code{method = "width"}: all columns named
#'   like a \code{classify_rk()} bin, e.g. \code{"1_2_on"}, \code{"1_1_off"}).
#' @param operation \code{"divide"} (default, \code{count / norm}), or a
#'   function \code{function(count, norm)} you supply for another operation.
#' @param suffix Suffix added to the names of the new normalized columns
#'   (default \code{"_norm"}). If \code{NULL}, the original columns are
#'   overwritten with the normalized values.
#' @param format Format of \code{counts}: \code{"wide"} (default) or
#'   \code{"long"}.
#' @param bin_col Long format only (required). Column of \code{counts} with
#'   the bin labels (bare or string), e.g. \code{bins}: matched against the
#'   \code{bin} column of \code{norm_df}.
#' @param count_col Long format only. Column of \code{counts} with the
#'   counts to normalize (bare or string; default \code{"n"}, as produced by
#'   \code{\link{count_rk}} with \code{format = "long"}).
#' @param keep_counts Logical, if \code{TRUE} (default) the original counts
#'   are kept next to the normalized ones. If \code{FALSE}, the original
#'   count columns (the bin columns in the wide format, \code{count_col} in
#'   the long format) are removed, so the result contains only the
#'   identifying columns and the normalized values. Has no effect if
#'   \code{suffix = NULL}, since the originals are already overwritten.
#' @param add_norm_col Long format only. Logical, if \code{TRUE} (default) a
#'   column named \code{norm} is added with the normalization constant
#'   applied to each row (the null proportion or the bin width, depending on
#'   \code{method}; \code{NA} where no constant is available). Ignored in
#'   the wide format.
#' @param method How to obtain the normalization constants: \code{"null"}
#'   (default), using \code{norm_df}, or \code{"width"}, using the width of
#'   the rk region of each bin.
#'
#' @return The \code{counts} data frame with the normalized values added (or
#'   overwritten if \code{suffix = NULL}). In the wide format, one new
#'   column per normalized bin (\code{<bin><suffix>}). In the long format, a
#'   new column \code{<count_col><suffix>} (e.g. \code{n_norm}) and, if
#'   \code{add_norm_col = TRUE}, a column \code{norm} with the constant used.
#'
#' @details
#' If no constant is available for a bin (with \code{method = "null"}: no
#' row in \code{norm_df} for that type + bin combination, e.g. a category
#' beyond the \code{n_max} used in \code{\link{compute_null_distributions}},
#' or a combination of \code{type_col} values not in \code{norm_df}; with
#' \code{method = "width"}: a label that is not a valid
#' \code{classify_rk()} bin), the resulting normalized value is \code{NA}.
#' In the wide format with \code{bin_cols = NULL}, bin columns with no row
#' at all in \code{norm_df} are not selected, so they are simply left
#' unnormalized (no \code{<bin><suffix>} column is added for them).
#' A bin that never occurred in the simulation has a constant of 0, which
#' produces \code{Inf}/\code{NaN}: check \code{norm_df} for bins with a
#' constant of 0 or close to 0.
#'
#' With \code{method = "width"}, the width of a bin is the total length of
#' the rk interval(s) it covers (see \code{\link{classify_rk}}), so the
#' result is a count per unit of rk (a density-like measure).
#'
#' @examples
#' y <- calculate_tk(rk_example, onset_col = onset, by = c(file_id, ind_id))
#' y <- calculate_rk(y, tk_col = tk, by = c(file_id, ind_id))
#' y$bins <- classify_rk(y$rk)
#'
#' # null distribution per sex (n_sim kept small for a quick example)
#' norm_df <- compute_null_distributions(y, tk_col = tk, type_col = sex,
#'                                       distribution = "lognormal",
#'                                       n_sim = 1e4, seed = 1)
#'
#' # wide counts
#' b <- count_rk(y, by = c(file_id, ind_id, sex), bin_col = bins)
#' b_norm <- normalize_rk_counts(b, norm_df, type_col = sex)
#'
#' # normalize by the width of the bins instead of the null distribution
#' b_width <- normalize_rk_counts(b, method = "width")
#'
#' # long counts
#' b_long <- count_rk(y, by = c(file_id, ind_id, sex), bin_col = bins,
#'                    format = "long")
#' b_long_norm <- normalize_rk_counts(b_long, norm_df, type_col = sex,
#'                                    format = "long", bin_col = bins)
#' head(b_long_norm)
#'
#' @export
normalize_rk_counts <- function(counts, norm_df = NULL, type_col = NULL, bin_cols = NULL,
                                operation = "divide", suffix = "_norm",
                                format = c("wide", "long"), bin_col = NULL,
                                count_col = "n", keep_counts = TRUE, add_norm_col = TRUE,
                                method = c("null", "width")) {
  
  format <- match.arg(format)
  method <- match.arg(method)
  
  stopifnot(is.data.frame(counts))
  if (method == "null") {
    if (is.null(norm_df)) {
      stop("With method = \"null\", 'norm_df' must be supplied (see compute_null_distributions()).")
    }
    stopifnot(is.data.frame(norm_df))
    if (!all(c("bin", "normalization_constant") %in% names(norm_df))) {
      stop("norm_df must have columns 'bin' and 'normalization_constant' (see compute_null_distributions()).")
    }
  }
  
  # bin widths do not depend on the type: type_col is ignored with method = "width"
  # (an empty selection also means "no type": do not match by type)
  type_col <- if (method == "null") .rk_select(rlang::enquo(type_col), counts, "type_col") else character(0)
  use_type <- length(type_col) > 0
  
  # columns of norm_df matched against the type_col columns of counts
  norm_type_cols <- character(0)
  if (use_type) {
    if (all(type_col %in% names(norm_df))) {
      norm_type_cols <- type_col
    } else {
      stop("norm_df has no column(s) named ",
           paste0("'", setdiff(type_col, names(norm_df)), "'", collapse = ", "),
           " (the columns used as type_col in compute_null_distributions()).")
    }
  }
  
  if (is.function(operation)) {
    op_fun <- operation
  } else {
    operation <- match.arg(operation, "divide")
    op_fun <- function(count, norm) count / norm
  }
  
  # ---------------------------------------------------------------- long
  if (format == "long") {
    if (!is.null(bin_cols)) {
      warning("'bin_cols' is ignored when format = \"long\".")
    }
    bin_quo <- rlang::enquo(bin_col)
    if (rlang::quo_is_null(bin_quo)) {
      stop("With format = \"long\", 'bin_col' (the column with the bin labels) must be specified.")
    }
    bin_col   <- .rk_select(bin_quo, counts, "bin_col", one = TRUE)
    count_col <- .rk_select(rlang::enquo(count_col), counts, "count_col", one = TRUE, numeric = TRUE)
    
    if (method == "width") {
      bins_present <- unique(as.character(counts[[bin_col]]))
      bins_present <- bins_present[!is.na(bins_present)]
      norm_df <- data.frame(
        bin = bins_present,
        normalization_constant = .rk_bin_widths(bins_present),
        stringsAsFactors = FALSE
      )
    }
    
    # keys: combination of type_col values (if any) + bin label
    counts_key_cols <- if (use_type) c(type_col, bin_col) else bin_col
    norm_key_cols   <- if (use_type) c(norm_type_cols, "bin") else "bin"
    key_counts <- .rk_make_key(counts, counts_key_cols)
    key_norm   <- .rk_make_key(norm_df, norm_key_cols)
    if (anyDuplicated(key_norm)) {
      stop("norm_df has more than one row for the same ",
           if (use_type) "type + bin" else "bin",
           if (!use_type) ": specify type_col to match by type.")
    }
    
    lookup <- stats::setNames(norm_df$normalization_constant, key_norm)
    norm_values <- unname(lookup[key_counts])
    
    new_col <- if (is.null(suffix)) count_col else paste0(count_col, suffix)
    if (add_norm_col) {
      if ("norm" %in% names(counts)) {
        warning("Column 'norm' already exists and will be overwritten.")
      }
      counts[["norm"]] <- norm_values
    }
    counts[[new_col]] <- op_fun(counts[[count_col]], norm_values)
    if (!keep_counts && new_col != count_col) {
      counts[[count_col]] <- NULL
    }
    return(counts)
  }
  
  # ---------------------------------------------------------------- wide
  if (method == "width") {
    candidates <- if (is.null(bin_cols)) {
      grep("^[0-9]+_[0-9]+_(on|off)$", names(counts), value = TRUE)
    } else {
      bin_cols
    }
    if (length(candidates) == 0) {
      stop("No bin columns found in counts (expected names like '1_2_on', '1_1_off').")
    }
    norm_df <- data.frame(
      bin = candidates,
      normalization_constant = .rk_bin_widths(candidates),
      stringsAsFactors = FALSE
    )
  }
  
  if (is.null(bin_cols)) {
    bin_cols <- intersect(names(counts), unique(norm_df$bin))
  }
  missing_bins <- setdiff(bin_cols, names(counts))
  if (length(missing_bins) > 0) {
    stop("Bin columns not found in counts: ", paste(missing_bins, collapse = ", "))
  }
  
  # key of the type_col combination of each row of counts (computed once)
  counts_type_key <- if (use_type) .rk_make_key(counts, type_col) else NULL
  
  for (bin in bin_cols) {
    norm_bin <- norm_df[norm_df$bin == bin, , drop = FALSE]

    if (use_type) {
      key_bin <- .rk_make_key(norm_bin, norm_type_cols)
      if (anyDuplicated(key_bin)) {
        stop("norm_df has more than one row for the same type + bin ('", bin, "').")
      }
      norm_values <- norm_bin$normalization_constant[match(counts_type_key, key_bin)]
    } else {
      if (nrow(norm_bin) > 1) {
        stop("norm_df has more than one row for bin '", bin, "': specify type_col to match by type.")
      }
      # no row for this bin in norm_df: no constant available, NA
      norm_const <- if (nrow(norm_bin) == 1) norm_bin$normalization_constant else NA_real_
      norm_values <- rep(norm_const, nrow(counts))
    }
    
    new_col <- if (is.null(suffix)) bin else paste0(bin, suffix)
    counts[[new_col]] <- op_fun(counts[[bin]], norm_values)
  }
  
  if (!keep_counts && !is.null(suffix)) {
    counts <- counts[, setdiff(names(counts), bin_cols), drop = FALSE]
  }
  
  counts
}