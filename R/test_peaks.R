#' Paired on vs off Wilcoxon test for each rhythmic category, repeated over multiple normalizations
#'
#' For each rhythmic category (e.g. \code{"1_1"}, \code{"2_1"}, \code{"3_1"},
#' ...) compares, with a paired Wilcoxon test, the normalized count of its
#' "on" bin against its "off" bin (columns \code{<category>_on} and
#' \code{<category>_off} in \code{counts}, typically produced by
#' \code{\link{count_rk}}). The test is repeated once for each normalization in
#' \code{norm_runs} (typically produced by \code{\link{simulate_null_runs}}),
#' to assess the stability of the result with respect to the Monte Carlo
#' simulation used for normalization.
#'
#' If \code{type_col} is supplied, the tests are run \strong{separately for
#' each combination of its columns}: the normalization constants are matched
#' to the rows by the values of the \code{type_col} column(s) (see
#' \code{\link{normalize_rk_counts}}) and the Wilcoxon test is then
#' computed within each combination. If \code{type_col} is \code{NULL}, a
#' single test per category is run on all the rows of \code{counts}. The
#' grouping of the tests can be changed (or extended to other columns) with
#' \code{by}.
#'
#' Normalization is performed once per run on the whole dataset; the rows
#' are split into groups only afterwards, for the tests.
#'
#' @param counts Wide data frame with one row per unit (e.g. produced by
#'   \code{\link{count_rk}}) and \code{<category>_on} /
#'   \code{<category>_off} columns.
#' @param norm_runs List of normalization data frames (the \code{type_col}
#'   column(s), if any, plus \code{bin} and \code{normalization_constant}),
#'   typically produced by \code{\link{simulate_null_runs}}.
#' @param type_col Optional column(s) of \code{counts} identifying the type
#'   (bare, string or \code{c(...)}), e.g. \code{sex}: the \code{type_col}
#'   used in \code{\link{simulate_null_runs}}. Leave \code{NULL} if the
#'   normalizations do not distinguish by type. When supplied, it also
#'   defines the groups in which the tests are run separately (unless
#'   \code{by} is given).
#' @param by Optional column(s) of \code{counts} defining the groups in
#'   which the tests are run separately (bare, string or \code{c(...)}). If
#'   \code{NULL} (default) it coincides with \code{type_col} (a single test
#'   if that is \code{NULL} too). Use it to split by other columns, e.g.
#'   \code{type_col = sex, by = c(sex, file_id)} (null distribution per sex,
#'   tests per sex x recording). Rows with \code{NA} in a \code{by} column
#'   are not tested.
#' @param peaks Categories to test (e.g. \code{c("1_1", "2_1")}). If
#'   \code{NULL} (default), all categories for which both the
#'   \code{<category>_on} and \code{<category>_off} columns exist in
#'   \code{counts}. Only columns named like a \code{\link{classify_rk}} bin
#'   (e.g. \code{"1_2_on"}, \code{"3_1_off"}) are considered.
#' @param alternative Passed to \code{wilcox.test()}: \code{"greater"}
#'   (default, on > off), \code{"less"} (on < off) or \code{"two.sided"}.
#' @param na_replace Value used to replace any non-finite values (\code{NA},
#'   \code{NaN}, \code{Inf}, \code{-Inf}) in the normalized counts before the
#'   test (default 0). They arise, in particular, when the normalization
#'   constant of a bin is 0 (bin never occurred in the simulation). It is
#'   \strong{not} used for bins with no normalization constant at all: those
#'   categories are not tested (see Details).
#' @param alpha Significance threshold used in the summary (default 0.05).
#' @param verbose If \code{TRUE} (default), prints a progress message every
#'   100 runs.
#'
#' @return A list with two elements:
#'   \item{p_values}{Data frame with one row per run x group (one row per
#'     run if there is no grouping). Columns: \code{run}, the grouping
#'     column(s) (if any) and one column \code{p_<category>} for each tested
#'     category.}
#'   \item{summary}{Data frame with one row per group x category (the
#'     grouping column(s) come first, if any): the tested hypothesis
#'     (\code{test}, e.g. \code{"1_1_on > 1_1_off"}), number of valid runs,
#'     number and proportion/percentage of significant runs
#'     (\code{p < alpha}).}
#'
#' @details
#' A category can be tested in a group only if, in the normalization of that
#' run, a constant exists for both \code{<category>_on} and
#' \code{<category>_off} for every row of the group (e.g. it does not exist
#' for \code{"1_7"} if the null distributions were computed with
#' \code{n_max = 5}, or for a type missing from \code{norm_runs}). Otherwise
#' the p-value of that run x group x category is \code{NA} (so it is not
#' counted in \code{n_run}) and a warning lists the categories (and groups)
#' that could not be tested.
#'
#' @examples
#' y <- calculate_tk(data_rktool, onset_col = onset, by = c(file_id, ind_id))
#' y <- calculate_rk(y, tk_col = tk, by = c(file_id, ind_id))
#' y$bins <- classify_rk(y$rk)
#' counts <- count_rk(y, by = c(file_id, ind_id, sex), bin_col = bins)
#'
#' # 10 null runs per sex, few simulated pairs: a quick example only (in a
#' # real analysis use at least the defaults of simulate_null_runs())
#' norm_runs <- simulate_null_runs(y, tk_col = tk, type_col = sex,
#'                                 distribution = "lognormal",
#'                                 n_sim = 1e4, runs = 10, seed = 1)
#'
#' # separate tests for each sex (type_col implies the split)
#' res <- test_peaks(counts, norm_runs, type_col = sex)
#' res$summary
#'
#' # a single null distribution and a single test per category on all rows
#' norm_all <- simulate_null_runs(y, tk_col = tk, distribution = "lognormal",
#'                                n_sim = 1e4, runs = 10, seed = 1)
#' res_all <- test_peaks(counts, norm_all)
#'
#' # null per sex, tests split by recording as well
#' res_file <- test_peaks(counts, norm_runs, type_col = sex,
#'                        by = c(sex, file_id))
#'
#' @export
test_peaks <- function(counts, norm_runs, type_col = NULL, by = NULL,
                       peaks = NULL, alternative = "greater",
                       na_replace = 0, alpha = 0.05, verbose = TRUE) {
  
  stopifnot(is.data.frame(counts), is.list(norm_runs), length(norm_runs) > 0)
  alternative <- match.arg(alternative, c("greater", "less", "two.sided"))

  # type_col: used to match the normalization constants (one or more columns)
  type_name <- .rk_select(rlang::enquo(type_col), counts, "type_col")

  # by: grouping for the tests; defaults to type_col
  by_quo  <- rlang::enquo(by)
  by_cols <- if (rlang::quo_is_null(by_quo)) type_name else .rk_select(by_quo, counts, "by")
  
  # only columns named like a classify_rk() bin
  on_cols  <- grep("^[0-9]+_[0-9]+_on$", names(counts), value = TRUE)
  off_cols <- grep("^[0-9]+_[0-9]+_off$", names(counts), value = TRUE)
  auto_peaks <- intersect(sub("_on$", "", on_cols), sub("_off$", "", off_cols))
  
  if (is.null(peaks)) {
    peaks <- auto_peaks
  } else {
    missing_peaks <- setdiff(peaks, auto_peaks)
    if (length(missing_peaks) > 0) {
      stop("Missing on/off columns for: ", paste(missing_peaks, collapse = ", "))
    }
  }
  if (length(peaks) == 0) {
    stop("No <category>_on / <category>_off column pair found in counts.")
  }
  
  # groups of rows on which each test is run
  if (length(by_cols) > 0) {
    grp_id <- interaction(as.list(counts[by_cols]), drop = TRUE, lex.order = TRUE)
    groups <- split(seq_len(nrow(counts)), grp_id)
    keys <- as.data.frame(counts[vapply(groups, `[`, integer(1), 1), by_cols, drop = FALSE])
  } else {
    groups <- list(seq_len(nrow(counts)))
    keys <- NULL
  }
  G <- length(groups)
  runs <- length(norm_runs)
  
  p_mat <- matrix(NA_real_, nrow = runs * G, ncol = length(peaks))
  bin_cols <- c(paste0(peaks, "_on"), paste0(peaks, "_off"))
  # group x peak combinations that could not be tested in at least one run
  excluded <- matrix(FALSE, nrow = G, ncol = length(peaks))

  for (r in seq_len(runs)) {
    # normalization constant of each row x bin (NA where norm_df has none);
    # an empty all_of() selection means "no type" in normalize_rk_counts()
    const <- normalize_rk_counts(
      counts, norm_runs[[r]],
      type_col = tidyselect::all_of(type_name),
      bin_cols = bin_cols, suffix = NULL,
      operation = function(count, norm) norm
    )
    # Normalized counts. Replace any non-finite value (NA, NaN, Inf, -Inf):
    # Inf/-Inf arise when the normalization constant is 0. Missing constants
    # are handled below (the peak is not tested).
    normalized <- lapply(stats::setNames(bin_cols, bin_cols), function(b) {
      x <- counts[[b]] / const[[b]]
      x[!is.finite(x)] <- na_replace
      x
    })

    for (g in seq_len(G)) {
      idx <- groups[[g]]
      for (k in seq_along(peaks)) {
        on  <- paste0(peaks[k], "_on")
        off <- paste0(peaks[k], "_off")
        if (anyNA(const[[on]][idx]) || anyNA(const[[off]][idx])) {
          excluded[g, k] <- TRUE
          next
        }
        test <- try(
          stats::wilcox.test(normalized[[on]][idx], normalized[[off]][idx],
                             paired = TRUE, alternative = alternative),
          silent = TRUE
        )
        p_mat[(r - 1) * G + g, k] <- if (inherits(test, "try-error")) NA_real_ else test$p.value
      }
    }
    if (verbose && r %% 100 == 0) message("Completed test ", r, "/", runs)
  }

  if (any(excluded)) {
    w <- which(excluded, arr.ind = TRUE)
    lab <- peaks[w[, 2]]
    if (!is.null(keys)) {
      grp_lab <- do.call(paste, c(lapply(keys, as.character), sep = ", "))
      lab <- paste0(lab, " [", grp_lab[w[, 1]], "]")
    }
    warning("No normalization constant in norm_runs for: ", paste(lab, collapse = "; "),
            ". These categories were not tested (p-value NA) in the affected run(s).")
  }
  
  colnames(p_mat) <- paste0("p_", peaks)
  p_values <- data.frame(run = rep(seq_len(runs), each = G))
  if (!is.null(keys)) {
    p_values <- cbind(p_values, keys[rep(seq_len(G), times = runs), , drop = FALSE])
  }
  p_values <- cbind(p_values, as.data.frame(p_mat))
  rownames(p_values) <- NULL
  
  # label of the tested hypothesis, consistent with 'alternative'
  op <- switch(alternative, greater = " > ", less = " < ", two.sided = " vs ")

  summary_df <- do.call(rbind, lapply(seq_len(G), function(g) {
    rows <- seq(g, by = G, length.out = runs)
    do.call(rbind, lapply(peaks, function(peak) {
      pcol  <- p_values[[paste0("p_", peak)]][rows]
      n_run <- sum(!is.na(pcol))
      sig   <- sum(pcol < alpha, na.rm = TRUE)
      out <- data.frame(
        test  = paste0(peak, "_on", op, peak, "_off"),
        n_run = n_run,
        sig   = sig,
        perc  = if (n_run > 0) sig / n_run * 100 else NA_real_,
        prop  = if (n_run > 0) sig / n_run else NA_real_
      )
      if (!is.null(keys)) out <- cbind(keys[g, , drop = FALSE], out)
      out
    }))
  }))
  rownames(summary_df) <- NULL
  
  list(p_values = p_values, summary = summary_df)
}
