#' Compute the null distribution of rhythmic bins from observed tk
#'
#' Estimates the parameters of a null distribution from the observed
#' intervals \code{tk}, simulates pairs of random intervals from that
#' distribution, computes the corresponding rhythmic ratios
#' \code{rk = tk / (tk + tk+1)}, classifies them with
#' \code{\link{classify_rk}} and returns the proportion of simulated pairs
#' falling in each rhythmic bin. These proportions are the normalization
#' constants used by \code{\link{normalize_rk_counts}} to compare observed
#' data against what would be expected by chance. The approach follows
#' Jadoul et al. (2025).
#'
#' If \code{type_col} is supplied, the procedure is repeated independently
#' for each combination of the values of the selected column(s) (e.g.
#' \code{type_col = c(sex, ind_id)} gives one null distribution for each
#' sex x individual combination present in the data).
#'
#' @param data A data frame containing the \code{tk_col} (and, optionally,
#'   \code{type_col}) columns.
#' @param tk_col Column with the observed intervals (bare or string).
#' @param type_col Optional column(s) identifying the groups for which to
#'   estimate a separate null distribution (bare, string or \code{c(...)}),
#'   e.g. \code{sex}. One null distribution is estimated for each
#'   combination of values present in the data. If \code{NULL} (default),
#'   a single null distribution is estimated on all \code{tk_col} values
#'   (the result then has no grouping columns).
#' @param distribution Null distribution of \code{tk}, one of
#'   \code{"uniform"} (parameters \code{min}, \code{max}),
#'   \code{"exponential"} (\code{rate = 1 / mean(tk)}),
#'   \code{"halfnormal"} (\code{sigma = mean(tk) * sqrt(pi / 2)}) or
#'   \code{"lognormal"} (\code{meanlog}, \code{sdlog} of \code{log(tk)}).
#'   Default \code{"uniform"}.
#' @param n_sim Number of simulated interval pairs per group (default 1e6).
#' @param n_max Most extreme category to include (\code{1_n_max} and
#'   \code{n_max_1}), passed to \code{\link{classify_rk}} (default 5).
#' @param type_levels Optional, fixes which groups to include and in what
#'   order (default: all combinations found in \code{data}, in order of
#'   appearance). With a single \code{type_col} it can be a vector of
#'   values; in general it can be a data frame with one column for each
#'   column of \code{type_col} (same names) and one row per combination.
#'   Ignored if \code{type_col} is \code{NULL}.
#' @param seed Optional seed for simulation reproducibility. The same seed
#'   is set before simulating each group.
#' @param tk_max Optional upper bound used to filter \code{tk_col} before
#'   fitting/simulating: values with \code{tk_col > tk_max} are dropped. If
#'   \code{NULL} (default), no upper filtering is applied.
#' @param plot Logical, if \code{TRUE} draws, for each group, two density
#'   plots: (1) the simulated (null) \code{tk} together with the observed
#'   \code{tk} (dashed red line), with the estimated parameters shown as
#'   subtitle; (2) the density of the resulting simulated \code{rk}. Default
#'   \code{FALSE} (so that repeated calls, e.g. from
#'   \code{\link{simulate_null_runs}}, do not produce plots). Densities are
#'   computed on at most 1e5 simulated values.
#'
#' @return A data frame with the columns of \code{type_col} (with the same
#'   names as in \code{data}; none if \code{type_col} is \code{NULL}),
#'   followed by \code{bin} and \code{normalization_constant} (simulated
#'   proportion of each bin). The parameters estimated from the observed
#'   \code{tk} are attached as the attribute \code{"params"}, a data frame
#'   with the \code{type_col} columns followed by \code{parameter} and
#'   \code{value}, retrievable with \code{attr(x, "params")}.
#'
#' @details
#' Before fitting/simulating, rows with non-positive or missing
#' \code{tk_col} values are always dropped (and, if \code{tk_max} is
#' supplied, rows with \code{tk_col > tk_max} as well); this filtering is
#' applied to the whole dataset, before any splitting by \code{type_col}.
#' Rows with \code{NA} in any \code{type_col} column are then dropped, with
#' a warning.
#' When \code{distribution = "uniform"}, the \code{min} parameter is always
#' estimated from the dataset filtered on \code{tk_col} (the same minimum is
#' used for every group), while \code{max} is estimated per group. The
#' parameters of the other distributions are estimated per group.
#'
#' With several \code{type_col} columns the number of groups is the number
#' of combinations present in the data, and some may contain very few
#' \code{tk} values: a warning is issued for groups with fewer than 10
#' values, since the parameters (and hence the null distribution) are then
#' poorly estimated.
#'
#' @references
#' Jadoul, Y., Tufarelli, T., Coissac, C., Gamba, M., & Ravignani, A.
#' (2025). Hidden assumptions of integer ratio analyses in bioacoustics and
#' music. \emph{Annals of the New York Academy of Sciences}, 1553(1),
#' 363-378.
#'
#' @examples
#' y <- calculate_tk(rk_example, onset_col = onset, by = c(file_id, ind_id))
#'
#' # one null distribution per sex (n_sim kept small for a quick example)
#' norm_df <- compute_null_distributions(
#'   y, tk_col = tk, type_col = sex,
#'   distribution = "lognormal", n_sim = 1e4, seed = 1
#' )
#' head(norm_df)
#' attr(norm_df, "params")   # parameters estimated from the observed tk
#'
#' # without grouping: a single null distribution over all tk, with plots
#' norm_df_all <- compute_null_distributions(
#'   y, tk_col = tk, distribution = "lognormal", n_sim = 1e4, seed = 1,
#'   plot = TRUE
#' )
#'
#' @export
compute_null_distributions <- function(data, tk_col, type_col = NULL,
                                       distribution = c("uniform", "exponential", "halfnormal", "lognormal"),
                                       n_sim = 1e6, n_max = 5, type_levels = NULL, seed = NULL,
                                       tk_max = NULL, plot = FALSE) {
  distribution <- match.arg(distribution)
  stopifnot(is.data.frame(data))

  tk_col    <- .rk_select(rlang::enquo(tk_col), data, "tk_col", one = TRUE, numeric = TRUE)
  type_cols <- .rk_select(rlang::enquo(type_col), data, "type_col")
  
  # Always drop missing/non-positive tk values, and (optionally) values
  # above tk_max, on the full dataset before any per-group splitting.
  tk_raw <- data[[tk_col]]
  keep <- !is.na(tk_raw) & tk_raw > 0
  if (!is.null(tk_max)) keep <- keep & tk_raw <= tk_max
  data <- data[keep, , drop = FALSE]
  if (nrow(data) == 0) {
    stop("No tk values > 0", if (!is.null(tk_max)) paste0(" and <= tk_max (", tk_max, ")") else "",
         " left in data after filtering.")
  }
  
  # Lower bound of the uniform distribution: estimated once on the full data
  global_min <- min(data[[tk_col]])
  
  # Short text describing a group (one-row data frame of type_col values)
  label_text <- function(label) {
    if (is.null(label)) return(NA_character_)
    vals <- vapply(label, function(x) as.character(x[1]), character(1))
    if (length(vals) == 1) unname(vals) else paste(names(vals), vals, sep = "=", collapse = ", ")
  }
  
  # Whole pipeline for one set of observed tk:
  # fit parameters -> simulate tk pairs -> rk -> classify -> proportions
  # `label` is NULL (no grouping) or a one-row data frame with the values of
  # the type_col columns identifying the group.
  null_for <- function(tk_values, label) {
    lab_txt <- label_text(label)
    
    # 1. estimate parameters from the observed tk
    params <- switch(distribution,
                     uniform     = list(min = global_min, max = max(tk_values)),
                     exponential = list(rate = 1 / mean(tk_values)),
                     halfnormal  = list(sigma = mean(tk_values) * sqrt(pi / 2)),
                     lognormal   = list(meanlog = mean(log(tk_values)), sdlog = stats::sd(log(tk_values)))
    )
    if (!all(is.finite(unlist(params)))) {
      stop("Could not estimate the '", distribution, "' parameters",
           if (!is.na(lab_txt)) paste0(" for group '", lab_txt, "'"),
           " (too few tk values?).")
    }
    
    # 2. random generator for tk under the null distribution
    draw <- function(n) {
      switch(distribution,
             uniform     = stats::runif(n, params$min, params$max),
             exponential = stats::rexp(n, rate = params$rate),
             halfnormal  = abs(stats::rnorm(n, mean = 0, sd = params$sigma)),
             lognormal   = stats::rlnorm(n, meanlog = params$meanlog, sdlog = params$sdlog)
      )
    }
    
    # 3. simulate pairs of consecutive intervals and compute rk
    if (!is.null(seed)) set.seed(seed)
    t1 <- draw(n_sim)
    t2 <- draw(n_sim)
    rk <- t1 / (t1 + t2)
    
    # optional plots: simulated vs observed tk, and simulated rk
    if (plot) {
      n_plot <- min(n_sim, 1e5)
      tk_sim <- t1[seq_len(n_plot)]
      rk_sim <- rk[seq_len(n_plot)]
      xmax <- max(stats::quantile(tk_sim, 0.99), stats::quantile(tk_values, 0.99))
      d_sim <- stats::density(tk_sim, from = 0, to = xmax)
      d_obs <- stats::density(tk_values, from = 0, to = xmax)
      ttl <- if (is.na(lab_txt)) "" else paste0(" (", lab_txt, ")")
      par_txt <- paste(names(params), signif(unlist(params), 3), sep = " = ", collapse = ", ")
      
      graphics::plot(d_sim, xlim = c(0, xmax), ylim = range(d_sim$y, d_obs$y), lwd = 2,
                     main = paste0("tk: null (", distribution, ") vs observed", ttl),
                     sub = par_txt, xlab = "tk", ylab = "density")
      graphics::lines(d_obs, lty = 2, col = "red")
      graphics::legend("topright", legend = c("simulated (null)", "observed"),
                       lty = c(1, 2), col = c("black", "red"), lwd = c(2, 1),
                       bty = "n", cex = 0.8)
      
      graphics::plot(stats::density(rk_sim, from = 0, to = 1), lwd = 2,
                     main = paste0("rk (simulated)", ttl),
                     xlab = "rk", ylab = "density")
    }
    
    # 4. classify rk and compute the proportion of each bin
    tab <- table(classify_rk(rk, n_max = n_max), useNA = "no")
    bins <- data.frame(
      bin = names(tab),
      normalization_constant = as.numeric(tab) / n_sim,
      stringsAsFactors = FALSE
    )
    par_df <- data.frame(
      parameter = names(params),
      value = as.numeric(unlist(params)),
      stringsAsFactors = FALSE
    )
    # prepend the grouping columns (same values repeated on every row)
    if (!is.null(label)) {
      bins   <- cbind(label[rep(1L, nrow(bins)), , drop = FALSE], bins)
      par_df <- cbind(label[rep(1L, nrow(par_df)), , drop = FALSE], par_df)
    }
    rownames(bins) <- NULL
    rownames(par_df) <- NULL
    list(bins = bins, params = par_df)
  }
  
  # collect bins + params from a list of null_for() results
  combine <- function(res_list) {
    out <- do.call(rbind, lapply(res_list, `[[`, "bins"))
    rownames(out) <- NULL
    par_all <- do.call(rbind, lapply(res_list, `[[`, "params"))
    rownames(par_all) <- NULL
    attr(out, "params") <- par_all
    out
  }
  
  if (length(type_cols) == 0) {
    # no grouping: a single null distribution over the whole dataset
    if (plot) {
      old_par <- graphics::par(mfrow = c(1, 2))
      on.exit(graphics::par(old_par), add = TRUE)
    }
    return(combine(list(null_for(data[[tk_col]], NULL))))
  }
  
  # rows with NA in a grouping column cannot be assigned to any group
  na_rows <- !stats::complete.cases(as.data.frame(data)[type_cols])
  if (any(na_rows)) {
    warning(sum(na_rows), " row(s) with NA in ", paste0("'", type_cols, "'", collapse = ", "),
            " dropped.")
    data <- data[!na_rows, , drop = FALSE]
    if (nrow(data) == 0) stop("No rows left in data after dropping NA in the type_col column(s).")
  }
  
  # groups: one row per combination of type_col values
  row_keys <- .rk_make_key(data, type_cols)
  if (is.null(type_levels)) {
    groups <- as.data.frame(data)[!duplicated(row_keys), type_cols, drop = FALSE]
  } else if (is.data.frame(type_levels)) {
    if (!all(type_cols %in% names(type_levels))) {
      stop("'type_levels' must have a column for each type_col column: ",
           paste(type_cols, collapse = ", "), ".")
    }
    groups <- as.data.frame(type_levels)[type_cols]
  } else if (length(type_cols) == 1) {
    groups <- stats::setNames(data.frame(type_levels, stringsAsFactors = FALSE), type_cols)
  } else {
    stop("With several type_col columns, 'type_levels' must be a data frame with one column for each of them (",
         paste(type_cols, collapse = ", "), ").")
  }
  rownames(groups) <- NULL
  group_keys <- .rk_make_key(groups, type_cols)
  
  if (plot) {
    old_par <- graphics::par(mfrow = c(min(nrow(groups), 3), 2))
    on.exit(graphics::par(old_par), add = TRUE)
  }
  
  results <- lapply(seq_len(nrow(groups)), function(i) {
    label <- groups[i, , drop = FALSE]
    tk_values <- data[[tk_col]][row_keys == group_keys[i]]
    if (length(tk_values) == 0) {
      stop("Group '", label_text(label), "' has no tk values > 0",
           if (!is.null(tk_max)) paste0(" and <= tk_max (", tk_max, ")") else "",
           " after filtering.")
    }
    if (length(tk_values) < 10) {
      warning("Group '", label_text(label), "' has only ", length(tk_values),
              " tk value(s): its null distribution is poorly estimated.")
    }
    null_for(tk_values, label)
  })
  
  combine(results)
}