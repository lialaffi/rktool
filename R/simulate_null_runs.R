#' Repeat the null distribution simulation multiple times
#'
#' Calls \code{\link{compute_null_distributions}} \code{runs} times,
#' generating that many independent sets of normalization constants. Used to
#' assess how much a result (e.g. a significance test) depends on the
#' specific Monte Carlo simulation used for normalization, rather than being
#' an artifact of random sampling. The null distributions follow the
#' approach of Jadoul et al. (2025).
#'
#' @inheritParams compute_null_distributions
#' @param runs Number of independent simulations to generate (default 100).
#' @param n_sim Number of simulated interval pairs in each run (default 1e5).
#'   Large values combined with many \code{runs} can be slow.
#' @param seed Optional seed for reproducibility: run \code{i} uses
#'   \code{seed + i}, so runs differ from each other but the whole list is
#'   reproducible.
#' @param verbose If \code{TRUE} (default), prints a progress message every
#'   100 simulations.
#' @return A list of length \code{runs}, each element a data frame as
#'   returned by \code{\link{compute_null_distributions}}.
#'
#' @details
#' Warnings raised by \code{compute_null_distributions()} (e.g. groups with
#' very few \code{tk} values) depend only on the data, not on the
#' simulation, so they are shown for the first run only and suppressed in
#' the following ones.
#'
#' @references
#' Jadoul, Y., Tufarelli, T., Coissac, C., Gamba, M., & Ravignani, A.
#' (2025). Hidden assumptions of integer ratio analyses in bioacoustics and
#' music. \emph{Annals of the New York Academy of Sciences}, 1553(1),
#' 363-378.
#'
#' @examples
#' y <- calculate_tk(data_rktool, onset_col = onset, by = c(file_id, ind_id))
#'
#' # 10 runs with few simulated pairs, for a quick example (in a real
#' # analysis use at least the defaults, runs = 100 and n_sim = 1e5)
#' norm_runs <- simulate_null_runs(y, tk_col = tk, type_col = sex,
#'                                 distribution = "lognormal",
#'                                 n_sim = 1e4, runs = 10, seed = 1)
#' length(norm_runs)
#' head(norm_runs[[1]])
#'
#' @export
simulate_null_runs <- function(data, tk_col, type_col = NULL,
                               distribution = c("uniform", "exponential", "halfnormal", "lognormal"),
                               n_sim = 1e5, n_max = 5, runs = 100,
                               type_levels = NULL, seed = NULL, verbose = TRUE,
                               tk_max = NULL) {
  distribution <- match.arg(distribution)
  if (length(runs) != 1 || is.na(runs) || runs < 1) stop("'runs' must be an integer >= 1.")

  one_run <- function(i) {
    run_seed <- if (!is.null(seed)) seed + i else NULL
    compute_null_distributions(
      data, tk_col = {{ tk_col }}, type_col = {{ type_col }},
      distribution = distribution, n_sim = n_sim, n_max = n_max, tk_max = tk_max,
      type_levels = type_levels, seed = run_seed
    )
  }

  lapply(seq_len(runs), function(i) {
    res <- if (i == 1) one_run(i) else suppressWarnings(one_run(i))
    if (verbose && i %% 100 == 0) message("Completed simulation ", i, "/", runs)
    res
  })
}