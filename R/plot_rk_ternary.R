#' Build triplets of consecutive intervals (tk1, tk2, tk3)
#'
#' Internal helper for \code{\link{plot_rk_ternary}}. Within each group
#' builds, for every interval \code{tk}, the triplet
#' \code{(tk1, tk2, tk3) = (tk, tk+1, tk+2)}. A triplet is valid only if the
#' three intervals are finite, > 0, consecutive, and (if \code{tk_max} is
#' supplied) all \code{<= tk_max}; otherwise it is \code{NA}.
#'
#' @param data A data frame.
#' @param tk_col Column with the intervals (bare or string).
#' @param by Optional grouping column(s) (bare, string or \code{c(...)}).
#'   Triplets never span two groups.
#' @param order_col Optional column used to sort rows within each group
#'   (bare or string). If \code{NULL} and \code{int_col} is given, rows are
#'   sorted by \code{int_col}; otherwise the row order is kept.
#' @param int_col Optional column with the progressive interval number
#'   within the group (bare or string). If supplied, the three numbers must
#'   be consecutive (\code{n, n+1, n+2}).
#' @param tk_max Optional upper threshold on the intervals.
#' @param drop_incomplete Logical, if \code{TRUE} (default) rows without a
#'   valid triplet are removed.
#' @return An ungrouped data frame with the original columns plus
#'   \code{tk1}, \code{tk2}, \code{tk3}.
#' @keywords internal
.calculate_tk_triplets <- function(data, tk_col, by = NULL, order_col = NULL,
                                   int_col = NULL, tk_max = NULL,
                                   drop_incomplete = TRUE) {
  tk_name    <- .rk_select(rlang::enquo(tk_col), data, "tk_col", one = TRUE, numeric = TRUE)
  by_cols    <- .rk_select(rlang::enquo(by), data, "by")
  order_cols <- .rk_select(rlang::enquo(order_col), data, "order_col")
  int_name   <- .rk_select(rlang::enquo(int_col), data, "int_col")
  if (length(int_name) > 1) stop("'int_col' must select at most one column.")
  
  if (any(c("tk1", "tk2", "tk3") %in% names(data))) {
    warning("Columns tk1/tk2/tk3 already exist and will be overwritten.")
  }
  
  sort_cols <- if (length(order_cols) > 0) order_cols else int_name
  if (length(c(by_cols, sort_cols)) > 0) {
    data <- dplyr::arrange(data, dplyr::across(dplyr::all_of(c(by_cols, sort_cols))))
  }
  
  data <- dplyr::group_by(data, dplyr::across(dplyr::all_of(by_cols)))
  data <- dplyr::mutate(
    data,
    tk1 = .data[[tk_name]],
    tk2 = dplyr::lead(.data[[tk_name]], 1),
    tk3 = dplyr::lead(.data[[tk_name]], 2),
    .ok = is.finite(.data$tk1) & is.finite(.data$tk2) & is.finite(.data$tk3) &
      .data$tk1 > 0 & .data$tk2 > 0 & .data$tk3 > 0
  )
  if (!is.null(tk_max)) {
    data <- dplyr::mutate(
      data,
      .ok = .data$.ok & .data$tk1 <= tk_max & .data$tk2 <= tk_max & .data$tk3 <= tk_max
    )
  }
  if (length(int_name) == 1) {
    data <- dplyr::mutate(
      data,
      .ok = .data$.ok &
        (dplyr::lead(.data[[int_name]], 1) - .data[[int_name]] == 1) &
        (dplyr::lead(.data[[int_name]], 2) - .data[[int_name]] == 2)
    )
  }
  data <- dplyr::ungroup(data)
  
  ok <- !is.na(data$.ok) & data$.ok
  data$tk1[!ok] <- NA_real_
  data$tk2[!ok] <- NA_real_
  data$tk3[!ok] <- NA_real_
  data$.ok <- NULL
  if (drop_incomplete) data <- data[ok, , drop = FALSE]
  data
}


#' Ternary plot of triplets of consecutive intervals
#'
#' Starting from a table of intervals \code{tk}, builds the triplets of
#' consecutive intervals \code{(tk1, tk2, tk3) = (tk, tk+1, tk+2)} within
#' each group (checking that they are consecutive and that all three respect
#' \code{tk_max}) and draws them in a ternary diagram. ggtern automatically
#' rescales the three values to proportions of their sum. Density is shown
#' as filled contour polygons (darker = denser), as points, or both.
#'
#' Requires the \pkg{ggtern} package (and \pkg{ggplot2}).
#'
#' @param data A data frame with one row per interval.
#' @param tk_col Column with the intervals (bare or string).
#' @param by Optional grouping column(s) (bare, string or \code{c(...)}),
#'   e.g. \code{c(file_id, ind_id)}. Triplets are built only within each
#'   combination of these columns. If \code{NULL} the whole dataset is a
#'   single group.
#' @param order_col Optional column used to sort rows within each group
#'   (bare or string), e.g. \code{onset}.
#' @param int_col Optional column with the progressive interval number
#'   within the group (bare or string). If supplied, the three intervals of
#'   a triplet must have consecutive numbers, so intervals missing from the
#'   data are not "glued" together.
#' @param tk_max Optional upper threshold: a triplet is used only if all
#'   three intervals are \code{<= tk_max}. Default \code{NULL}.
#' @param geom \code{"density"} (default), \code{"point"} or \code{"both"}.
#' @param facet_by Optional column(s) to facet by (bare, string or
#'   \code{c(...)}), e.g. \code{sex}.
#' @param color_by Optional column mapped to colour/fill (bare or string),
#'   e.g. \code{sex}. If \code{NULL} (default) a fixed \code{fill} is used.
#' @param bins Number of density levels (default 8).
#' @param fill Fill colour when \code{color_by} is \code{NULL}.
#' @param point_alpha Transparency of the points (default 0.1).
#' @param labels Vertex labels (top, left, right). Default
#'   \code{c("tk1", "tk2", "tk3")}.
#' @param title Optional plot title.
#'
#' @return A \code{ggtern} / \code{ggplot} object.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("ggtern", quietly = TRUE)) {
#'   y <- calculate_tk(rk_example, onset_col = onset, by = c(file_id, ind_id))
#'
#'   # triplets within each recording x individual, intervals <= 2 s
#'   print(plot_rk_ternary(y, tk_col = tk, by = c(file_id, ind_id),
#'                         order_col = onset, tk_max = 2, facet_by = sex))
#'
#'   print(plot_rk_ternary(y, tk_col = tk, by = c(file_id, ind_id),
#'                         order_col = onset, tk_max = 2, geom = "both",
#'                         color_by = sex))
#' }
#' }
#'
#' @export
plot_rk_ternary <- function(data, tk_col, by = NULL, order_col = NULL,
                            int_col = NULL, tk_max = NULL,
                            geom = c("density", "point", "both"),
                            facet_by = NULL, color_by = NULL,
                            bins = 8, fill = "yellow", point_alpha = 0.1,
                            labels = c("tk1", "tk2", "tk3"), title = NULL) {
  geom <- match.arg(geom)
  if (!requireNamespace("ggtern", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_rk_ternary() requires the 'ggtern' and 'ggplot2' packages to be installed.")
  }
  stopifnot(is.data.frame(data), length(labels) == 3)

  # facet / color columns are resolved on the original data (before the
  # triplets are built, which keeps all the original columns)
  facet_cols <- .rk_select(rlang::enquo(facet_by), data, "facet_by")
  color_cols <- .rk_select(rlang::enquo(color_by), data, "color_by")
  if (length(color_cols) > 1) stop("'color_by' must select at most one column.")
  
  # triplets of consecutive tk (sub-function)
  data <- .calculate_tk_triplets(
    data, tk_col = {{ tk_col }}, by = {{ by }}, order_col = {{ order_col }},
    int_col = {{ int_col }}, tk_max = tk_max, drop_incomplete = TRUE
  )
  if (nrow(data) == 0) stop("No valid triplets of consecutive tk to plot.")
  
  aes_args <- list(x = quote(tk1), y = quote(tk2), z = quote(tk3))
  if (length(color_cols) == 1) {
    aes_args$colour <- rlang::sym(color_cols)
    aes_args$fill   <- rlang::sym(color_cols)
  }
  p <- ggtern::ggtern(data, do.call(ggplot2::aes, aes_args))
  
  if (geom %in% c("density", "both")) {
    level <- NULL  # computed by stat_density_tern(); avoids a check NOTE
    dens_args <- list(
      mapping = ggplot2::aes(alpha = ggplot2::after_stat(level)),
      geom = "polygon", bins = bins
    )
    if (length(color_cols) == 0) {
      dens_args$colour <- "black"
      dens_args$fill   <- fill
    }
    p <- p + do.call(ggtern::stat_density_tern, dens_args) +
      ggplot2::guides(alpha = "none")
  }
  if (geom %in% c("point", "both")) {
    pt_args <- list(alpha = point_alpha, size = 0.8)
    if (length(color_cols) == 0) pt_args$colour <- "grey30"
    p <- p + do.call(ggplot2::geom_point, pt_args)
  }
  
  if (length(facet_cols) > 0) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(!!!rlang::syms(facet_cols)))
  }
  
  p +
    ggplot2::labs(title = title, T = labels[1], L = labels[2], R = labels[3]) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text  = ggplot2::element_text(size = 8),
      axis.title = ggplot2::element_text(size = 9),
      plot.title = ggplot2::element_text(size = 15),
      strip.text = ggplot2::element_text(size = 12)
    )
}

