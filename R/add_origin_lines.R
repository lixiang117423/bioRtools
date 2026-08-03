#' Add Dashed Reference Lines Through the Origin
#'
#' @description
#' Convenience helper returning a vertical and a horizontal dashed reference
#' line at (0, 0), the two-line combo commonly used on PCA / PCoA / RDA / t-SNE
#' score plots. Add it directly with `+`, e.g. `p + add_origin_lines()`.
#'
#' @param x_intercept Numeric position on the x axis for the vertical line
#'   (default: 0).
#' @param y_intercept Numeric position on the y axis for the horizontal line
#'   (default: 0).
#' @param linetype Character string for line type (default: "dashed").
#' @param color Character string for line color. `NULL` (default) inherits the
#'   ggplot2 / theme default, so the lines render exactly like a bare
#'   [ggplot2::geom_vline()] / [ggplot2::geom_hline()].
#' @param alpha Numeric transparency between 0 and 1. `NULL` (default)
#'   inherits the ggplot2 default.
#' @param linewidth Numeric line width. `NULL` (default) inherits the ggplot2
#'   default.
#' @param ... Additional arguments passed to both [ggplot2::geom_vline()] and
#'   [ggplot2::geom_hline()].
#'
#' @return A list of two ggplot2 layer objects (`geom_vline` then `geom_hline`)
#'   that can be added to a plot with `+`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#'
#' p <- ggplot(iris, aes(Sepal.Length, Sepal.Width)) +
#'   geom_point() +
#'   add_origin_lines()
#' print(p)
#' }
#'
#' @author Xiang LI <lixiang117423@gmail.com>
add_origin_lines <- function(x_intercept = 0, y_intercept = 0,
                             linetype = "dashed", color = NULL,
                             alpha = NULL, linewidth = NULL, ...) {
  # 1. Input validation
  if (!is.numeric(x_intercept) || length(x_intercept) != 1 || is.na(x_intercept)) {
    stop("'x_intercept' must be a single numeric value, got: ",
         paste(class(x_intercept), collapse = "/"))
  }
  if (!is.numeric(y_intercept) || length(y_intercept) != 1 || is.na(y_intercept)) {
    stop("'y_intercept' must be a single numeric value, got: ",
         paste(class(y_intercept), collapse = "/"))
  }
  if (!is.character(linetype) || length(linetype) != 1 || is.na(linetype)) {
    stop("'linetype' must be a single character value, got: ",
         paste(class(linetype), collapse = "/"))
  }
  if (!is.null(color) &&
      (!is.character(color) || length(color) != 1 || is.na(color))) {
    stop("'color' must be a single character value or NULL, got: ",
         paste(class(color), collapse = "/"))
  }
  if (!is.null(alpha) &&
      (!is.numeric(alpha) || length(alpha) != 1 || is.na(alpha) ||
       alpha < 0 || alpha > 1)) {
    stop("'alpha' must be a single numeric value between 0 and 1 or NULL, got: ",
         paste(class(alpha), collapse = "/"))
  }
  if (!is.null(linewidth) &&
      (!is.numeric(linewidth) || length(linewidth) != 1 || is.na(linewidth) ||
       linewidth <= 0)) {
    stop("'linewidth' must be a single positive numeric value or NULL, got: ",
         paste(class(linewidth), collapse = "/"))
  }

  # 2. Build layers; NULL style args are omitted so ggplot2 defaults apply
  vline_args <- list(xintercept = x_intercept, linetype = linetype)
  hline_args <- list(yintercept = y_intercept, linetype = linetype)
  if (!is.null(color)) {
    vline_args$color <- color
    hline_args$color <- color
  }
  if (!is.null(alpha)) {
    vline_args$alpha <- alpha
    hline_args$alpha <- alpha
  }
  if (!is.null(linewidth)) {
    vline_args$linewidth <- linewidth
    hline_args$linewidth <- linewidth
  }
  extra_args <- list(...)
  if (length(extra_args) > 0) {
    vline_args <- c(vline_args, extra_args)
    hline_args <- c(hline_args, extra_args)
  }

  list(
    do.call(ggplot2::geom_vline, vline_args),
    do.call(ggplot2::geom_hline, hline_args)
  )
}
