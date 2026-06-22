#' pheatmap-style Heatmap Color Scale
#'
#' Discrete-tile heatmap fill scale using the classic blue-white-red
#' diverging palette inherited from \code{pheatmap::pheatmap}'s default
#' (\code{colorRampPalette(c("#0C6291", "white", "#A63446"))}). Best suited
#' for diverging data centered around zero (correlations, log fold changes,
#' z-scores). For strictly positive sequential data (0-1 scaled values,
#' counts), prefer \code{\link{scale_fill_sci_c}} with a sequential palette.
#'
#' @param n Integer; number of colors in the interpolated gradient
#'   (default: 100).
#' @param alpha Numeric in (0, 1]; transparency (default: 1).
#' @param reverse Logical; reverse the palette (default: FALSE).
#' @param ... Additional arguments passed to \code{ggplot2::scale_fill_gradientn()}
#'   (e.g., \code{values}, \code{limits}, \code{oob}, \code{na.value}).
#'
#' @return A ggplot scale object.
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#'
#' # Diverging data (z-scores)
#' df <- expand.grid(x = 1:5, y = 1:5)
#' df$z <- rnorm(nrow(df))
#'
#' ggplot(df, aes(x, y, fill = z)) +
#'   geom_tile() +
#'   scale_fill_heatmap()
#'
#' # Reverse direction
#' ggplot(df, aes(x, y, fill = z)) +
#'   geom_tile() +
#'   scale_fill_heatmap(reverse = TRUE)
#'
#' # Align midpoint to non-zero value via values (passed to gradientn)
#' ggplot(df, aes(x, y, fill = z)) +
#'   geom_tile() +
#'   scale_fill_heatmap(values = scales::rescale(c(min(df$z), 0, max(df$z))))
#' }
#'
scale_fill_heatmap <- function(n = 100L, alpha = 1, reverse = FALSE, ...) {

  if (!is.numeric(n) || length(n) != 1L || n < 2L || n != round(n)) {
    stop("'n' must be a single integer >= 2")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0L || alpha > 1L) {
    stop("'alpha' must be a single numeric value in (0, 1]")
  }
  if (!is.logical(reverse) || length(reverse) != 1L) {
    stop("'reverse' must be a single TRUE or FALSE")
  }

  colors <- grDevices::colorRampPalette(c("#0C6291", "white", "#A63446"))(n)

  if (alpha < 1L) {
    cols_rgb <- grDevices::col2rgb(colors)
    colors <- rgb(
      cols_rgb[1L, ], cols_rgb[2L, ], cols_rgb[3L, ],
      alpha = alpha * 255L,
      maxColorValue = 255L
    )
  }

  if (reverse) colors <- rev(colors)

  ggplot2::scale_fill_gradientn(colours = colors, ...)
}

#' @rdname scale_fill_heatmap
#' @export
scale_color_heatmap <- function(n = 100L, alpha = 1, reverse = FALSE, ...) {
  scale_fill_heatmap(n = n, alpha = alpha, reverse = reverse, ...)
}

#' @rdname scale_fill_heatmap
#' @export
scale_colour_heatmap <- scale_color_heatmap
