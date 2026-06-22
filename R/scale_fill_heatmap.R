#' pheatmap-style Heatmap Color Scale
#'
#' Discrete-tile heatmap fill scale using a blue-red diverging palette
#' inspired by \code{pheatmap::pheatmap}'s default. Middle color defaults
#' to light gray (\code{#F0F0F0}) instead of pure white so middle-valued
#' cells remain visible on white plot backgrounds. Pass
#' \code{mid_color = "white"} to restore the original pheatmap look.
#'
#' Best suited for diverging data centered around zero (correlations, log
#' fold changes, z-scores). For strictly positive sequential data (0-1
#' scaled values, counts), prefer \code{\link{scale_fill_sci_c}} with a
#' sequential palette.
#'
#' @param n Integer; number of colors in the interpolated gradient
#'   (default: 100).
#' @param alpha Numeric in (0, 1]; transparency (default: 1).
#' @param reverse Logical; reverse the palette (default: FALSE).
#' @param mid_color Middle color between low (\code{#0C6291}) and high
#'   (\code{#A63446}). Default: \code{"#F0F0F0"} (light gray, visible on
#'   white backgrounds). Use \code{"white"} for original pheatmap style,
#'   or any other color (e.g., \code{"yellow"} for RdYlBu style).
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
#' # Original pheatmap style (pure white middle)
#' ggplot(df, aes(x, y, fill = z)) +
#'   geom_tile() +
#'   scale_fill_heatmap(mid_color = "white")
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
scale_fill_heatmap <- function(n = 100L, alpha = 1, reverse = FALSE,
                               mid_color = "#F0F0F0", ...) {

  if (!is.numeric(n) || length(n) != 1L || n < 2L || n != round(n)) {
    stop("'n' must be a single integer >= 2")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0L || alpha > 1L) {
    stop("'alpha' must be a single numeric value in (0, 1]")
  }
  if (!is.logical(reverse) || length(reverse) != 1L) {
    stop("'reverse' must be a single TRUE or FALSE")
  }
  if (!is.character(mid_color) || length(mid_color) != 1L || is.na(mid_color)) {
    stop("'mid_color' must be a single character string")
  }

  colors <- grDevices::colorRampPalette(c("#0C6291", mid_color, "#A63446"))(n)

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
scale_color_heatmap <- function(n = 100L, alpha = 1, reverse = FALSE,
                                mid_color = "#F0F0F0", ...) {
  scale_fill_heatmap(n = n, alpha = alpha, reverse = reverse,
                     mid_color = mid_color, ...)
}

#' @rdname scale_fill_heatmap
#' @export
scale_colour_heatmap <- scale_color_heatmap
