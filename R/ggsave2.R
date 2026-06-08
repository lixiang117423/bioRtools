#' Save a ggplot with Sensible Defaults
#'
#' A thin wrapper around \code{\link[ggplot2]{ggsave}} with \code{width = 8}
#' and \code{height = 6} as defaults. All other parameters are passed through
#' unchanged and support full autocomplete in RStudio.
#'
#' @param filename File name to create on disk.
#' @param plot Plot to save, defaults to last plot displayed.
#' @param device Device to use. Can be a function or one of "eps", "ps",
#'   "tex", "pdf", "jpeg", "tiff", "png", "bmp", "svg", or "wmf".
#' @param path Path to save plot to (combined with filename).
#' @param scale Multiplicative scaling factor (default: 1).
#' @param width Plot width in \code{units} (default: 8).
#' @param height Plot height in \code{units} (default: 6).
#' @param units Units for width and height: "in", "cm", "mm", or "px".
#' @param dpi Plot resolution, or "retina"/"print"/"screen".
#' @param limitsize When TRUE (default), stops if dimensions exceed 50 inches.
#' @param bg Background colour, e.g. "white".
#' @param ... Other arguments passed to graphics device.
#'
#' @return Invisibly returns the input \code{plot}.
#' @export
#' @seealso \code{\link[ggplot2]{ggsave}}
#'
#' @examples
#' \dontrun{
#' p <- ggplot2::ggplot(mtcars, ggplot2::aes(mpg, disp)) + ggplot2::geom_point()
#' ggsave2("plot.pdf", p)
#' }
ggsave2 <- function(filename, plot = ggplot2::last_plot(),
                    device = NULL, path = NULL, scale = 1,
                    width = 8, height = 6, units = "in",
                    dpi = 600, limitsize = TRUE, bg = NULL, ...) {
  ggplot2::ggsave(
    filename  = filename,
    plot      = plot,
    device    = device,
    path      = path,
    scale     = scale,
    width     = width,
    height    = height,
    units     = units,
    dpi       = dpi,
    limitsize = limitsize,
    bg        = bg,
    ...
  )
}
