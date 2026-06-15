#' Nature/High-Impact Journal Publication Theme
#'
#' A minimalist ggplot2 theme tuned for Nature-family and other high-impact
#' journal figures at publication column width (89 mm single column, 183 mm
#' double column). Built on \code{theme_classic()} with thin axis lines,
#' Arial typography, and a tight font hierarchy calibrated for dense
#' multi-panel composites.
#'
#' @param base_size Numeric. Base font size in points (default: 6.5,
#'   matching Nature final-text regime for dense composites). Increase to
#'   8-9 for less dense single-panel figures.
#' @param base_family Character. Font family (default: "Arial", the Nature
#'   standard). Use "Helvetica" as an equivalent on systems without Arial.
#' @param line_width Numeric. Axis line and tick width in points
#'   (default: 0.35, the Nature-standard thin line).
#'
#' @return A ggplot2 theme object.
#'
#' @details
#' Design rationale (from Nature 2026 sampled figures):
#' \itemize{
#'   \item \strong{Thin lines}: 0.35pt axis/tick lines keep panels from
#'     dominating dense composites.
#'   \item \strong{Font hierarchy}: plot.title (+0.5) > axis.title (base)
#'     > axis.text (-0.5) > legend.title (-0.3) > legend.text (-0.7),
#'     creating clear hierarchy without large jumps.
#'   \item \strong{No grid}: sparse y-ticks guide the eye; grids add visual
#'     noise at small scales.
#'   \item \strong{Frameless legend}: legends sit quietly beside data.
#' }
#'
#' @seealso \code{\link{save_pub}} for publication export (SVG/PDF/TIFF),
#'   \code{\link{theme_prism}} for preview-style theme
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' library(ggplot2)
#'
#' ggplot(iris, aes(x = Species, y = Sepal.Length, fill = Species)) +
#'   geom_boxplot() +
#'   bioRtools::theme_nature()
#'
#' # Larger base for single-panel figures
#' ggplot(mtcars, aes(x = wt, y = mpg, color = cyl)) +
#'   geom_point() +
#'   bioRtools::theme_nature(base_size = 9)
theme_nature <- function(base_size = 6.5, base_family = "Arial",
                         line_width = 0.35) {

  if (!is.numeric(base_size) || length(base_size) != 1 || base_size <= 0) {
    stop("'base_size' must be a single positive number")
  }
  if (!is.character(base_family) || length(base_family) != 1) {
    stop("'base_family' must be a single character string")
  }
  if (!is.numeric(line_width) || length(line_width) != 1 || line_width <= 0) {
    stop("'line_width' must be a single positive number")
  }

  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = line_width, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = line_width, colour = "black"),
      axis.title = ggplot2::element_text(size = base_size),
      axis.text = ggplot2::element_text(size = base_size - 0.5, colour = "black"),
      legend.title = ggplot2::element_text(size = base_size - 0.3),
      legend.text = ggplot2::element_text(size = base_size - 0.7),
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = base_size - 0.3, face = "bold"),
      strip.background = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = base_size + 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = base_size),
      plot.caption = ggplot2::element_text(size = base_size - 1),
      panel.grid = ggplot2::element_blank()
    )
}
