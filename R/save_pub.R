#' Export Plot as Publication Bundle (SVG + PDF + TIFF)
#'
#' Writes a ggplot object to multiple publication-ready files in one call:
#' editable SVG for post-hoc Illustrator/Inkscape work, embedded-font PDF
#' for journal submission portals, and high-DPI TIFF for journals that
#' require bitmap formats. Defaults match Nature-family double-column
#' dimensions (183 x 120 mm) and Arial typography.
#'
#' @param plot ggplot object (or any object printable by \code{print()}).
#' @param filename Character. Output path **without extension**; each
#'   format's extension is appended automatically. If a path with an
#'   extension is supplied, it is stripped.
#' @param width_mm Numeric. Figure width in millimeters (default: 183,
#'   Nature double-column width). Use 89 for single-column figures.
#' @param height_mm Numeric. Figure height in millimeters (default: 120).
#' @param dpi Numeric. Resolution for TIFF/PNG raster output
#'   (default: 600). 300 is acceptable for line art; 600 for dense panels.
#' @param formats Character vector. Formats to export; subset of
#'   \code{c("svg", "pdf", "tiff")} (default: all three).
#' @param base_family Character. Font family embedded in PDF
#'   (default: "Arial"). Must be a font registered on the system.
#'
#' @return Invisibly returns \code{filename}. The full paths of written
#'   files are attached as attribute \code{"files"}.
#'
#' @details
#' \strong{Why three formats?}
#' \itemize{
#'   \item \strong{SVG} via \code{svglite} — text stored as \code{<text>}
#'     nodes (not bezier outlines), so labels remain editable in vector
#'     editors. Always save SVG first.
#'   \item \strong{PDF} via \code{grDevices::cairo_pdf} — embeds fonts
#'     natively, accepted by most journal submission systems.
#'   \item \strong{TIFF} via \code{ragg::agg_tiff} — anti-aliased bitmap
#'     for journals that require TIFF (some legacy Nature protocols).
#' }
#'
#' @note All mm dimensions are converted to inches internally
#'   (\code{mm / 25.4}).
#'
#' @seealso \code{\link{theme_nature}} for the matching theme,
#'   \code{\link{write_data}} for single-format export
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#'
#' p <- ggplot(iris, aes(Species, Sepal.Length, fill = Species)) +
#'   geom_boxplot() +
#'   bioRtools::theme_nature()
#'
#' # Single-column Nature figure
#' bioRtools::save_pub(p, "~/figures/fig1", width_mm = 89)
#' # writes: fig1.svg, fig1.pdf, fig1.tiff
#'
#' # PDF only
#' bioRtools::save_pub(p, "~/figures/fig1", formats = "pdf")
#' }
save_pub <- function(plot, filename, width_mm = 183, height_mm = 120,
                     dpi = 600, formats = c("svg", "pdf", "tiff"),
                     base_family = "Arial") {

  if (missing(plot)) stop("'plot' is required")
  if (!is.character(filename) || length(filename) != 1) {
    stop("'filename' must be a single character string")
  }
  for (arg in c("width_mm", "height_mm")) {
    val <- get(arg)
    if (!is.numeric(val) || length(val) != 1 || val <= 0) {
      stop("'", arg, "' must be a single positive number")
    }
  }
  if (!is.numeric(dpi) || length(dpi) != 1 || dpi <= 0) {
    stop("'dpi' must be a single positive number")
  }
  formats <- match.arg(formats, c("svg", "pdf", "tiff"), several.ok = TRUE)

  if (grepl("\\.[a-zA-Z0-9]+$", filename)) {
    filename <- sub("\\.[a-zA-Z0-9]+$", "", filename)
    warning("'filename' had an extension; stripped to '", filename, "'")
  }

  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  written <- character(0)

  if ("svg" %in% formats) {
    out <- paste0(filename, ".svg")
    svglite::svglite(out, width = width_in, height = height_in)
    print(plot)
    grDevices::dev.off()
    written <- c(written, out)
  }

  if ("pdf" %in% formats) {
    out <- paste0(filename, ".pdf")
    grDevices::cairo_pdf(out, width = width_in, height = height_in,
                         family = base_family)
    print(plot)
    grDevices::dev.off()
    written <- c(written, out)
  }

  if ("tiff" %in% formats) {
    out <- paste0(filename, ".tiff")
    ragg::agg_tiff(out, width = width_in, height = height_in,
                   units = "in", res = dpi)
    print(plot)
    grDevices::dev.off()
    written <- c(written, out)
  }

  attr(filename, "files") <- written
  invisible(filename)
}
