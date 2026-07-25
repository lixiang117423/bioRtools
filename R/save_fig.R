#' Save a Figure for Journal Submission
#'
#' @description Thin wrapper around [ggplot2::ggsave()] with publication-sensible
#'   defaults. Vector formats use font-embedding devices so text stays editable
#'   (not outlined) when opened in Illustrator — what the Nature figure guide
#'   asks for (embedded TrueType fonts, vector artwork). Use this for the
#'   final export step instead of [ggplot2::ggsave()].
#'
#'   - `pdf` via [grDevices::cairo_pdf] (fonts embedded, text stays editable),
#'     `png` via [ragg::agg_png], `tiff` via [ragg::agg_tiff], `jpeg` via
#'     [ragg::agg_jpeg] (high-quality anti-aliasing).
#'   - `eps` via [grDevices::cairo_ps] (font-embedded PostScript).
#'   - `svg` via the ggsave default device.
#'
#'   Pass `formats = c("pdf", "png")` to emit both a vector (for editing) and a
#'   raster (for preview) copy in one call.
#'
#' @param plot A ggplot object. Defaults to [ggplot2::last_plot()].
#' @param filename Output path. If `formats` is `NULL`, its extension picks the
#'   format. If `formats` is given, the extension (if any) is stripped and one
#'   file is written per format, e.g. `"fig"` + `c("pdf","png")` -> `fig.pdf`,
#'   `fig.png`.
#' @param formats Character vector of formats to write (e.g. `c("pdf", "png")`).
#'   `NULL` (default) writes a single file using `filename`'s extension.
#' @param width,height,in Plot dimensions. `NA` (default) lets ggsave pick from
#'   the current device.
#' @param dpi Dots per inch for raster output (default 300; raise for Nature's
#'   450-dpi preference).
#' @param scale Multiplicative scaling factor (default 1).
#' @param ... Passed to [ggplot2::ggsave()].
#'
#' @return `filename` invisibly (the base name when `formats` is given).
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' p <- ggplot(iris, aes(Sepal.Length, Sepal.Width, colour = Species)) +
#'   geom_point() + bioRtools::theme_prism()
#'
#' # single vector PDF (fonts embedded, not outlined)
#' save_fig(p, "fig1.pdf")
#'
#' # both vector + raster in one call
#' save_fig(p, "fig1", formats = c("pdf", "png"), width = 180, height = 70,
#'          units = "mm", dpi = 450)
#' }
save_fig <- function(plot = ggplot2::last_plot(),
                     filename,
                     formats = NULL,
                     width = NA,
                     height = NA,
                     units = "in",
                     dpi = 300,
                     scale = 1,
                     ...) {
  if (is.null(formats)) {
    formats <- tolower(tools::file_ext(filename))
    if (!nzchar(formats)) {
      stop("'filename' has no extension and 'formats' is NULL; ",
        "give an extension (e.g. \"fig.pdf\") or pass ",
        "formats = c(\"pdf\", \"png\")",
        call. = FALSE
      )
    }
  }

  base <- sub("\\.[^.]*$", "", filename)

  # Vector formats use font-embedding devices (text stays editable, not outlined);
  # raster formats use ragg for high-quality anti-aliasing.
  devices <- list(
    pdf   = grDevices::cairo_pdf,
    png   = ragg::agg_png,
    tiff  = ragg::agg_tiff,
    jpeg  = ragg::agg_jpeg,
    eps   = grDevices::cairo_ps
  )
  default_device_fmts <- "svg"

  for (fmt in formats) {
    key <- if (fmt == "jpg") "jpeg" else fmt
    device <- devices[[key]]
    if (is.null(device) && !(fmt %in% default_device_fmts)) {
      stop("Unsupported format: '", fmt,
        "'. Supported: pdf, png, tiff, eps, svg, jpg",
        call. = FALSE
      )
    }
    args <- list(
      filename = paste0(base, ".", fmt), plot = plot,
      width = width, height = height, units = units,
      dpi = dpi, scale = scale, ...
    )
    if (!is.null(device)) args$device <- device
    do.call(ggplot2::ggsave, args)
  }

  invisible(filename)
}
