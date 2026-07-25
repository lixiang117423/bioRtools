#' Wong (Okabe-Ito) Colour-Blind-Friendly Palette
#'
#' @description Qualitative 8-colour palette from Wong, B. (2011) *Points of
#'   view: Colour blindness.* Nature Methods 8, 441 — the accessible palette
#'   recommended by the Nature figure guide. The eight colours are mutually
#'   distinguishable for viewers with deuteranopia, protanopia, and tritanopia.
#'
#'   `wong_palette()` returns the colours; `scale_colour_wong()` /
#'   `scale_color_wong()` / `scale_fill_wong()` are drop-in ggplot2 discrete
#'   scales usable on any plot (`+ bioRtools::scale_fill_wong()`).
#'
#' @param n Integer, number of colours. If `NULL` (default) all 8 are returned.
#'   For `1 <= n <= 8` the first `n` colours are returned; for `n > 8` the
#'   palette is recycled (interpolation is deliberately avoided, as it can
#'   produce colours outside the colourblind-safe set). `n <= 0` returns
#'   `character(0)`.
#' @param ... Passed to [ggplot2::discrete_scale()].
#'
#' @return `wong_palette()` returns a character vector of hex colours. The
#'   `scale_*_wong()` functions return a ggplot2 discrete scale object.
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @name wong_palette
#'
#' @examples
#' library(ggplot2)
#'
#' wong_palette()
#' wong_palette(3)
#'
#' ggplot(iris, aes(Species, Sepal.Length, fill = Species)) +
#'   geom_boxplot() +
#'   bioRtools::scale_fill_wong() +
#'   bioRtools::theme_prism()
NULL

# Wong (2011) Nature Methods 8:441 — eight-colour colourblind-safe set.
wong_colors <- c(
  "#000000", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
)

#' @rdname wong_palette
#' @export
wong_palette <- function(n = NULL) {
  if (is.null(n)) return(wong_colors)
  if (n <= 0) return(character(0))
  if (n <= length(wong_colors)) wong_colors[seq_len(n)] else rep(wong_colors, length.out = n)
}

#' @rdname wong_palette
#' @export
scale_colour_wong <- function(...) {
  ggplot2::discrete_scale("colour", palette = wong_palette, ...)
}

#' @rdname wong_palette
#' @export
scale_color_wong <- scale_colour_wong

#' @rdname wong_palette
#' @export
scale_fill_wong <- function(...) {
  ggplot2::discrete_scale("fill", palette = wong_palette, ...)
}
