# Bioinfo colour palettes =====================================================
# 10 base hues, each in a dark/light pair; dark tones are spent first so
# small group counts get maximum contrast. Sequential ramps reuse the same
# hue names so a group shown in categorical "blue" can keep that hue in
# continuous layers.

bioinfo_cat20 <- c(
  blue_d   = "#497EAD", blue_l   = "#B1D3F0",
  red_d    = "#B85449", red_l    = "#F0B5AF",
  teal_d   = "#479E8D", teal_l   = "#ADE6DA",
  purple_d = "#795F9E", purple_l = "#CDBCE6",
  orange_d = "#D98E4C", orange_l = "#F0CAAA",
  pink_d   = "#D17997", pink_l   = "#F0C2D1",
  gold_d   = "#CCA752", gold_l   = "#F0DCAF",
  navy_d   = "#30426B", navy_l   = "#8694B2",
  olive_d  = "#7F9456", olive_l  = "#CDDBB2",
  brown_d  = "#855D49", brown_l  = "#CCB0A3"
)

bioinfo_cat10 <- unname(bioinfo_cat20[seq(1, 20, by = 2)])

bioinfo_seq_anchors <- list(
  blue   = c("#EBF3FA", "#2F5A80"),
  red    = c("#FAECEB", "#8A3930"),
  teal   = c("#EBFAF7", "#2D7063"),
  purple = c("#F1EBFA", "#533E70"),
  orange = c("#FAF2EB", "#AB6B33"),
  pink   = c("#FAEBF0", "#A35670"),
  gold   = c("#FAF5EB", "#9E7F37"),
  navy   = c("#EBEFFA", "#1F2C4C"),
  olive  = c("#F5FAEB", "#566636"),
  brown  = c("#FAF0EB", "#573A2B")
)

bioinfo_div_anchors <- list(
  soft   = c("#2C5F8A", "#7FA8D9", "#F5F5F5", "#F2A38B", "#B23A3A"),
  strong = c("#08519C", "#F7F7F7", "#B2182B")
)

bioinfo_multiseq_anchors <- c("#3A1F5D", "#2E6E8E", "#3FA796", "#F2E94E")

# Ramps are meaningless with fewer than two colours.
bioinfo_check_n <- function(n, min = 2L) {
  if (!is.numeric(n) || length(n) != 1L || is.na(n)) {
    stop("'n' must be a single number, got: ", paste(class(n), collapse = "/"))
  }
  if (n < min || n != round(n)) {
    stop("'n' must be a whole number >= ", min)
  }
  invisible(n)
}

bioinfo_check_alpha <- function(alpha) {
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
    alpha <= 0 || alpha > 1) {
    stop("'alpha' must be a single numeric value in (0, 1]")
  }
  invisible(alpha)
}

bioinfo_apply_alpha <- function(cols, alpha) {
  if (alpha == 1) return(cols)
  mat <- grDevices::col2rgb(cols)
  grDevices::rgb(
    mat[1L, ], mat[2L, ], mat[3L, ],
    alpha = round(alpha * 255L), maxColorValue = 255L
  )
}

bioinfo_style_cols <- function(cols, alpha = 1, reverse = FALSE) {
  if (reverse) cols <- rev(cols)
  bioinfo_apply_alpha(cols, alpha)
}

bioinfo_discrete_scale <- function(aesthetic, alpha, ...) {
  ggplot2::discrete_scale(
    aesthetic,
    palette = function(n) bioinfo_apply_alpha(get_cat_colors(n), alpha),
    ...
  )
}

#' Bioinfo colour palette getters
#'
#' Colour schemes distilled from Nature / Science / Nature Communications
#' figure styles. `get_cat_colors()` returns categorical colours; the other
#' getters return colour ramps for sequential, diverging, and multi-hue
#' (viridis-like) sequential data. Use [preview_bioinfo_palettes()] to view
#' every palette at once.
#'
#' @param n Number of colours. `get_cat_colors()` accepts any positive whole
#'   number: n <= 10 returns dark tones first, n <= 20 interleaves light
#'   tones, larger values interpolate. The ramp getters require n >= 2.
#' @param name Sequential hue family: one of "blue", "red", "teal",
#'   "purple", "orange", "pink", "gold", "navy", "olive", "brown" — the same
#'   hue names used by the categorical palette.
#' @param palette Diverging variant: "soft" (five-step blue-white-red) or
#'   "strong" (saturated three-step blue-white-red).
#'
#' @return An unnamed character vector of hex colours. The diverging ramps
#'   have a near-white midpoint, suited to z-scores, log2 fold changes, and
#'   correlation matrices; the multi-hue ramp is for intensity-only heatmaps.
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @name bioinfo_colors
#'
#' @examples
#' get_cat_colors(3)
#' get_seq_colors("teal", 5)
#' get_div_colors(5)
#' get_multiseq_colors(4)
#'
#' \dontrun{
#' # Diverging heatmap centred on 0
#' pheatmap::pheatmap(mat, color = get_div_colors(100),
#'   breaks = seq(-2, 2, length.out = 101))
#'
#' # Single-hue ramp matching a group's categorical colour
#' pheatmap::pheatmap(mat, color = get_seq_colors("teal", 100))
#' }
NULL

#' @rdname bioinfo_colors
#' @export
#'
#' @return A character vector of hexadecimal color codes.
#'
#' @examples
#' \dontrun{
#' head(get_cat_colors(3))
#' }
get_cat_colors <- function(n) {
  if (!is.numeric(n) || length(n) != 1L || is.na(n)) {
    stop("'n' must be a single number, got: ", paste(class(n), collapse = "/"))
  }
  if (n <= 0) return(character(0))
  if (n != round(n)) stop("'n' must be a whole number, got: ", n)
  if (n <= 10L) {
    bioinfo_cat10[seq_len(n)]
  } else if (n <= 20L) {
    unname(bioinfo_cat20[seq_len(n)])
  } else {
    grDevices::colorRampPalette(bioinfo_cat20, space = "Lab")(n)
  }
}

#' @rdname bioinfo_colors
#' @export
#'
#' @return A character vector of hexadecimal color codes.
#'
#' @examples
#' \dontrun{
#' head(get_seq_colors("teal", 5))
#' }
get_seq_colors <- function(name = "blue", n = 100) {
  choices <- names(bioinfo_seq_anchors)
  if (!name %in% choices) {
    stop("'name' must be one of: ", paste(choices, collapse = ", "), call. = FALSE)
  }
  bioinfo_check_n(n)
  grDevices::colorRampPalette(bioinfo_seq_anchors[[name]], space = "Lab")(n)
}

#' @rdname bioinfo_colors
#' @export
#'
#' @return A character vector of hexadecimal color codes.
#'
#' @examples
#' \dontrun{
#' head(get_div_colors(5))
#' }
get_div_colors <- function(n = 100, palette = c("soft", "strong")) {
  palette <- match.arg(palette)
  bioinfo_check_n(n)
  grDevices::colorRampPalette(bioinfo_div_anchors[[palette]], space = "Lab")(n)
}

#' @rdname bioinfo_colors
#' @export
#'
#' @return A character vector of hexadecimal color codes.
#'
#' @examples
#' \dontrun{
#' head(get_multiseq_colors(4))
#' }
get_multiseq_colors <- function(n = 100) {
  bioinfo_check_n(n)
  grDevices::colorRampPalette(bioinfo_multiseq_anchors, space = "Lab")(n)
}

# Discrete categorical scales =================================================

#' Discrete bioinfo colour scale
#'
#' Maps discrete groups onto the bioinfo categorical palette via
#' [get_cat_colors()]: dark tones first for up to 10 groups, light tones
#' interleaved up to 20, interpolated beyond.
#'
#' @param alpha Transparency in (0, 1] applied to every colour.
#' @param ... Passed to [ggplot2::discrete_scale()].
#'
#' @return A ggplot2 discrete scale object.
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export scale_color_bioinfo
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#'
#' ggplot(iris, aes(Sepal.Length, Sepal.Width, colour = Species)) +
#'   geom_point(size = 2) +
#'   scale_color_bioinfo()
#' }
scale_color_bioinfo <- function(alpha = 1, ...) {
  bioinfo_check_alpha(alpha)
  bioinfo_discrete_scale("colour", alpha, ...)
}

#' @rdname scale_color_bioinfo
#' @export scale_colour_bioinfo
#'
#' @return A ggplot2 Scale object (discrete) that can be added to a ggplot with \code{+}.
#'
#' @examples
#' \dontrun{
#' df_bioinfo <- data.frame(
#'   x = rep(1:3, times = 3),
#'   y = stats::rnorm(9),
#'   group = rep(c("A", "B", "C"), each = 3)
#' )
#' ggplot2::ggplot(df_bioinfo, ggplot2::aes(x = x, y = y, colour = group)) +
#'   ggplot2::geom_point(size = 3) +
#'   scale_colour_bioinfo()
#' }
scale_colour_bioinfo <- scale_color_bioinfo

#' @rdname scale_color_bioinfo
#' @export scale_fill_bioinfo
#'
#' @return A ggplot2 Scale object (discrete) that can be added to a ggplot with \code{+}.
#'
#' @examples
#' \dontrun{
#' df_bioinfo <- data.frame(
#'   x = rep(1:3, times = 3),
#'   y = stats::rnorm(9),
#'   group = rep(c("A", "B", "C"), each = 3)
#' )
#' ggplot2::ggplot(df_bioinfo, ggplot2::aes(x = x, y = y, fill = group)) +
#'   ggplot2::geom_col() +
#'   scale_fill_bioinfo()
#' }
scale_fill_bioinfo <- function(alpha = 1, ...) {
  bioinfo_check_alpha(alpha)
  bioinfo_discrete_scale("fill", alpha, ...)
}

# Single-hue sequential scales =================================================

#' Single-hue sequential bioinfo scales
#'
#' Continuous gradient in one of the ten bioinfo hue families
#' ([get_seq_colors()]); pairing it with the matching categorical colour
#' keeps a group's hue consistent across plot layers.
#'
#' @param name Hue family: "blue", "red", "teal", "purple", "orange",
#'   "pink", "gold", "navy", "olive", or "brown".
#' @param alpha Transparency in (0, 1] applied to every colour.
#' @param reverse Logical; reverse the ramp (default FALSE).
#' @param ... Passed to [ggplot2::scale_color_gradientn()].
#'
#' @return A ggplot2 continuous scale object.
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export scale_color_bioinfo_seq
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#'
#' ggplot(faithfuld, aes(waiting, eruptions, fill = density)) +
#'   geom_tile() +
#'   scale_fill_bioinfo_seq("blue")
#' }
scale_color_bioinfo_seq <- function(name = "blue", alpha = 1, reverse = FALSE, ...) {
  bioinfo_check_alpha(alpha)
  cols <- bioinfo_style_cols(
    get_seq_colors(name, 256L), alpha = alpha, reverse = reverse
  )
  ggplot2::scale_color_gradientn(colours = cols, ...)
}

#' @rdname scale_color_bioinfo_seq
#' @export scale_colour_bioinfo_seq
#'
#' @return A ggplot2 Scale object (continuous) for gradient color/fill mappings.
#'
#' @examples
#' \dontrun{
#' df_seq <- data.frame(x = 1:20, y = stats::rnorm(20))
#' ggplot2::ggplot(df_seq, ggplot2::aes(x = x, y = y, colour = y)) +
#'   ggplot2::geom_point(size = 3) +
#'   scale_colour_bioinfo_seq("teal")
#' }
scale_colour_bioinfo_seq <- scale_color_bioinfo_seq

#' @rdname scale_color_bioinfo_seq
#' @export scale_fill_bioinfo_seq
#'
#' @return A ggplot2 Scale object (continuous) for gradient color/fill mappings.
#'
#' @examples
#' \dontrun{
#' df_seq <- expand.grid(x = 1:5, y = 1:5)
#' df_seq$z <- stats::rnorm(25)
#' ggplot2::ggplot(df_seq, ggplot2::aes(x, y, fill = z)) +
#'   ggplot2::geom_tile() +
#'   scale_fill_bioinfo_seq("blue")
#' }
scale_fill_bioinfo_seq <- function(name = "blue", alpha = 1, reverse = FALSE, ...) {
  bioinfo_check_alpha(alpha)
  cols <- bioinfo_style_cols(
    get_seq_colors(name, 256L), alpha = alpha, reverse = reverse
  )
  ggplot2::scale_fill_gradientn(colours = cols, ...)
}

# Diverging scales =============================================================

#' Diverging bioinfo scales
#'
#' Continuous gradient through a near-white midpoint ([get_div_colors()]),
#' suited to data centred on zero: z-scores, log2 fold changes, correlations.
#' To pin the midpoint to a non-zero value, pass `values` as described in
#' [ggplot2::scale_fill_gradientn()]. For a fixed pheatmap-style scheme see
#' [scale_fill_heatmap()].
#'
#' @param palette Diverging variant: "soft" (five-step) or "strong"
#'   (saturated three-step).
#' @param alpha Transparency in (0, 1] applied to every colour.
#' @param reverse Logical; reverse the ramp (default FALSE).
#' @param ... Passed to [ggplot2::scale_color_gradientn()].
#'
#' @return A ggplot2 continuous scale object.
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export scale_color_bioinfo_div
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#'
#' df <- expand.grid(x = 1:5, y = 1:5)
#' df$z <- rnorm(25)
#'
#' ggplot(df, aes(x, y, fill = z)) +
#'   geom_tile() +
#'   scale_fill_bioinfo_div()
#'
#' # Centre the midpoint on 0 explicitly
#' ggplot(df, aes(x, y, fill = z)) +
#'   geom_tile() +
#'   scale_fill_bioinfo_div(values = scales::rescale(c(min(df$z), 0, max(df$z))))
#' }
scale_color_bioinfo_div <- function(palette = c("soft", "strong"), alpha = 1,
                                    reverse = FALSE, ...) {
  palette <- match.arg(palette)
  bioinfo_check_alpha(alpha)
  cols <- bioinfo_style_cols(
    get_div_colors(256L, palette), alpha = alpha, reverse = reverse
  )
  ggplot2::scale_color_gradientn(colours = cols, ...)
}

#' @rdname scale_color_bioinfo_div
#' @export scale_colour_bioinfo_div
#'
#' @return A ggplot2 Scale object (continuous) for gradient color/fill mappings.
#'
#' @examples
#' \dontrun{
#' df_div <- data.frame(x = 1:20, y = stats::rnorm(20))
#' ggplot2::ggplot(df_div, ggplot2::aes(x = x, y = y, colour = y)) +
#'   ggplot2::geom_point(size = 3) +
#'   scale_colour_bioinfo_div()
#' }
scale_colour_bioinfo_div <- scale_color_bioinfo_div

#' @rdname scale_color_bioinfo_div
#' @export scale_fill_bioinfo_div
#'
#' @return A ggplot2 Scale object (continuous) for gradient color/fill mappings.
#'
#' @examples
#' \dontrun{
#' df_div <- expand.grid(x = 1:5, y = 1:5)
#' df_div$z <- stats::rnorm(25)
#' ggplot2::ggplot(df_div, ggplot2::aes(x, y, fill = z)) +
#'   ggplot2::geom_tile() +
#'   scale_fill_bioinfo_div()
#' }
scale_fill_bioinfo_div <- function(palette = c("soft", "strong"), alpha = 1,
                                   reverse = FALSE, ...) {
  palette <- match.arg(palette)
  bioinfo_check_alpha(alpha)
  cols <- bioinfo_style_cols(
    get_div_colors(256L, palette), alpha = alpha, reverse = reverse
  )
  ggplot2::scale_fill_gradientn(colours = cols, ...)
}

# Palette preview ==============================================================

#' Preview bioinfo palettes
#'
#' Draws every palette in the bioinfo family with [scales::show_col()] on the
#' current graphics device — handy for picking a hue family before plotting.
#'
#' @param family Which palettes to draw: "all" (default), "cat", "seq"
#'   (includes the multi-hue ramp), or "div".
#'
#' @return `NULL`, invisibly; called for its side effect of drawing.
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' \dontrun{
#' preview_bioinfo_palettes()
#' preview_bioinfo_palettes(family = "seq")
#' }
preview_bioinfo_palettes <- function(family = c("all", "cat", "seq", "div")) {
  family <- match.arg(family)

  rows <- list()
  if (family %in% c("all", "cat")) {
    rows[["categorical (20)"]] <- get_cat_colors(20)
  }
  if (family %in% c("all", "seq")) {
    for (nm in names(bioinfo_seq_anchors)) {
      rows[[paste0("seq: ", nm)]] <- get_seq_colors(nm, 10)
    }
    rows[["seq: multi-hue"]] <- get_multiseq_colors(10)
  }
  if (family %in% c("all", "div")) {
    rows[["div: soft"]] <- get_div_colors(11)
    rows[["div: strong"]] <- get_div_colors(11, "strong")
  }

  old_par <- graphics::par(mfrow = c(length(rows), 1L), mar = c(1.5, 4, 1.5, 1))
  on.exit(graphics::par(old_par), add = TRUE)

  for (i in seq_along(rows)) {
    scales::show_col(rows[[i]])
    graphics::mtext(names(rows)[i], side = 2, las = 1, line = 0.5, cex = 0.7)
  }
  invisible(NULL)
}
