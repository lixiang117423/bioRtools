#' Internal geometry helpers for ideogram / synteny plotting
#'
#' Pure coordinate maths (no ggplot2). Shared by \code{\link{plot_ideogram}}
#' and \code{\link{plot_synteny}}. All coordinates follow the original
#' RIdeogram SVG convention: pixels at 90 DPI, y increasing downward.
#'
#' @name utils_ideogram
#' @keywords internal
NULL

# Pixel-per-mm and pixel-per-cm at 90 DPI (90 / 25.4 = 3.543307).
ideo_mpx     <- 3.543307
ideo_cm_px   <- 35.43307
ideo_max_chr <- 150   # the longest chromosome is drawn 150 mm tall

# --- single-species chromosome layout ---------------------------------------

# Compute x/y positions for each chromosome (vertical, bottom-aligned).
# Adds: chr_index, x_left, x_right, y_top, y_bottom, chr_width, and
# y_ce_start / y_ce_end when centromere columns are present.
ideo_layout_chr <- function(karyotype, width) {
  n <- nrow(karyotype)
  chr_width <- width / (2.6 * n) * ideo_mpx
  max_end <- max(karyotype$End)
  i <- seq_len(n)

  out <- data.frame(
    chr_index = i,
    Chr       = karyotype$Chr,
    chr_end   = karyotype$End,
    x_left    = 20 * ideo_mpx + (i - 1) * 2.6 * chr_width,
    chr_width = chr_width,
    stringsAsFactors = FALSE
  )
  out$x_right  <- out$x_left + chr_width
  # body top/bottom already exclude the cap radius (caps stay white)
  out$y_top    <- (25 + ideo_max_chr * (1 - karyotype$End / max_end)) * ideo_mpx + chr_width / 2
  out$y_bottom <- (25 + ideo_max_chr) * ideo_mpx - chr_width / 2

  if (all(c("CE_start", "CE_end") %in% names(karyotype))) {
    out$y_ce_start <- (25 + ideo_max_chr * (1 - (karyotype$End - karyotype$CE_start) / max_end)) * ideo_mpx
    out$y_ce_end   <- (25 + ideo_max_chr * (1 - (karyotype$End - karyotype$CE_end)   / max_end)) * ideo_mpx
  }
  out
}

# Map a genomic position on a chromosome (given its End length) to a y pixel.
ideo_y_of_pos <- function(pos, chr_end, max_end) {
  (25 + ideo_max_chr * (1 - (chr_end - pos) / max_end)) * ideo_mpx
}

# Sample one semicircle cap. side = "top" | "bottom"; returns n x 2 matrix.
ideo_cap_arc <- function(x_left, x_right, y_anchor, side = "top", n = 24) {
  r <- (x_right - x_left) / 2
  cx <- (x_left + x_right) / 2
  if (side == "top") {
    theta <- seq(pi, 2 * pi, length.out = n)   # 180° -> 360° through 270° (apex up)
  } else {
    theta <- seq(0, pi, length.out = n)        # 0° -> 180° through 90° (apex down)
  }
  cbind(cx + r * cos(theta), y_anchor + r * sin(theta))
}

# Build the chromosome outline polygon (with optional centromere pinch).
# Returns a data.frame of (x, y) vertices in path order, plus a group id.
ideo_chr_outline <- function(layout_row, n_arc = 24) {
  xl <- layout_row$x_left
  xr <- layout_row$x_right
  yt <- layout_row$y_top
  yb <- layout_row$y_bottom

  top_cap <- ideo_cap_arc(xl, xr, yt, "top", n_arc)
  bot_cap <- ideo_cap_arc(xl, xr, yb, "bottom", n_arc)

  if (!is.null(layout_row$y_ce_start)) {
    ycs <- layout_row$y_ce_start
    yce <- layout_row$y_ce_end
    verts <- rbind(
      c(xl, yt), top_cap, c(xr, yt),
      c(xr, ycs), c(xl, yce), c(xl, yb),
      bot_cap, c(xr, yb),
      c(xr, yce), c(xl, ycs)
    )
  } else {
    verts <- rbind(
      c(xl, yt), top_cap, c(xr, yt),
      c(xr, yb), bot_cap, c(xl, yb)
    )
  }
  data.frame(x = verts[, 1], y = verts[, 2])
}

# White centromere band polygon (masks the heatmap across the centromere so
# the pinch reads as white), drawn at full chromosome width between CE_start
# and CE_end. Only meaningful for 5-column karyotypes.
ideo_centromere_band <- function(layout_row) {
  if (is.null(layout_row$y_ce_start)) return(NULL)
  data.frame(
    x = c(layout_row$x_left, layout_row$x_right, layout_row$x_right, layout_row$x_left),
    y = c(layout_row$y_ce_start, layout_row$y_ce_start, layout_row$y_ce_end, layout_row$y_ce_end)
  )
}

# --- 1D vertical repel (ported from RIdeogram) ------------------------------

# Push the closest pair of sorted y values apart by `force` until every gap is
# >= force (or the iteration cap is hit). Faithful port of the original.
ideo_repel_1d <- function(y, force, tag = 1L, cap = 3000L) {
  if (length(y) <= 1L) return(y)
  if (min(diff(y)) >= force || tag >= cap) return(y)
  sp <- which.min(diff(y))
  y[sp]     <- y[sp] - force
  y[sp + 1] <- y[sp + 1] + force
  ideo_repel_1d(sort(y), force, tag + 1L, cap)
}

# --- quadratic Bezier sampler (synteny ribbons) -----------------------------

# Quadratic Bezier curve from p0 through control c1 to p2. Returns n x 2.
ideo_qbezier <- function(p0, c1, p2, n = 30) {
  t <- seq(0, 1, length.out = n)
  cbind((1 - t)^2 * p0[1] + 2 * (1 - t) * t * c1[1] + t^2 * p2[1],
        (1 - t)^2 * p0[2] + 2 * (1 - t) * t * c1[2] + t^2 * p2[2])
}

# Sample the synteny ribbon edge: a quadratic-then-smooth-quadratic (Q ... T)
# curve from p_start to p_end, with explicit first control c1 and the second
# control obtained by reflecting c1 about the midpoint (T semantics).
ideo_qt_edge <- function(p_start, c1, p_mid, p_end, n = 40) {
  seg1 <- ideo_qbezier(p_start, c1, p_mid, ceiling(n / 2))
  c2 <- 2 * p_mid - c1                        # T reflects the previous control
  seg2 <- ideo_qbezier(p_mid, c2, p_end, ceiling(n / 2))
  rbind(seg1, seg2[-1, , drop = FALSE])
}

# --- dual-species synteny layout --------------------------------------------

# Compute x position (px) of a genomic position on a horizontal chromosome
# within one species' row. chr_offsets: named vector of left-edge x per Chr.
# widths_total: total px available for the species row; total_len: sum of End.
ideo_synteny_x <- function(chr, pos, chr_offsets, chr_lens, widths_total, total_len) {
  chr_offsets[[as.character(chr)]] +
    pos * widths_total / total_len +
    (match(as.character(chr), names(chr_offsets)) - 1) * 0.1 * ideo_cm_px
}

# Build per-species horizontal chromosome layout: cumulative x offsets and
# per-chr pixel widths, scaled to fit widths_total with 0.1 cm gaps.
ideo_synteny_row_layout <- function(chr_ends, widths_total) {
  chrs <- names(chr_ends)
  total_len <- sum(chr_ends)
  n <- length(chrs)
  widths <- chr_ends * widths_total / total_len
  # 0.1 cm gap between chromosomes; start at 3.5 cm
  gaps <- c(0, rep(0.1 * ideo_cm_px, n - 1))
  starts <- 3.5 * ideo_cm_px + cumsum(c(0, head(widths, -1))) + cumsum(gaps)
  list(starts = starts, widths = widths, total_len = total_len)
}

# Common ggplot finalization for ideogram plots: reverse y (SVG convention),
# fix aspect ratio, strip axes. layout is currently unused but kept for
# potential future limit computation.
ideo_finish <- function(p, layout = NULL) {
  p <- p + ggplot2::scale_y_reverse()
  p <- p + ggplot2::coord_fixed(ratio = 1, clip = "off")
  p <- p + ggplot2::theme_void(base_size = 9)
  p <- p + ggplot2::theme(legend.position = "right")
  p
}

# --- ternary (3-species) layout --------------------------------------------

# Rotate an n x 2 matrix of points around a pivot (degrees, SVG convention).
ideo_rotate <- function(xy, angle_deg, pivot) {
  a <- angle_deg * pi / 180
  dx <- xy[, 1] - pivot[1]
  dy <- xy[, 2] - pivot[2]
  cbind(pivot[1] + dx * cos(a) - dy * sin(a),
        pivot[2] + dx * sin(a) + dy * cos(a))
}

# Fixed layout parameters for the triangular arrangement: species 1 horizontal
# along the bottom, species 2 rotated -60 deg (upper-left), species 3 +60 deg
# (upper-right). All share a baseline; species 2/3 rotate around their pivots.
ideo_ternary_params <- function() {
  baseline <- (2.5 + 11.25833) * ideo_cm_px
  list(
    baseline = baseline,
    pivot2 = c(2 * ideo_cm_px, baseline),
    pivot3 = c(15 * ideo_cm_px, baseline)
  )
}

# Local horizontal layout for one species' chromosomes (mirrors species 1 of
# the original: x from 3 cm, 0.04 cm gaps). Returns per-chr left edges, widths,
# the bp->px scale, and total bp.
ideo_ternary_species <- function(kat) {
  n <- nrow(kat)
  scale <- (11 - (n - 1) * 0.04) * ideo_cm_px
  total <- sum(kat$End)
  widths <- kat$End * scale / total
  starts <- 3 * ideo_cm_px +
    cumsum(c(0, head(widths, -1))) +
    cumsum(c(0, rep(0.04 * ideo_cm_px, n - 1)))
  list(starts = starts, widths = widths, scale = scale, total = total)
}

# Global (x, y) of a genomic position on a species' chromosome. The point is
# placed in the local horizontal frame then rotated per species.
ideo_ternary_anchor <- function(sp_idx, chr_idx, pos, species_layouts, params) {
  lay <- species_layouts[[sp_idx]]
  local_x <- lay$starts[chr_idx] + pos * lay$scale / lay$total
  pt <- cbind(local_x, params$baseline)
  switch(as.character(sp_idx),
    "2" = ideo_rotate(pt, -60, params$pivot2),
    "3" = ideo_rotate(pt,  60, params$pivot3),
    pt)
}

