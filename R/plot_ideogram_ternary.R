#' Plot a Three-Species (Ternary) Synteny Idiogram
#'
#' Draw three species' chromosomes in a triangular arrangement — species 1
#' horizontal along the bottom, species 2 rotated -60 degrees (upper-left),
#' species 3 rotated +60 degrees (upper-right) — with ribbons connecting
#' syntenic blocks for each of the three species pairs. A ggplot2
#' reimplementation of the ternary mode of RIdeogram's \code{ideogram()}.
#'
#' The chromosome layout matches the original RIdeogram coordinates exactly.
#' Ribbons use a clean self-contained geometry (anchors land on the rotated
#' chromosome edges, connected by smooth curves) rather than the original's
#' hand-tuned artistic offsets, so ribbon curvature is close but not
#' pixel-identical to the original.
#'
#' @param karyotype Data frame with columns \code{Chr, Start, End, fill,
#'   species, size, color}: exactly three species in display order (first
#'   species on the bottom, second upper-left, third upper-right).
#' @param synteny Data frame in RIdeogram's ternary format: columns
#'   \code{Species_1, Start_2, End_2, Species_2, Start_1, End_1, fill, type}.
#'   \code{Species_1}/\code{Species_2} are 1-based chromosome indices; the
#'   pair of species a row connects is given by \code{type} (1 = species 1-2,
#'   2 = species 1-3, 3 = species 2-3). \code{fill} is a hex colour without
#'   \code{#}, or the literal \code{"gradient"} to blend the two species'
#'   colours along the ribbon.
#' @param width Overall layout width in mm (default 170).
#'
#' @return A \link[ggplot2]{ggplot} object.
#'
#' @seealso \code{\link{plot_ideogram}}, \code{\link{plot_ideogram_synteny}}
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' \dontrun{
#' data(df.ideo.ternary_karyotype)
#' data(df.ideo.ternary)
#' plot_ideogram_ternary(df.ideo.ternary_karyotype, df.ideo.ternary)
#' }
plot_ideogram_ternary <- function(karyotype, synteny, width = 170) {

  if (!is.data.frame(karyotype) || nrow(karyotype) == 0L) {
    stop("'karyotype' must be a non-empty data frame")
  }
  need <- c("Chr", "Start", "End", "fill", "species", "size", "color")
  if (!all(need %in% names(karyotype))) {
    stop("'karyotype' must have columns Chr, Start, End, fill, species, size, color")
  }
  species <- unique(karyotype$species)
  if (length(species) != 3L) {
    stop("plot_ideogram_ternary() supports exactly three species; found ", length(species))
  }
  if (!is.data.frame(synteny) || nrow(synteny) == 0L) {
    stop("'synteny' must be a non-empty data frame")
  }

  kats <- lapply(species, function(s) karyotype[karyotype$species == s, ])
  layouts <- lapply(kats, ideo_ternary_species)
  params <- ideo_ternary_params()

  p <- ggplot2::ggplot()

  # ribbons first so chromosomes sit on top of their edges -------------------
  ribbons <- build_ternary_ribbons(synteny, layouts, params, kats)
  for (hex in unique(ribbons$fill_hex)) {
    p <- p + ggplot2::geom_polygon(
      data = ribbons[ribbons$fill_hex == hex, ],
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$ribbon_id),
      fill = hex, colour = hex, linewidth = 0.05)
  }

  # chromosomes: transparent outline + coloured inner band -------------------
  for (sp in 1:3) {
    chr_poly <- build_ternary_chr_poly(kats[[sp]], layouts[[sp]], sp, params)
    p <- p + ggplot2::geom_polygon(
      data = chr_poly$outer,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$group),
      fill = NA, color = "black", linewidth = 0.2)
    p <- p + ggplot2::geom_polygon(
      data = chr_poly$inner,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$group),
      fill = chr_poly$inner_fill, color = chr_poly$inner_fill, linewidth = 0.1)
  }

  # species-name labels at the triangle corners ------------------------------
  sp_lab <- build_ternary_species_labels(species, kats, layouts, params)
  for (i in seq_len(nrow(sp_lab))) {
    p <- p + ggplot2::geom_text(
      data = sp_lab[i, ],
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      size = sp_lab$size_mm[i], color = sp_lab$color[i], angle = sp_lab$angle[i])
  }

  ideo_finish(p, NULL)
}

# Chromosome polygons for one species: outer outline (h=8) and inner band
# (h=4, centred), rotated per species.
build_ternary_chr_poly <- function(kat, lay, sp_idx, params) {
  bl <- params$baseline
  rot <- switch(as.character(sp_idx),
    "2" = function(m) ideo_rotate(m, -60, params$pivot2),
    "3" = function(m) ideo_rotate(m,  60, params$pivot3),
    function(m) m)

  outer_list <- inner_list <- vector("list", nrow(kat))
  for (i in seq_len(nrow(kat))) {
    x <- lay$starts[i]; w <- lay$widths[i]
    outer_local <- matrix(c(x, bl, x + w, bl, x + w, bl + 8, x, bl + 8), nrow = 4, byrow = TRUE)
    inner_local <- matrix(c(x, bl + 2, x + w, bl + 2, x + w, bl + 6, x, bl + 6), nrow = 4, byrow = TRUE)
    og <- rot(outer_local); ig <- rot(inner_local)
    outer_list[[i]] <- data.frame(x = og[, 1], y = og[, 2], group = i)
    inner_list[[i]] <- data.frame(x = ig[, 1], y = ig[, 2], group = i)
  }
  list(outer = do.call(rbind, outer_list),
       inner = do.call(rbind, inner_list),
       inner_fill = paste0("#", kat$fill[1]))
}

# All ternary ribbons. Solid ribbons become one polygon; "gradient" ribbons
# are subdivided into thin slices interpolating the two species' colours.
build_ternary_ribbons <- function(synteny, layouts, params, kats) {
  pair <- list(`1` = c(1, 2), `2` = c(1, 3), `3` = c(2, 3))  # type -> (sp_a, sp_b)
  n <- nrow(synteny)
  out <- vector("list", n)
  for (i in seq_len(n)) {
    sp <- pair[[as.character(synteny$type[i])]]
    sp_a <- sp[1]; sp_b <- sp[2]
    a1 <- ideo_ternary_anchor(sp_a, synteny$Species_1[i], synteny$Start_2[i], layouts, params)
    a2 <- ideo_ternary_anchor(sp_a, synteny$Species_1[i], synteny$End_2[i],   layouts, params)
    b1 <- ideo_ternary_anchor(sp_b, synteny$Species_2[i], synteny$Start_1[i], layouts, params)
    b2 <- ideo_ternary_anchor(sp_b, synteny$Species_2[i], synteny$End_1[i],   layouts, params)
    fill <- synteny$fill[i]
    out[[i]] <- if (fill == "gradient") {
      gradient_ribbon(a1, a2, b1, b2,
                      paste0("#", kats[[sp_a]]$fill[1]),
                      paste0("#", kats[[sp_b]]$fill[1]), i)
    } else {
      quad_ribbon(a1, a2, b1, b2, paste0("#", fill), i)
    }
  }
  do.call(rbind, out)
}

# A single solid-colour quad ribbon with slightly bowed cross-edges.
quad_ribbon <- function(a1, a2, b1, b2, hex, id, n_curve = 12) {
  edge1 <- ideo_qt_edge(a1, a1, (a1 + b1) / 2, b1, n_curve)
  edge2 <- ideo_qt_edge(a2, a2, (a2 + b2) / 2, b2, n_curve)
  poly <- rbind(edge1, edge2[nrow(edge2):1, ])
  data.frame(x = poly[, 1], y = poly[, 2], ribbon_id = id, fill_hex = hex)
}

# A gradient ribbon: N slices interpolated between the two species' colours.
gradient_ribbon <- function(a1, a2, b1, b2, hex_a, hex_b, id, n_slice = 12) {
  ca <- grDevices::col2rgb(hex_a); cb <- grDevices::col2rgb(hex_b)
  ramp <- colorRampSubdiv(ca, cb, n_slice)
  pts_top <- (seq_len(n_slice + 1) - 1) / n_slice   # a-side interpolation fractions
  rows <- vector("list", n_slice)
  for (k in seq_len(n_slice)) {
    t0 <- pts_top[k]; t1 <- pts_top[k + 1]
    p1 <- a1 + (b1 - a1) * t0; p2 <- a1 + (b1 - a1) * t1
    p3 <- a2 + (b2 - a2) * t1; p4 <- a2 + (b2 - a2) * t0
    poly <- rbind(p1, p2, p3, p4)
    rows[[k]] <- data.frame(x = poly[, 1], y = poly[, 2],
                            ribbon_id = id, fill_hex = ramp[k])
  }
  do.call(rbind, rows)
}

# n colours interpolating from rgb matrix ca to cb (grDevices::colorRampPalette
# would reorder; do it manually to keep direction).
colorRampSubdiv <- function(ca, cb, n) {
  cols <- vector("character", n)
  for (k in seq_len(n)) {
    t <- (k - 1) / (n - 1)
    rgb <- round(ca + (cb - ca) * t)
    cols[k] <- grDevices::rgb(rgb[1, 1], rgb[2, 1], rgb[3, 1], maxColorValue = 255)
  }
  cols
}

# Species-name labels, placed near each species' chromosome row, rotated to
# match. size is the original pt size -> mm.
build_ternary_species_labels <- function(species, kats, layouts, params) {
  bl <- params$baseline
  lab <- data.frame(label = species, size_pt = sapply(kats, function(k) k$size[1]),
                    color = paste0("#", sapply(kats, function(k) k$color[1])))
  mid <- function(lay) (lay$starts[1] + tail(lay$starts, 1) + tail(lay$widths, 1)) / 2
  pts <- list(
    c(mid(layouts[[1]]), bl + 8 + 0.8 * ideo_cm_px),                                  # species 1 below
    ideo_rotate(cbind(mid(layouts[[2]]), bl + 8 + 0.4 * ideo_cm_px), -60, params$pivot2),
    ideo_rotate(cbind(mid(layouts[[3]]), bl + 8 + 0.4 * ideo_cm_px),  60, params$pivot3))
  lab$x <- sapply(pts, function(p) p[1])
  lab$y <- sapply(pts, function(p) p[2])
  lab$angle <- c(0, -60, 60)
  lab$size_mm <- lab$size_pt * 0.3528
  lab
}
