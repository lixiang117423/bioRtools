#' Plot a Single-Species Idiogram
#'
#' Draw chromosome idiograms (vertical, bottom-aligned) with optional overlaid
#' heatmaps and annotation tracks. A ggplot2 reimplementation of the
#' single-species mode of RIdeogram's \code{ideogram()}; returns a ggplot
#' object instead of writing an SVG file.
#'
#' @param karyotype Data frame. Either 3 columns \code{Chr, Start, End}, or
#'   5 columns \code{Chr, Start, End, CE_start, CE_end} (centromere; drawn as
#'   a constriction).
#' @param overlaid Optional data frame with columns \code{Chr, Start, End,
#'   Value}, drawn as a heatmap filling each chromosome body (e.g. gene
#'   density).
#' @param label Optional data frame. Required columns depend on
#'   \code{label_type}:
#'   \describe{
#'     \item{marker}{\code{Type, Shape, Chr, Start, End, color}.
#'       \code{Shape} is one of \code{triangle}, \code{box}, \code{circle};
#'       \code{color} is a hex code without \code{#}. Shapes are placed to the
#'       right of each chromosome with leader lines, vertically repelled to
#'       avoid overlap.}
#'     \item{heatmap}{\code{Chr, Start, End, Value}. A second heatmap track
#'       drawn to the right of the chromosomes (uses \code{color_label}).}
#'     \item{line}{\code{Chr, Start, End, Value, color} (one line) or
#'       \code{Chr, Start, End, Value_1, Value_2, color_1, color_2} (two).}
#'     \item{polygon}{Same columns as \code{line}, drawn as a filled area.}
#'   }
#' @param label_type One of \code{"marker"}, \code{"heatmap"}, \code{"line"},
#'   \code{"polygon"}. Ignored (and \code{label} unused) when \code{label}
#'   is \code{NULL}.
#' @param color_overlaid Colors for the \code{overlaid} heatmap (passed to
#'   \code{scale_fill_gradientn}).
#' @param color_label Colors for the \code{label_type = "heatmap"} track.
#' @param width Overall plot width in mm (default 170). Larger values widen
#'   every chromosome proportionally.
#'
#' @return A \link[ggplot2]{ggplot} object.
#'
#' @seealso \code{\link{plot_ideogram_synteny}}
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' \dontrun{
#' data(df.ideo.karyotype)
#' data(df.ideo.gene_density)
#' data(df.ideo.rna_marker)
#'
#' plot_ideogram(df.ideo.karyotype)
#' plot_ideogram(df.ideo.karyotype, overlaid = df.ideo.gene_density)
#' plot_ideogram(df.ideo.karyotype, overlaid = df.ideo.gene_density,
#'              label = df.ideo.rna_marker, label_type = "marker")
#' }
plot_ideogram <- function(karyotype, overlaid = NULL, label = NULL,
                          label_type = NULL,
                          color_overlaid = c("#4575b4", "#ffffbf", "#d73027"),
                          color_label = c("#b35806", "#f7f7f7", "#542788"),
                          width = 170) {

  if (!is.data.frame(karyotype) || nrow(karyotype) == 0L) {
    stop("'karyotype' must be a non-empty data frame")
  }
  need <- c("Chr", "Start", "End")
  if (!all(need %in% names(karyotype))) {
    stop("'karyotype' must have columns Chr, Start, End ",
         "(optionally CE_start, CE_end)")
  }
  has_centromere <- all(c("CE_start", "CE_end") %in% names(karyotype))

  layout <- ideo_layout_chr(karyotype, width)
  max_end <- max(karyotype$End)
  cw <- layout$chr_width[1]
  layout_cols <- c("Chr", "chr_end", "x_left", "x_right", "y_top", "y_bottom", "chr_width")

  p <- ggplot2::ggplot()

  # 1. chromosome bodies (white background; caps/centromere read as white) ----
  body_poly <- do.call(rbind, lapply(seq_len(nrow(layout)), function(k) {
    v <- ideo_chr_outline(layout[k, ])
    data.frame(v, chr_index = layout$chr_index[k])
  }))
  p <- p + ggplot2::geom_polygon(
    data = body_poly,
    ggplot2::aes(x = .data$x, y = .data$y, group = .data$chr_index),
    fill = "white", color = NA)

  # 2. overlaid heatmap -------------------------------------------------------
  if (!is.null(overlaid)) {
    od <- merge_validate(overlaid, layout[, layout_cols], "overlaid", by = "Chr")
    od$ymin <- ideo_y_of_pos(od$Start, od$chr_end, max_end)
    od$ymax <- ideo_y_of_pos(od$End,   od$chr_end, max_end)
    od$ymin <- pmax(od$ymin, od$y_top)      # keep rounded caps white
    od$ymax <- pmin(od$ymax, od$y_bottom)
    p <- p + ggplot2::geom_rect(
      data = od,
      ggplot2::aes(xmin = .data$x_left, xmax = .data$x_right,
                   ymin = .data$ymin, ymax = .data$ymax, fill = .data$Value))
    p <- p + ggplot2::scale_fill_gradientn(
      colours = color_overlaid,
      guide = ggplot2::guide_colorbar(title = "Low → High", title.position = "top",
                                      frame.colour = "black", ticks.colour = "black"))
  }

  # 3. centromere white band (over heatmap, under outline) -------------------
  if (has_centromere) {
    ce_bands <- do.call(rbind, lapply(seq_len(nrow(layout)), function(k) {
      v <- ideo_centromere_band(layout[k, ])
      if (is.null(v)) return(NULL)
      data.frame(v, chr_index = layout$chr_index[k])
    }))
    if (!is.null(ce_bands)) {
      p <- p + ggplot2::geom_polygon(
        data = ce_bands,
        ggplot2::aes(x = .data$x, y = .data$y, group = .data$chr_index),
        fill = "white", color = NA)
    }
  }

  # 4. label track ------------------------------------------------------------
  if (!is.null(label)) {
    label_type <- match.arg(label_type, c("marker", "heatmap", "line", "polygon"))
    p <- switch(label_type,
      marker  = add_label_marker(p, label, layout, max_end, cw,
                                 has_overlaid = !is.null(overlaid)),
      heatmap = add_label_heatmap(p, label, layout, max_end, cw, color_label,
                                  has_overlaid = !is.null(overlaid)),
      line    = add_label_line(p, label, layout, max_end, cw, polygon = FALSE),
      polygon = add_label_line(p, label, layout, max_end, cw, polygon = TRUE))
  }

  # 5. chromosome outlines (on top) ------------------------------------------
  p <- p + ggplot2::geom_polygon(
    data = body_poly,
    ggplot2::aes(x = .data$x, y = .data$y, group = .data$chr_index),
    fill = NA, color = "grey50", linewidth = 0.4)

  # 6. chromosome names (below each chromosome) ------------------------------
  chr_text <- data.frame(
    x = (layout$x_left + layout$x_right) / 2 - nchar(as.character(layout$Chr)) * 2.2,
    y = (ideo_max_chr + 25) * ideo_mpx + 15,
    label = as.character(layout$Chr)
  )
  p <- p + ggplot2::geom_text(
    data = chr_text, ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
    size = 9 * 0.3528, color = "black")  # 9 pt in mm

  ideo_finish(p, layout)
}

# --- merge with a friendly error if a Chr has no layout ---------------------
merge_validate <- function(df, layout_cols_df, arg_name, by) {
  if (!by %in% names(df)) stop("'", arg_name, "' must have a '", by, "' column")
  matched <- intersect(unique(df[[by]]), unique(layout_cols_df[[by]]))
  dropped <- setdiff(unique(df[[by]]), unique(layout_cols_df[[by]]))
  if (length(dropped)) {
    warning("'", arg_name, "' contains chromosome(s) not in 'karyotype$Chr': ",
            paste(head(dropped, 10), collapse = ", "),
            if (length(dropped) > 10) " ...",
            " (these rows are dropped)", call. = FALSE)
  }
  merged <- merge(df, layout_cols_df, by = by, sort = FALSE)
  if (nrow(merged) == 0L) {
    stop("No chromosome in '", arg_name, "' matches 'karyotype$Chr'")
  }
  merged
}

# --- marker track -----------------------------------------------------------
add_label_marker <- function(p, label, layout, max_end, cw, has_overlaid = FALSE) {
  if (has_overlaid) p <- p + ggnewscale::new_scale_fill()
  ld <- merge_validate(label, layout[, c("Chr", "chr_end", "x_right", "x_left", "chr_width")],
                       "label", by = "Chr")
  ld$mid <- (ld$Start + ld$End) / 2
  ld$y0  <- ideo_y_of_pos(ld$mid, ld$chr_end, max_end)
  ld$x_marker <- ld$x_right + ld$chr_width / 2

  # per-chromosome vertical repel (sort by y0 within chr, push apart)
  ld <- ld[order(ld$Chr, ld$y0), ]
  ld$y <- ld$y0
  for (ch in unique(ld$Chr)) {
    sel <- ld$Chr == ch
    if (sum(sel) > 1L) {
      ld$y[sel] <- ideo_repel_1d(ld$y0[sel], cw / 3)
    }
  }

  shape_map  <- c(triangle = 24, box = 22, circle = 21)
  ld$shape_pch <- shape_map[ld$Shape]
  type_color <- unique(ld[, c("Type", "color")])
  type_color$fill <- paste0("#", type_color$color)
  size_mm <- (cw / 2) / ideo_mpx

  # leader lines
  p <- p + ggplot2::geom_segment(
    data = ld,
    ggplot2::aes(x = .data$x_right, y = .data$y0,
                 xend = .data$x_marker, yend = .data$y),
    color = paste0("#", ld$color), linewidth = 0.25)
  # markers
  p <- p + ggplot2::geom_point(
    data = ld,
    ggplot2::aes(x = .data$x_marker, y = .data$y, fill = .data$Type,
                 shape = .data$Shape),
    size = size_mm, color = "black", stroke = 0.3) +
    ggplot2::scale_shape_manual(values = shape_map) +
    ggplot2::scale_fill_manual(values = stats::setNames(type_color$fill, type_color$Type))
  p
}

# --- second heatmap track (to the right of the chromosomes) ----------------
add_label_heatmap <- function(p, label, layout, max_end, cw, color_label, has_overlaid) {
  if (has_overlaid) p <- p + ggnewscale::new_scale_fill()
  ld <- merge_validate(label, layout[, c("Chr", "chr_end", "x_left", "x_right",
                                         "y_top", "y_bottom", "chr_width")],
                       "label", by = "Chr")
  offset <- 1.2 * cw
  ld$ymin <- pmax(ideo_y_of_pos(ld$Start, ld$chr_end, max_end), ld$y_top)
  ld$ymax <- pmin(ideo_y_of_pos(ld$End,   ld$chr_end, max_end), ld$y_bottom)
  p <- p + ggplot2::geom_rect(
    data = ld,
    ggplot2::aes(xmin = .data$x_left + offset, xmax = .data$x_right + offset,
                 ymin = .data$ymin, ymax = .data$ymax, fill = .data$Value)) +
    ggplot2::scale_fill_gradientn(
      colours = color_label,
      guide = ggplot2::guide_colorbar(title = "Low → High", title.position = "top",
                                      frame.colour = "black", ticks.colour = "black"))

  # second chromosome outline, offset right
  poly2 <- do.call(rbind, lapply(seq_len(nrow(layout)), function(k) {
    v <- ideo_chr_outline(layout[k, ])
    data.frame(v, chr_index = layout$chr_index[k])
  }))
  poly2$x <- poly2$x + offset
  p <- p + ggplot2::geom_polygon(
    data = poly2,
    ggplot2::aes(x = .data$x, y = .data$y, group = .data$chr_index),
    fill = NA, color = "grey50", linewidth = 0.4)
  p
}

# --- line / polygon track ---------------------------------------------------
add_label_line <- function(p, label, layout, max_end, cw, polygon = FALSE) {
  two_line <- all(c("Value_1", "Value_2") %in% names(label))
  if (two_line) {
    series <- list(c("Value_1", "color_1"), c("Value_2", "color_2"))
  } else {
    series <- list(c("Value", "color"))
  }
  offset <- 1.2 * cw
  for (s in series) {
    ld <- merge_validate(label, layout[, c("Chr", "chr_end", "x_left", "chr_width")],
                         "label", by = "Chr")
    ld$val <- ld[[s[1]]]
    ld$x <- ld$x_left + offset + ld$val * cw / max(ld$val, na.rm = TRUE)
    ld$y <- ideo_y_of_pos((ld$Start + ld$End) / 2, ld$chr_end, max_end)
    col_hex <- paste0("#", label[[s[2]]][1])
    ld <- ld[order(ld$Chr, ld$Start), ]
    if (polygon) {
      p <- p + ggplot2::geom_polygon(
        data = ld, ggplot2::aes(x = .data$x, y = .data$y, group = .data$Chr),
        fill = col_hex, color = col_hex, linewidth = 0.25, alpha = 0.6)
    } else {
      p <- p + ggplot2::geom_path(
        data = ld, ggplot2::aes(x = .data$x, y = .data$y, group = .data$Chr),
        color = col_hex, linewidth = 0.4)
    }
  }
  p
}

# --- common finalization lives in utils_ideogram.R (ideo_finish) ------------
