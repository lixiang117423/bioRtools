#' Plot a Dual-Species Synteny Idiogram
#'
#' Draw two species' chromosomes as horizontal rows with quadratic-Bezier
#' ribbons connecting syntenic regions. A ggplot2 reimplementation of the
#' dual-species mode of RIdeogram's \code{ideogram()}; returns a ggplot
#' object instead of writing an SVG file.
#'
#' @param karyotype Data frame with columns \code{Chr, Start, End, fill,
#'   species, size, color}: exactly two species, in the order they should
#'   appear (first species on top, second on the bottom). \code{fill} is the
#'   chromosome body colour (hex without \code{#}); \code{size}/\code{color}
#'   control the species-label font size and colour.
#' @param synteny Data frame with columns \code{Species_1, Start_1, End_1,
#'   Species_2, Start_2, End_2, fill}. \code{Species_1}/\code{Species_2} are
#'   1-based indices into the corresponding species' chromosome block;
#'   \code{fill} is the ribbon colour (hex without \code{#}).
#' @param width Overall layout width in mm (default 170). The horizontal
#'   extent of each species' chromosome row scales with this value.
#'
#' @return A \link[ggplot2]{ggplot} object.
#'
#' @seealso \code{\link{plot_ideogram}} for the single-species idiogram;
#'   \code{\link{plot_synteny}} for an unrelated gene-level synteny plot.
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' \dontrun{
#' data(df.ideo.synteny_karyotype)
#' data(df.ideo.synteny)
#' plot_ideogram_synteny(df.ideo.synteny_karyotype, df.ideo.synteny)
#' }
plot_ideogram_synteny <- function(karyotype, synteny, width = 170) {

  if (!is.data.frame(karyotype) || nrow(karyotype) == 0L) {
    stop("'karyotype' must be a non-empty data frame")
  }
  need <- c("Chr", "Start", "End", "fill", "species", "size", "color")
  if (!all(need %in% names(karyotype))) {
    stop("'karyotype' must have columns Chr, Start, End, fill, species, size, color")
  }
  species <- unique(karyotype$species)
  if (length(species) != 2L) {
    stop("plot_ideogram_synteny() supports exactly two species; found ", length(species))
  }
  if (!is.data.frame(synteny) || nrow(synteny) == 0L) {
    stop("'synteny' must be a non-empty data frame")
  }

  kat1 <- karyotype[karyotype$species == species[1], ]
  kat2 <- karyotype[karyotype$species == species[2], ]

  total_cm <- width / 10
  widths_total <- function(n) (total_cm - 1.5 - (n - 1) * 0.1) * ideo_cm_px
  rl1 <- ideo_synteny_row_layout(setNames(kat1$End, kat1$Chr), widths_total(nrow(kat1)))
  rl2 <- ideo_synteny_row_layout(setNames(kat2$End, kat2$Chr), widths_total(nrow(kat2)))

  # row y positions (outer rect top edges) and ribbon anchors
  y1 <- 2.5 * ideo_cm_px
  y2 <- (2.5 + 0.5 + 5) * ideo_cm_px
  y_top_bot <- (2.5 + 0.5) * ideo_cm_px         # bottom edge of top row
  y_bot_top <- (2.5 + 0.5 + 5) * ideo_cm_px     # top edge of bottom row
  h_outer <- 17.71654
  h_inner <- 11.71654

  p <- ggplot2::ggplot()

  # ribbons first, so chromosomes sit on top of their edges ------------------
  ribbons <- build_synteny_ribbons(synteny, rl1, rl2, y_top_bot, y_bot_top)
  for (hex in unique(ribbons$fill_hex)) {
    p <- p + ggplot2::geom_polygon(
      data = ribbons[ribbons$fill_hex == hex, ],
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$ribbon_id),
      fill = hex, colour = hex, linewidth = 0.1)
  }

  # chromosomes: outer (light grey) + inner (coloured) -----------------------
  for (row in list(
    list(rl = rl1, kat = kat1, y = y1),
    list(rl = rl2, kat = kat2, y = y2)
  )) {
    outer_df <- data.frame(
      xmin = row$rl$starts, xmax = row$rl$starts + row$rl$widths,
      ymin = row$y,        ymax = row$y + h_outer)
    inner_fill <- paste0("#", row$kat$fill)
    inner_df <- data.frame(
      xmin = row$rl$starts, xmax = row$rl$starts + row$rl$widths,
      ymin = row$y + 3,     ymax = row$y + 3 + h_inner)
    p <- p + ggplot2::geom_rect(data = outer_df,
        ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                     ymin = .data$ymin, ymax = .data$ymax),
        fill = "#f7f7f7", color = "black", linewidth = 0.5)
    p <- p + ggplot2::geom_rect(data = inner_df,
        ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                     ymin = .data$ymin, ymax = .data$ymax),
        fill = inner_fill, color = inner_fill, linewidth = 0.5)
  }

  # chromosome-name text -----------------------------------------------------
  lab1 <- data.frame(
    x = (rl1$starts + rl1$starts + rl1$widths) / 2 - nchar(as.character(kat1$Chr)) * 2.1,
    y = y1 - 0.1 * ideo_cm_px, label = as.character(kat1$Chr))
  lab2 <- data.frame(
    x = (rl2$starts + rl2$starts + rl2$widths) / 2 - nchar(as.character(kat2$Chr)) * 2.1,
    y = y2 + (0.1 + 0.5 + 0.15) * ideo_cm_px, label = as.character(kat2$Chr))
  for (lab in list(lab1, lab2)) {
    p <- p + ggplot2::geom_text(data = lab,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      size = 9 * 0.3528, color = "black")
  }

  # species-name text (right-justified on the left margin, just below each row)
  sp_y <- c(y1, y2) + h_outer - 3
  sp_kat <- list(kat1, kat2)
  for (k in 1:2) {
    p <- p + ggplot2::geom_text(
      data = data.frame(x = 2 * ideo_cm_px, y = sp_y[k], label = species[k]),
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      size = sp_kat[[k]]$size[1] * 0.3528,
      color = paste0("#", sp_kat[[k]]$color[1]), hjust = 1)
  }

  ideo_finish(p, NULL)
}

# Build all ribbon polygons. Each ribbon becomes ~80 rows tagged with a
# ribbon_id and a fill_hex so they can be grouped/coloured efficiently.
build_synteny_ribbons <- function(synteny, rl1, rl2, y1, y2) {
  n <- nrow(synteny)
  out_rows <- vector("list", n)
  for (i in seq_len(n)) {
    sp1 <- synteny$Species_1[i]
    sp2 <- synteny$Species_2[i]
    x1 <- rl1$starts[sp1] + synteny$Start_1[i] * rl1$widths[sp1] / rl1$total_len
    x4 <- rl1$starts[sp1] + synteny$End_1[i]   * rl1$widths[sp1] / rl1$total_len
    x2 <- rl2$starts[sp2] + synteny$Start_2[i] * rl2$widths[sp2] / rl2$total_len
    x3 <- rl2$starts[sp2] + synteny$End_2[i]   * rl2$widths[sp2] / rl2$total_len
    edge_top <- ideo_qt_edge(c(x1, y1), c(x1, y1 + 1.5 * ideo_cm_px),
                             c((x1 + x2) / 2, (y1 + y2) / 2), c(x2, y2))
    edge_bot <- ideo_qt_edge(c(x3, y2), c(x3, y2 - 1.5 * ideo_cm_px),
                             c((x3 + x4) / 2, (y1 + y2) / 2), c(x4, y1))
    poly <- rbind(edge_top, edge_bot)
    out_rows[[i]] <- data.frame(x = poly[, 1], y = poly[, 2],
                                ribbon_id = i,
                                fill_hex = paste0("#", synteny$fill[i]))
  }
  do.call(rbind, out_rows)
}
