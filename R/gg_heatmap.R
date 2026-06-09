#' ggplot2-based Heatmap
#'
#' A ggplot2 rewrite of \code{pheatmap::pheatmap()} with a similar interface.
#' Supports hierarchical clustering, dendrograms, row/column annotations,
#' gap lines, and most common heatmap features.
#'
#' @param mat Matrix to plot as a heatmap.
#' @param color Color palette as a character vector.
#'   Default is a 100-color gradient from blue to white to red.
#' @param breaks Numeric vector of breaks for the color scale. If NULL,
#'   automatically generated from the data range.
#' @param scale One of "none", "row", or "column". Whether to scale values
#'   before plotting. Default is "none".
#' @param cluster_rows Logical or hclust object. If TRUE, cluster rows using
#'   hierarchical clustering. If an hclust object, use it directly.
#' @param cluster_cols Logical or hclust object. Same as cluster_rows but for
#'   columns.
#' @param clustering_distance_rows Distance method for row clustering.
#'   Default is "euclidean".
#' @param clustering_distance_cols Distance method for column clustering.
#'   Default is "euclidean".
#' @param clustering_method Agglomeration method for hierarchical clustering.
#'   Default is "complete".
#' @param cutree_rows Number of row clusters to create via cutree, generating
#'   gap lines. Default is NA (no cutting).
#' @param cutree_cols Number of column clusters to create via cutree.
#'   Default is NA (no cutting).
#' @param treeheight_row Relative height of the row dendrogram in the layout.
#'   Default is 1.
#' @param treeheight_col Relative height of the column dendrogram in the layout.
#'   Default is 1.
#' @param show_rownames Logical, whether to display row names. Default is TRUE.
#' @param show_colnames Logical, whether to display column names. Default is TRUE.
#' @param annotation_row Data frame for row annotations. Row names must match
#'   \code{mat} row names.
#' @param annotation_col Data frame for column annotations. Row names must match
#'   \code{mat} column names.
#' @param annotation_colors Named list of color vectors for annotation variables.
#'   Names must match column names in annotation_row/annotation_col.
#' @param annotation_legend Logical, whether to show annotation legends.
#'   Default is TRUE.
#' @param main Character string for the plot title. Default is NULL.
#' @param fontsize Base font size in points. Default is 10.
#' @param fontsize_row Font size for row labels. Default is fontsize.
#' @param fontsize_col Font size for column labels. Default is fontsize.
#' @param cellwidth Fixed cell width in points. Default is NULL (auto).
#' @param cellheight Fixed cell height in points. Default is NULL (auto).
#' @param angle_col Column label rotation angle. Default is 270.
#' @param display_numbers Logical or matrix. If TRUE, display cell values.
#'   If a matrix, display those values instead.
#' @param number_format Format string for cell numbers. Default is "\%.2f".
#' @param number_color Color for cell number text. Default is "grey30".
#' @param gaps_row Integer vector of row positions for gap lines.
#' @param gaps_col Integer vector of column positions for gap lines.
#' @param border_color Cell border color. Default is "grey60".
#' @param na_col Color for NA values. Default is "\#DDDDDD".
#' @param legend Logical, whether to show the color legend. Default is TRUE.
#' @param filename File path to save the plot. Default is NULL (no saving).
#' @param width Plot width in inches for saving. Default is 7.
#' @param height Plot height in inches for saving. Default is 7.
#' @param silent Logical, whether to suppress plot display. Default is FALSE.
#'
#' @return Invisibly returns a list containing:
#'   \describe{
#'     \item{plot}{The patchwork ggplot object}
#'     \item{hc_row}{Row hclust object (or NULL)}
#'     \item{hc_col}{Column hclust object (or NULL)}
#'     \item{row_ord}{Row order after clustering}
#'     \item{col_ord}{Column order after clustering}
#'     \item{mat_ord}{The ordered matrix}
#'   }
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' set.seed(42)
#' mat <- matrix(rnorm(60), nrow = 6,
#'   dimnames = list(paste0("Gene", 1:6), paste0("Sample", 1:10)))
#'
#' # Basic usage
#' if (interactive()) {
#'   bioRtools::gg_heatmap(mat)
#' }
gg_heatmap <- function(
    mat,
    color             = grDevices::colorRampPalette(c("#0C6291", "white", "#A63446"))(100),
    breaks            = NULL,
    scale             = c("none", "row", "column"),
    cluster_rows      = TRUE,
    cluster_cols      = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    cutree_rows       = NA,
    cutree_cols       = NA,
    treeheight_row    = 1,
    treeheight_col    = 1,
    show_rownames     = TRUE,
    show_colnames     = TRUE,
    annotation_row    = NULL,
    annotation_col    = NULL,
    annotation_colors = NULL,
    annotation_legend = TRUE,
    main              = NULL,
    fontsize          = 10,
    fontsize_row      = fontsize,
    fontsize_col      = fontsize,
    cellwidth         = NULL,
    cellheight        = NULL,
    angle_col         = 270,
    display_numbers   = FALSE,
    number_format     = "%.2f",
    number_color      = "grey30",
    gaps_row          = NULL,
    gaps_col          = NULL,
    border_color      = "grey60",
    na_col            = "#DDDDDD",
    legend            = TRUE,
    filename          = NULL,
    width             = 7,
    height            = 7,
    silent            = FALSE
) {
  scale <- match.arg(scale)
  mat <- as.matrix(mat)

  # ── Scale ──────────────────────────────────────────────────────────────────
  if (scale == "row") {
    mat <- t(base::apply(mat, 1, function(x) {
      s <- stats::sd(x, na.rm = TRUE)
      if (s == 0) s <- 1
      (x - base::mean(x, na.rm = TRUE)) / s
    }))
  } else if (scale == "column") {
    mat <- base::apply(mat, 2, function(x) {
      s <- stats::sd(x, na.rm = TRUE)
      if (s == 0) s <- 1
      (x - base::mean(x, na.rm = TRUE)) / s
    })
  }

  # ── Clustering ─────────────────────────────────────────────────────────────
  if (isTRUE(cluster_rows) || inherits(cluster_rows, "hclust")) {
    hc_row <- if (inherits(cluster_rows, "hclust")) cluster_rows else {
      stats::hclust(stats::dist(mat, method = clustering_distance_rows),
                    method = clustering_method)
    }
    row_ord <- hc_row$order
  } else {
    hc_row <- NULL
    row_ord <- seq_len(nrow(mat))
    treeheight_row <- 0
  }

  if (isTRUE(cluster_cols) || inherits(cluster_cols, "hclust")) {
    hc_col <- if (inherits(cluster_cols, "hclust")) cluster_cols else {
      stats::hclust(stats::dist(t(mat), method = clustering_distance_cols),
                    method = clustering_method)
    }
    col_ord <- hc_col$order
  } else {
    hc_col <- NULL
    col_ord <- seq_len(ncol(mat))
    treeheight_col <- 0
  }

  mat_ord <- mat[row_ord, col_ord, drop = FALSE]

  # ── Cutree gaps ────────────────────────────────────────────────────────────
  if (!is.na(cutree_rows) && !is.null(hc_row)) {
    cl <- stats::cutree(hc_row, k = cutree_rows)[hc_row$order]
    gaps_row <- which(diff(cl) != 0)
  }
  if (!is.na(cutree_cols) && !is.null(hc_col)) {
    cl <- stats::cutree(hc_col, k = cutree_cols)[hc_col$order]
    gaps_col <- which(diff(cl) != 0)
  }

  # ── Dimensions and labels ──────────────────────────────────────────────────
  nr <- nrow(mat_ord)
  nc <- ncol(mat_ord)
  rn <- rownames(mat_ord)
  cn <- colnames(mat_ord)
  if (is.null(rn)) rn <- paste0("R", seq_len(nr))
  if (is.null(cn)) cn <- paste0("C", seq_len(nc))

  # ── Breaks ─────────────────────────────────────────────────────────────────
  if (is.null(breaks)) {
    rng <- range(mat_ord, na.rm = TRUE)
    breaks <- seq(rng[1], rng[2], length.out = length(color) + 1)
  }

  # ── Long format ────────────────────────────────────────────────────────────
  df <- data.frame(
    row   = factor(rep(rn, times = nc), levels = rev(rn)),
    col   = factor(rep(cn, each  = nr), levels = cn),
    value = as.vector(mat_ord),
    stringsAsFactors = FALSE
  )

  # ── Main heatmap ───────────────────────────────────────────────────────────
  p_heat <- ggplot2::ggplot(df, ggplot2::aes(x = col, y = row, fill = value)) +
    ggplot2::geom_tile(color = border_color, linewidth = 0.3, na.rm = TRUE) +
    ggplot2::scale_fill_gradientn(
      colors   = color,
      values   = scales::rescale(breaks),
      limits   = range(breaks),
      na.value = na_col,
      name     = ""
    ) +
    ggplot2::scale_x_discrete(expand = c(0, 0), position = "top") +
    ggplot2::scale_y_discrete(expand = c(0, 0))

  # Display numbers
  if (isTRUE(display_numbers)) {
    df$label <- sprintf(number_format, df$value)
    p_heat <- p_heat +
      ggplot2::geom_text(
        ggplot2::aes(label = label),
        color = number_color,
        size = fontsize * 0.3,
        na.rm = TRUE,
        show.legend = FALSE
      )
  } else if (is.matrix(display_numbers) || is.data.frame(display_numbers)) {
    nums <- as.matrix(display_numbers)[row_ord, col_ord, drop = FALSE]
    df$label <- sprintf(number_format, as.vector(nums))
    p_heat <- p_heat +
      ggplot2::geom_text(
        ggplot2::aes(label = label),
        color = number_color,
        size = fontsize * 0.3,
        na.rm = TRUE,
        show.legend = FALSE
      )
  }

  # Gap lines
  if (!is.null(gaps_row) && length(gaps_row) > 0) {
    p_heat <- p_heat +
      ggplot2::geom_hline(
        yintercept = nr - gaps_row + 0.5,
        color = "white", linewidth = 1.5
      )
  }
  if (!is.null(gaps_col) && length(gaps_col) > 0) {
    p_heat <- p_heat +
      ggplot2::geom_vline(
        xintercept = gaps_col + 0.5,
        color = "white", linewidth = 1.5
      )
  }

  # Axis text angle
  hjust <- if (angle_col == 0) 0.5 else 1
  vjust <- if (angle_col %in% c(270, 315)) 0.5 else 1

  p_heat <- p_heat +
    ggplot2::theme_minimal(base_size = fontsize) +
    ggplot2::theme(
      axis.text.x  = if (show_colnames)
        ggplot2::element_text(angle = angle_col, hjust = hjust, vjust = vjust,
                              size = fontsize_col)
      else ggplot2::element_blank(),
      axis.text.y  = if (show_rownames)
        ggplot2::element_text(size = fontsize_row, hjust = 1)
      else ggplot2::element_blank(),
      axis.ticks   = ggplot2::element_blank(),
      axis.title   = ggplot2::element_blank(),
      panel.grid   = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      legend.position = if (legend) "right" else "none",
      plot.margin = ggplot2::margin(2, 2, 2, 2)
    )

  if (!is.null(main)) {
    p_heat <- p_heat +
      ggplot2::ggtitle(main) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = fontsize + 2)
      )
  }

  # ── Dendrograms ────────────────────────────────────────────────────────────
  p_col_tree <- if (!is.null(hc_col) && treeheight_col > 0)
    heatmap_col_dendro(hc_col) else NULL
  p_row_tree <- if (!is.null(hc_row) && treeheight_row > 0)
    heatmap_row_dendro(hc_row) else NULL

  # ── Annotations ────────────────────────────────────────────────────────────
  ann_col_plots <- if (!is.null(annotation_col))
    heatmap_annotations(annotation_col, col_ord, cn, "col",
                        annotation_colors, annotation_legend, na_col)
  else list()

  ann_row_plots <- if (!is.null(annotation_row))
    heatmap_annotations(annotation_row, row_ord, rev(rn), "row",
                        annotation_colors, annotation_legend, na_col)
  else list()

  # ── Build layout ───────────────────────────────────────────────────────────
  final_plot <- heatmap_layout(
    p_heat, p_col_tree, p_row_tree,
    ann_col_plots, ann_row_plots,
    treeheight_col, treeheight_row
  )

  # ── Output ─────────────────────────────────────────────────────────────────
  if (!is.null(cellwidth) || !is.null(cellheight)) {
    width  <- if (!is.null(cellwidth))  cellwidth  * nc / 72 else width
    height <- if (!is.null(cellheight)) cellheight * nr / 72 else height
  }

  if (!is.null(filename)) {
    ggplot2::ggsave(filename, plot = final_plot, width = width, height = height)
  } else if (!silent) {
    print(final_plot)
  }

  invisible(list(
    plot    = final_plot,
    hc_row  = hc_row,
    hc_col  = hc_col,
    row_ord = row_ord,
    col_ord = col_ord,
    mat_ord = mat_ord
  ))
}


# ── Internal helpers ─────────────────────────────────────────────────────────

# Extract dendrogram line segments from an hclust object.
# Returns a data.frame with columns x, y, xend, yend.
# @param hc An hclust object.
# @noRd
heatmap_dendro_segments <- function(hc) {
  merge  <- hc$merge
  height <- hc$height
  n <- length(hc$order)

  leaf_x <- numeric(n)
  leaf_x[hc$order] <- seq_len(n)

  node_x <- numeric(nrow(merge))
  node_y <- numeric(nrow(merge))

  get_x <- function(i) if (i < 0) leaf_x[-i] else node_x[i]
  get_y <- function(i) if (i < 0) 0 else node_y[i]

  segs <- vector("list", nrow(merge))
  for (i in seq_len(nrow(merge))) {
    l <- merge[i, 1]; r <- merge[i, 2]
    lx <- get_x(l); rx <- get_x(r)
    ly <- get_y(l); ry <- get_y(r)
    my <- height[i]

    segs[[i]] <- data.frame(
      x    = c(lx, lx, rx),
      y    = c(ly, my, my),
      xend = c(lx, rx, rx),
      yend = c(my, my, ry)
    )
    node_x[i] <- (lx + rx) / 2
    node_y[i] <- my
  }
  do.call(rbind, segs)
}

# Column dendrogram (placed above heatmap, leaves pointing down).
# @param hc An hclust object.
# @noRd
heatmap_col_dendro <- function(hc) {
  n <- length(hc$order)
  segs <- heatmap_dendro_segments(hc)
  ggplot2::ggplot(segs) +
    ggplot2::geom_segment(
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      linewidth = 0.4, lineend = "round"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0.5, n + 0.5), expand = ggplot2::expansion(add = 0)
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = 0)) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(2, 2, 0, 2))
}

# Row dendrogram (placed left of heatmap, leaves pointing right).
# @param hc An hclust object.
# @noRd
heatmap_row_dendro <- function(hc) {
  n <- length(hc$order)
  segs <- heatmap_dendro_segments(hc)
  ggplot2::ggplot(segs) +
    ggplot2::geom_segment(
      ggplot2::aes(x = y, y = x, xend = yend, yend = xend),
      linewidth = 0.4, lineend = "round"
    ) +
    ggplot2::scale_y_reverse(
      limits = c(n + 0.5, 0.5),
      expand = ggplot2::expansion(add = 0)
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(add = 0)) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 2, 2, 2))
}

# Build annotation strip plots.
# Returns a list of ggplot objects, one per annotation variable.
# @param ann_df Annotation data.frame.
# @param ord Row/column order after clustering.
# @param labels Factor levels matching heatmap axis.
# @param orientation "col" for column annotations, "row" for row annotations.
# @param user_colors Named list of user-provided color vectors.
# @param show_legend Logical, whether to show legend.
# @param na_col Color for NA values.
# @noRd
heatmap_annotations <- function(ann_df, ord, labels, orientation,
                                user_colors, show_legend, na_col) {
  if (is.null(ann_df)) return(list())

  ann_df <- ann_df[ord, , drop = FALSE]

  pal <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
           "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62",
           "#E5C494", "#B3B3B3")

  plots <- lapply(names(ann_df), function(var) {
    vals <- ann_df[[var]]
    if (orientation == "col") {
      adf <- data.frame(id = factor(labels, levels = labels), val = vals)
    } else {
      adf <- data.frame(id = factor(labels, levels = labels), val = vals)
    }

    if (is.numeric(vals)) {
      p <- if (orientation == "col") {
        ggplot2::ggplot(adf, ggplot2::aes(x = id, y = 1, fill = val)) +
          ggplot2::geom_tile(color = NA)
      } else {
        ggplot2::ggplot(adf, ggplot2::aes(x = 1, y = id, fill = val)) +
          ggplot2::geom_tile(color = NA)
      }
      p <- p +
        ggplot2::scale_fill_distiller(
          palette = "Blues", direction = 1, name = var, na.value = na_col
        )
    } else {
      uv <- unique(na.omit(vals))
      clrs <- if (!is.null(user_colors[[var]])) {
        user_colors[[var]]
      } else {
        stats::setNames(pal[seq_along(uv)], uv)
      }

      p <- if (orientation == "col") {
        ggplot2::ggplot(adf, ggplot2::aes(x = id, y = 1, fill = val)) +
          ggplot2::geom_tile(color = NA)
      } else {
        ggplot2::ggplot(adf, ggplot2::aes(x = 1, y = id, fill = val)) +
          ggplot2::geom_tile(color = NA)
      }
      p <- p +
        ggplot2::scale_fill_manual(
          values = clrs, name = var, na.value = na_col, drop = FALSE
        )
    }

    if (orientation == "col") {
      p <- p + ggplot2::scale_x_discrete(expand = c(0, 0))
    } else {
      p <- p + ggplot2::scale_y_discrete(expand = c(0, 0))
    }

    p + ggplot2::theme_void() +
      ggplot2::theme(
        legend.position = if (show_legend) "right" else "none",
        plot.margin = ggplot2::margin(1, 1, 1, 1)
      )
  })

  plots
}

# Assemble the final layout using patchwork design strings.
# The layout is a grid:
#   Row 1..k  : [empty] [col_tree / col_ann]   (above heatmap)
#   Last row  : [row_tree] [row_ann] [heatmap]  (main row)
# @noRd
heatmap_layout <- function(p_heat, p_col_tree, p_row_tree,
                           ann_col_plots, ann_row_plots,
                           treeheight_col, treeheight_row) {
  has_col_tree <- !is.null(p_col_tree)
  has_row_tree <- !is.null(p_row_tree)
  n_col_ann <- length(ann_col_plots)
  n_row_ann <- length(ann_row_plots)

  # If nothing extra, just return the heatmap
  if (!has_col_tree && !has_row_tree && n_col_ann == 0 && n_row_ann == 0) {
    return(p_heat)
  }

  # Count grid dimensions
  left_cols <- has_row_tree + n_row_ann  # columns left of heatmap
  total_cols <- left_cols + 1            # +1 for heatmap
  top_rows <- has_col_tree + n_col_ann   # rows above heatmap
  total_rows <- top_rows + 1             # +1 for heatmap row

  # Assign plot letters and build design matrix
  design_mat <- matrix("#", nrow = total_rows, ncol = total_cols)
  plots <- list()
  idx <- 1L

  # Heatmap always at bottom-right
  heat_letter <- LETTERS[idx]
  design_mat[total_rows, total_cols] <- heat_letter
  plots[[idx]] <- p_heat
  idx <- idx + 1L

  # Row tree at bottom-left
  if (has_row_tree) {
    design_mat[total_rows, 1] <- LETTERS[idx]
    plots[[idx]] <- p_row_tree
    idx <- idx + 1L
  }

  # Row annotations next to row tree
  for (i in seq_len(n_row_ann)) {
    col_pos <- has_row_tree + i
    design_mat[total_rows, col_pos] <- LETTERS[idx]
    plots[[idx]] <- ann_row_plots[[i]]
    idx <- idx + 1L
  }

  # Column tree at top of heatmap column
  row_pos <- 1L
  if (has_col_tree) {
    design_mat[row_pos, total_cols] <- LETTERS[idx]
    plots[[idx]] <- p_col_tree
    row_pos <- row_pos + 1L
    idx <- idx + 1L
  }

  # Column annotations below col tree
  for (i in seq_len(n_col_ann)) {
    design_mat[row_pos, total_cols] <- LETTERS[idx]
    plots[[idx]] <- ann_col_plots[[i]]
    row_pos <- row_pos + 1L
    idx <- idx + 1L
  }

  # Build design string
  design_str <- paste(apply(design_mat, 1, paste, collapse = ""), collapse = "\n")

  # Build widths and heights
  widths <- c()
  if (has_row_tree) widths <- c(widths, treeheight_row)
  widths <- c(widths, rep(0.4, n_row_ann))
  widths <- c(widths, 10)  # heatmap

  heights <- c()
  if (has_col_tree) heights <- c(heights, treeheight_col)
  heights <- c(heights, rep(0.3, n_col_ann))
  heights <- c(heights, 10)  # heatmap

  Reduce(`+`, plots) +
    patchwork::plot_layout(
      design = design_str,
      widths = widths,
      heights = heights,
      guides = "collect"
    )
}
