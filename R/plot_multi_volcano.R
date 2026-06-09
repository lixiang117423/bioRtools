#' Multi-Group Volcano Plot (Manhattan Style)
#'
#' Creates a Manhattan-style volcano plot for visualizing differential analysis
#' results across multiple groups. Auto-detects input format: works with the
#' original cluster format (Seurat-style) and \code{\link{find_dams_deseq2}}
#' output directly.
#'
#' @param data A data frame. Two formats are auto-detected:
#'   \itemize{
#'     \item Original: columns \code{cluster}, \code{gene},
#'       \code{avg_log2FC}, \code{p_val_adj}
#'     \item \code{find_dams_deseq2} output: columns \code{comparison},
#'       \code{feature_id}, \code{log2FoldChange}, \code{padj}
#'   }
#' @param groupCol Column for grouping (default: auto-detect).
#' @param featureCol Column for feature labels (default: auto-detect).
#' @param log2FCCol Column for log2 fold change (default: auto-detect).
#' @param padjCol Column for adjusted p-values (default: auto-detect).
#' @param log2FC_thr Absolute log2FC threshold for significance coloring
#'   (default: 1).
#' @param padj_thr Adjusted p-value threshold for significance coloring
#'   (default: 0.05).
#' @param topN Number of top features to label per group per direction
#'   (default: 5).
#' @param label_size Font size for feature labels (default: 2.5).
#' @param group_colors Fill colors for background bands, or a palette name
#'   from \code{\link{scale_fill_research}} (e.g. "default",
#'   "nature_genomics"). Default: "default". Recycled if needed.
#' @param up_color Color for up-regulated points (default: "#0073C2FF").
#' @param down_color Color for down-regulated points (default: "#EE0000FF").
#' @param band_alpha Alpha for center bands (default: 0.3).
#' @param y_limits Y-axis limits. If NULL, auto-computed.
#' @param y_breaks Y-axis breaks. If NULL, auto-computed.
#' @param legend_position Legend position (default: \code{c(0.08, 0.9)}).
#'
#' @return A ggplot object.
#'
#' @examples
#' # With original cluster format
#' \dontrun{
#' df <- read.delim("data.tsv")
#' plot_multi_volcano(df)
#' }
#'
#' # With find_dams_deseq2 output
#' \dontrun{
#' res <- find_dams_deseq2(mat, meta, groupCol = "treatment")
#' plot_multi_volcano(res)
#' }
#'
#' @seealso \code{\link{find_dams_deseq2}}, \code{\link{find_dams_lefse}}
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
plot_multi_volcano <- function(data,
                               groupCol = NULL,
                               featureCol = NULL,
                               log2FCCol = NULL,
                               padjCol = NULL,
                               log2FC_thr = 1,
                               padj_thr = 0.05,
                               topN = 5,
                               label_size = 2.5,
                               group_colors = "default",
                               up_color = "#0073C2FF",
                               down_color = "#EE0000FF",
                               band_alpha = 0.3,
                               y_limits = NULL,
                               y_breaks = NULL,
                               legend_position = c(0.08, 0.9)) {

  # --- Auto-detect columns ---------------------------------------------------
  cols <- detect_volcano_cols(data, groupCol, featureCol, log2FCCol, padjCol)

  # --- Prepare data ----------------------------------------------------------
  df <- data.frame(
    group   = as.character(data[[cols$groupCol]]),
    feature = as.character(data[[cols$featureCol]]),
    log2FC  = as.numeric(data[[cols$log2FCCol]]),
    padj    = as.numeric(data[[cols$padjCol]]),
    stringsAsFactors = FALSE
  )
  df <- df[!is.na(df$log2FC), ]

  lvls <- unique(df$group)
  df$group <- factor(df$group, levels = lvls)
  n_groups <- length(lvls)

  # Classify points
  df$type <- ifelse(
    abs(df$log2FC) >= log2FC_thr & df$padj < padj_thr,
    ifelse(df$log2FC > 0, "Up", "Down"),
    "NS"
  )
  df$type <- factor(df$type, levels = c("NS", "Down", "Up"))

  # Resolve colors: palette name or explicit vector
  if (is.character(group_colors) && length(group_colors) == 1) {
    pal <- tryCatch(
      bioRtools:::academic_db$research[[group_colors]],
      error = function(e) NULL
    )
    if (!is.null(pal)) group_colors <- unname(pal)
  }
  group_colors <- rep_len(group_colors, n_groups)

  # --- Background rectangles -------------------------------------------------
  bg_center <- data.frame(
    group = lvls,
    xmin  = seq_along(lvls) - 0.48,
    xmax  = seq_along(lvls) + 0.48,
    ymin  = -0.5,
    ymax  = 0.5,
    stringsAsFactors = FALSE
  )
  bg_center$group <- factor(bg_center$group, levels = lvls)

  bg_vert <- do.call(rbind, lapply(split(df$log2FC, df$group), function(v) {
    c(min(v, na.rm = TRUE) - 0.1, max(v, na.rm = TRUE) + 0.1)
  }))
  bg_vert <- data.frame(
    group = lvls,
    ymin  = bg_vert[, 1],
    ymax  = bg_vert[, 2],
    xmin  = seq_along(lvls) - 0.48,
    xmax  = seq_along(lvls) + 0.48,
    stringsAsFactors = FALSE
  )
  bg_vert$group <- factor(bg_vert$group, levels = lvls)

  # --- Labels ----------------------------------------------------------------
  sig <- df[df$padj < padj_thr & abs(df$log2FC) >= log2FC_thr, ]
  if (nrow(sig) > 0) {
    label_df <- do.call(rbind, lapply(split(sig, sig$group), function(sub) {
      ups   <- utils::head(sub[order(-sub$log2FC), ], topN)
      downs <- utils::head(sub[order(sub$log2FC), ], topN)
      unique(rbind(ups, downs))
    }))
    rownames(label_df) <- NULL
    label_df$group <- factor(label_df$group, levels = lvls)
  } else {
    label_df <- df[0, ]
    label_df$group <- factor(label_df$group, levels = lvls)
  }

  # --- Y-axis ----------------------------------------------------------------
  if (is.null(y_limits)) {
    yr <- range(df$log2FC, na.rm = TRUE)
    pad <- diff(yr) * 0.05
    y_limits <- c(yr[1] - pad, yr[2] + pad)
  }
  if (is.null(y_breaks)) {
    y_breaks <- ggplot2::waiver()
  }

  # --- Build plot ------------------------------------------------------------
  ggplot2::ggplot(df, ggplot2::aes(x = group, y = log2FC)) +
    ggplot2::geom_rect(
      data = bg_vert,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE, fill = "grey95"
    ) +
    ggplot2::geom_rect(
      data = bg_center,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                   fill = group),
      inherit.aes = FALSE, alpha = band_alpha
    ) +
    ggplot2::scale_fill_manual(values = group_colors, guide = "none") +
    # NS points: grey, no legend
    ggplot2::geom_point(
      data = df[df$type == "NS", ],
      ggplot2::aes(x = group, y = log2FC),
      color = "grey70", size = 0.8, stroke = 0, show.legend = FALSE
    ) +
    # Significant points: colored, with legend
    ggplot2::geom_jitter(
      data = df[df$type != "NS", ],
      ggplot2::aes(x = group, y = log2FC, color = type),
      stroke = 0
    ) +
    ggplot2::scale_color_manual(
      values = c("Down" = down_color, "Up" = up_color),
      labels = c("Down", "Up")
    ) +
    ggrepel::geom_text_repel(
      data = label_df,
      ggplot2::aes(label = feature),
      size = label_size, color = "black", box.padding = 0.2,
      max.overlaps = 20
    ) +
    ggplot2::geom_text(
      data = bg_center,
      ggplot2::aes(x = group, y = 0, label = as.character(group)),
      size = 3, color = "white", show.legend = FALSE
    ) +
    ggplot2::scale_y_continuous(
      limits = y_limits, breaks = y_breaks,
      guide = ggprism::guide_prism_offset_minor()
    ) +
    ggplot2::labs(x = NULL, y = "average log2FC") +
    ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(size = 5, shape = 19))
    ) +
    ggprism::theme_prism(base_line_size = 0.3) +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      axis.text.x  = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.line.x  = ggplot2::element_blank(),
      axis.text.y  = ggplot2::element_text(size = 10),
      axis.title   = ggplot2::element_text(size = 11),
      legend.position = legend_position,
      legend.key = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(margin = ggplot2::margin(l = 0))
    )
}


#' @rdname plot_multi_volcano
#' @keywords internal
detect_volcano_cols <- function(data, groupCol, featureCol, log2FCCol, padjCol) {
  nms <- names(data)

  if (is.null(groupCol)) {
    groupCol <- intersect(c("comparison", "cluster"), nms)
    if (length(groupCol) == 0) stop("Cannot auto-detect group column. Specify 'groupCol'.")
    groupCol <- groupCol[1]
  }
  if (is.null(featureCol)) {
    featureCol <- intersect(c("feature_id", "gene"), nms)
    if (length(featureCol) == 0) stop("Cannot auto-detect feature column. Specify 'featureCol'.")
    featureCol <- featureCol[1]
  }
  if (is.null(log2FCCol)) {
    log2FCCol <- intersect(c("log2FoldChange", "avg_log2FC"), nms)
    if (length(log2FCCol) == 0) stop("Cannot auto-detect log2FC column. Specify 'log2FCCol'.")
    log2FCCol <- log2FCCol[1]
  }
  if (is.null(padjCol)) {
    padjCol <- intersect(c("padj", "p_val_adj"), nms)
    if (length(padjCol) == 0) stop("Cannot auto-detect padj column. Specify 'padjCol'.")
    padjCol <- padjCol[1]
  }

  list(groupCol = groupCol, featureCol = featureCol,
       log2FCCol = log2FCCol, padjCol = padjCol)
}
