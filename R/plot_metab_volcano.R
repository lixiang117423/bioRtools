#' Multi-Group Volcano Plot for Metabolomics Differential Results
#'
#' Manhattan-style multi-group volcano plot tailored for metabolomics
#' differential analysis results (e.g. the \code{differential_analysis} slot of
#' \code{\link{opls_analysis}}). Metabolomics comparisons typically contain a
#' large proportion of non-significant features; drawing them as fixed
#' \code{geom_point()} would collapse them into a single vertical line. This
#' function therefore jitters both the NS and the Up/Down features horizontally
#' (the y positions stay exact, so log2FC values are never distorted) and draws
#' the significant Up/Down features as larger points on top of the NS layer.
#'
#' Features whose adjusted p-value is \code{NA} (e.g. untestable
#' constant/zero-variance metabolites) are kept visible as NS instead of being
#' silently dropped, which is common in metabolomics matrices.
#'
#' The Up/Down/NS coloring follows the \code{regulation} column when it is
#' present (as in \code{opls_analysis} differential output, where it is
#' VIP-based by default), so the colored points match the Up/Down counts of
#' the differential table exactly. Only when no regulation column is available
#' are \code{log2fc_thr} and \code{padj_thr} used to classify points.
#'
#' @param data A data frame. Two formats are auto-detected:
#'   \itemize{
#'     \item \code{opls_analysis} differential output: columns
#'       \code{comparison}, \code{feature_id}, \code{log2FoldChange}, \code{padj},
#'       \code{regulation}
#'     \item Seurat-style cluster table: columns \code{cluster}, \code{gene},
#'       \code{avg_log2FC}, \code{p_val_adj}
#'   }
#' @param group_col Column for grouping (default: auto-detect).
#' @param feature_col Column for feature labels (default: auto-detect).
#' @param log2fc_col Column for log2 fold change (default: auto-detect).
#' @param padj_col Column for adjusted p-values (default: auto-detect).
#' @param regulation_col Column with Up/Down/NS labels (default:
#'   auto-detect \code{"regulation"}). When present, this column decides the
#'   coloring; \code{log2fc_thr} and \code{padj_thr} then only apply when no
#'   regulation column is available.
#' @param log2fc_thr Absolute log2FC threshold for significance coloring
#'   (default: 1). Only used when no regulation column is available.
#' @param padj_thr Adjusted p-value threshold for significance coloring
#'   (default: 0.05). Only used when no regulation column is available.
#' @param top_n Number of top features to label per group per direction
#'   (default: 5).
#' @param label_size Font size for feature labels (default: 2.5).
#' @param jitter_width Horizontal jitter width in x-axis units for NS points
#'   (default: 0.4). Vertical jitter is always 0 so log2FC values stay exact.
#' @param sig_jitter_width Horizontal jitter width for Up/Down points
#'   (default: same as \code{jitter_width}). Set to 0 to keep them on the
#'   exact group-center line. Vertical jitter is always 0.
#' @param jitter_seed Seed for reproducible jitter positions (default: 123).
#'   Set to \code{NULL} for random jitter. The global RNG state is not
#'   modified.
#' @param ns_color Fill color for NS points (default: "grey75").
#' @param ns_size Point size for NS points (default: 0.8).
#' @param ns_alpha Alpha for NS points (default: 0.6).
#' @param sig_size Point size for significant points (default: 1.8).
#' @param sig_alpha Alpha for significant points (default: 1).
#' @param up_color Color for up-regulated points (default: "#0073C2FF").
#' @param down_color Color for down-regulated points (default: "#EE0000FF").
#' @param group_colors Fill colors for background bands, or a palette name
#'   from \code{academic_db$research} (e.g. "default",
#'   "nature_genomics"). Default: "default". Recycled if needed.
#' @param band_alpha Alpha for center bands (default: 0.3).
#' @param y_limits Y-axis limits. If NULL, auto-computed.
#' @param y_breaks Y-axis breaks. If NULL, auto-computed.
#' @param legend_position Legend position inside the panel
#'   (default: \code{c(0.08, 0.9)}).
#'
#' @return A ggplot object.
#'
#' @seealso \code{\link{plot_multi_volcano}} for the RNA-seq-oriented version,
#'   \code{\link{opls_analysis}} for the metabolomics OPLS-DA differential
#'   results feeding this plot.
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' \dontrun{
#' # opls_analysis differential output feeds the plot directly
#' opls_res <- opls_analysis(
#'   data = df_metab, sample = df_sample,
#'   group_col = "group", ref_group = "T1"
#' )
#' p <- plot_metab_volcano(opls_res$differential_analysis)
#' p
#' }
plot_metab_volcano <- function(data,
                               group_col = NULL,
                               feature_col = NULL,
                               log2fc_col = NULL,
                               padj_col = NULL,
                               regulation_col = NULL,
                               log2fc_thr = 1,
                               padj_thr = 0.05,
                               top_n = 5,
                               label_size = 2.5,
                               jitter_width = 0.4,
                               sig_jitter_width = jitter_width,
                               jitter_seed = 123,
                               ns_color = "grey75",
                               ns_size = 0.8,
                               ns_alpha = 0.6,
                               sig_size = 1.8,
                               sig_alpha = 1,
                               up_color = "#0073C2FF",
                               down_color = "#EE0000FF",
                               group_colors = "default",
                               band_alpha = 0.3,
                               y_limits = NULL,
                               y_breaks = NULL,
                               legend_position = c(0.08, 0.9)) {

  # 1. Input validation
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  if (!is.null(regulation_col) &&
    (!is.character(regulation_col) || length(regulation_col) != 1)) {
    stop("'regulation_col' must be NULL or a single column name")
  }
  if (!is.numeric(log2fc_thr) || length(log2fc_thr) != 1 || log2fc_thr < 0) {
    stop("'log2fc_thr' must be a single non-negative number")
  }
  if (!is.numeric(padj_thr) || length(padj_thr) != 1 ||
    padj_thr <= 0 || padj_thr >= 1) {
    stop("'padj_thr' must be a single number between 0 and 1")
  }
  if (!is.numeric(top_n) || length(top_n) != 1 ||
    top_n < 0 || top_n != round(top_n)) {
    stop("'top_n' must be a single non-negative integer")
  }
  if (!is.numeric(jitter_width) || length(jitter_width) != 1 ||
    jitter_width < 0) {
    stop("'jitter_width' must be a single non-negative number")
  }
  if (!is.numeric(sig_jitter_width) || length(sig_jitter_width) != 1 ||
    sig_jitter_width < 0) {
    stop("'sig_jitter_width' must be a single non-negative number")
  }
  if (!is.null(jitter_seed) &&
    (!is.numeric(jitter_seed) || length(jitter_seed) != 1)) {
    stop("'jitter_seed' must be NULL or a single number")
  }
  for (p in c("ns_alpha", "sig_alpha", "band_alpha")) {
    v <- get(p)
    if (!is.numeric(v) || length(v) != 1 || v < 0 || v > 1) {
      stop("'", p, "' must be a single number between 0 and 1")
    }
  }
  for (p in c("ns_size", "sig_size", "label_size")) {
    v <- get(p)
    if (!is.numeric(v) || length(v) != 1 || v < 0) {
      stop("'", p, "' must be a single non-negative number")
    }
  }

  # --- Auto-detect columns ---------------------------------------------------
  cols <- detect_volcano_cols(data, group_col, feature_col, log2fc_col, padj_col)

  # --- Resolve regulation column --------------------------------------------
  # opls_analysis differential output carries a regulation column (VIP-based
  # by default); when present it drives the coloring so the plot matches the
  # Up/Down/NS counts of the differential table exactly.
  if (is.null(regulation_col)) {
    if ("regulation" %in% colnames(data)) regulation_col <- "regulation"
  } else if (!regulation_col %in% colnames(data)) {
    stop("'regulation_col' (", regulation_col, ") not found in 'data'")
  }

  # --- Prepare data ----------------------------------------------------------
  keep <- !is.na(data[[cols$log2fc_col]])
  df <- data.frame(
    group   = as.character(data[[cols$group_col]]),
    feature = as.character(data[[cols$feature_col]]),
    log2FC  = as.numeric(data[[cols$log2fc_col]]),
    padj    = as.numeric(data[[cols$padj_col]]),
    stringsAsFactors = FALSE
  )[keep, ]

  lvls <- unique(df$group)
  df$group <- factor(df$group, levels = lvls)
  n_groups <- length(lvls)

  # Classify points
  if (!is.null(regulation_col)) {
    reg <- as.character(data[[regulation_col]])[keep]
    df$type <- ifelse(!is.na(reg) & reg %in% c("Up", "Down"), reg, "NS")
  } else {
    # NA padj (untestable features) stay visible as NS instead of being
    # silently dropped.
    df$type <- ifelse(
      !is.na(df$padj) & abs(df$log2FC) >= log2fc_thr & df$padj < padj_thr,
      ifelse(df$log2FC > 0, "Up", "Down"),
      "NS"
    )
  }
  df$type <- factor(df$type, levels = c("NS", "Down", "Up"))

  # Resolve colors: palette name or explicit vector
  if (is.character(group_colors) && length(group_colors) == 1) {
    pal <- tryCatch(
      academic_db$research[[group_colors]],
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
  sig <- df[df$type != "NS", ]
  if (nrow(sig) > 0) {
    label_df <- do.call(rbind, lapply(split(sig, sig$group), function(sub) {
      ups   <- utils::head(sub[order(-sub$log2FC), ], top_n)
      downs <- utils::head(sub[order(sub$log2FC), ], top_n)
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
    # NS points: horizontal jitter only (height = 0 keeps log2FC exact),
    # drawn first so the Up/Down layer below always stays on top.
    ggplot2::geom_jitter(
      data = df[df$type == "NS", ],
      ggplot2::aes(x = group, y = log2FC),
      position = ggplot2::position_jitter(width = jitter_width, height = 0,
                                          seed = jitter_seed),
      color = ns_color, size = ns_size, alpha = ns_alpha, stroke = 0,
      show.legend = FALSE
    ) +
    # Significant points: horizontal jitter too (height = 0 keeps log2FC
    # exact), drawn after NS so they always stay on top, with legend.
    ggplot2::geom_jitter(
      data = df[df$type != "NS", ],
      ggplot2::aes(x = group, y = log2FC, color = type),
      position = ggplot2::position_jitter(width = sig_jitter_width, height = 0,
                                          seed = jitter_seed),
      size = sig_size, alpha = sig_alpha, stroke = 0
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
      limits = y_limits, breaks = y_breaks
    ) +
    ggplot2::labs(x = NULL, y = "log2FC") +
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
      legend.position.inside = legend_position,
      legend.key = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(margin = ggplot2::margin(l = 0))
    )
}
