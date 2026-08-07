#' CIM/QTL LOD Score Plot
#'
#' @description
#' Create a Manhattan-style LOD score plot for Composite Interval Mapping
#' (CIM) or other QTL scan results. Computes cumulative physical positions
#' across chromosomes, applies alternating chromosome coloring, and
#' optionally draws threshold lines.
#'
#' @param data A data.frame, data.table, or file path (character string).
#'   If a file path, it is read with \code{data.table::fread()}.
#'   Must contain chromosome, position, and LOD score columns.
#' @param threshold Numeric value for a genome-wide horizontal threshold
#'   line. Default is NULL (no line).
#' @param chr_col Character. Column name for chromosome. Default is "chr".
#' @param pos_col Character. Column name for physical position.
#'   Default is "pos".
#' @param lod_col Character. Column name for LOD scores. Default is "lod".
#' @param pos_scale Numeric. Multiplier applied to the position column to
#'   convert to base pairs. Default is 1e6 for Mbp → bp.
#'   Set to 1 if positions are already in bp.
#' @param lod_offset Numeric. Value added to LOD scores to avoid log-scale
#'   issues with zero values. Default is 0.1.
#' @param chr_prefix Character. Prefix added to chromosome names for
#'   correct sort order (e.g., "Chr1" > "Chr10" > "Chr2" without padding).
#'   Stripped from display labels. Default is "Chr".
#' @param title Character. Plot title. Default is NULL.
#' @param x_lab Character. X-axis label. Default is "Chromosome".
#' @param y_lab Character. Y-axis label. Default is "LOD".
#' @param colors Character vector. Alternating chromosome colors, recycled
#'   across chromosomes. Default is c("grey60", "grey30").
#' @param threshold_color Character. Threshold line color. Default is "red".
#' @param threshold_linetype Character. Threshold line type.
#'   Default is "dashed".
#' @param threshold_linewidth Numeric. Threshold line width.
#'   Default is 0.8.
#' @param y_max Numeric. Y-axis upper limit. Default is 100.
#' @param y_breaks Numeric vector. Y-axis tick positions.
#'   Default is seq(0, 100, 20).
#' @param line_size Numeric. Line width for LOD curves. Default is 0.5.
#'
#' @return A list containing:
#' \describe{
#'   \item{plot}{A ggplot object of the LOD score plot.}
#'   \item{data}{The processed data frame with cumulative positions.}
#'   \item{chr_centers}{Data frame of chromosome center positions for
#'     custom axis labeling.}
#' }
#'
#' @details
#' The function performs the following preprocessing steps:
#' \itemize{
#'   \item Adds \code{chr_prefix} to chromosome names for natural sort order
#'   \item Scales positions by \code{pos_scale}
#'   \item Filters rows with NA LOD scores
#'   \item Adds \code{lod_offset} to LOD scores
#'   \item Computes cumulative positions for genome-wide x-axis
#'   \item Strips \code{chr_prefix} from display labels
#' }
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage with a data frame
#' lod_data <- data.frame(
#'   chr = rep(1:5, each = 100),
#'   pos = rep(seq(0, 99), 5) * 1e6,
#'   lod = abs(rnorm(500, 5, 2))
#' )
#' result <- plot_cim_lod(lod_data, threshold = 3.0)
#' result$plot
#'
#' # From file
#' result <- plot_cim_lod("cim_lod_data.tsv",
#'   threshold = 3.5,
#'   y_max = 80,
#'   y_breaks = seq(0, 80, 10)
#' )
#' result$plot
#' }
plot_cim_lod <- function(data,
                         threshold = NULL,
                         chr_col = "chr",
                         pos_col = "pos",
                         lod_col = "lod",
                         pos_scale = 1e6,
                         lod_offset = 0.1,
                         chr_prefix = "Chr",
                         title = NULL,
                         x_lab = "Chromosome",
                         y_lab = "LOD",
                         colors = c("grey60", "grey30"),
                         threshold_color = "red",
                         threshold_linetype = "dashed",
                         threshold_linewidth = 0.8,
                         y_max = 100,
                         y_breaks = seq(0, 100, 20),
                         line_size = 0.5) {

  # ── Read file if path provided ──────────────────────────────────────────────
  if (is.character(data) && length(data) == 1L) {
    if (!file.exists(data)) {
      stop("File not found: ", data)
    }
    data <- data.table::fread(data)
  }

  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame, data.table, or a file path (character)")
  }

  # ── Check required columns ──────────────────────────────────────────────────
  required_cols <- c(chr_col, pos_col, lod_col)
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0L) {
    stop("Required columns missing: ", paste(missing_cols, collapse = ", "))
  }

  # ── Validate numeric parameters ─────────────────────────────────────────────
  if (!is.numeric(pos_scale) || pos_scale <= 0) {
    stop("'pos_scale' must be a positive number")
  }

  if (!is.numeric(y_max) || y_max <= 0) {
    stop("'y_max' must be a positive number")
  }

  # ── Preprocess data ─────────────────────────────────────────────────────────
  df <- data
  df[["chr"]] <- paste0(chr_prefix, df[[chr_col]])
  df[["pos"]] <- as.numeric(df[[pos_col]]) * pos_scale
  df[["lod"]] <- as.numeric(df[[lod_col]])

  df <- df %>%
    dplyr::filter(!is.na(.data$lod)) %>%
    dplyr::mutate(lod = .data$lod + lod_offset)

  if (nrow(df) == 0L) {
    stop("No data remaining after filtering NA LOD values")
  }

  # ── Compute cumulative positions ────────────────────────────────────────────
  chr_len <- df %>%
    dplyr::group_by(.data$chr) %>%
    dplyr::summarise(chr_len = max(.data$pos), .groups = "drop")

  chr_offset <- chr_len %>%
    dplyr::mutate(
      offset = dplyr::lag(cumsum(.data$chr_len), default = 0)
    ) %>%
    dplyr::select("chr", "offset")

  df_plot <- df %>%
    dplyr::left_join(chr_offset, by = "chr") %>%
    dplyr::mutate(pos_cum = .data$pos + .data$offset) %>%
    dplyr::mutate(chr = stringr::str_remove(.data$chr, chr_prefix))

  # ── Compute chromosome centers for axis labels ──────────────────────────────
  chr_center <- df_plot %>%
    dplyr::group_by(.data$chr) %>%
    dplyr::summarise(
      center = (min(.data$pos_cum) + max(.data$pos_cum)) / 2,
      .groups = "drop"
    )

  # ── Build color mapping ─────────────────────────────────────────────────────
  n_chr <- dplyr::n_distinct(df_plot$chr)
  chr_colors <- rep(colors, length.out = n_chr)
  names(chr_colors) <- unique(df_plot$chr)

  # ── Build ggplot ────────────────────────────────────────────────────────────
  p <- ggplot2::ggplot(
    df_plot,
    ggplot2::aes(x = .data$pos_cum, y = .data$lod, color = .data$chr)
  ) +
    ggplot2::geom_line(linewidth = line_size) +
    ggplot2::scale_x_continuous(
      labels = chr_center$chr,
      breaks = chr_center$center
    ) +
    ggplot2::scale_y_continuous(
      expand = c(0, 0),
      limits = c(0, y_max),
      breaks = y_breaks
    ) +
    ggplot2::scale_color_manual(values = chr_colors) +
    ggplot2::labs(
      title = title,
      x     = x_lab,
      y     = y_lab
    ) +
    ggprism::theme_prism() +
    ggplot2::theme(legend.position = "none")

  # ── Add threshold line if requested ─────────────────────────────────────────
  if (!is.null(threshold)) {
    if (!is.numeric(threshold) || length(threshold) != 1L) {
      stop("'threshold' must be a single numeric value")
    }
    p <- p +
      ggplot2::geom_hline(
        yintercept = threshold,
        color = threshold_color,
        linetype = threshold_linetype,
        linewidth = threshold_linewidth
      )
  }

  # ── Return ──────────────────────────────────────────────────────────────────
  list(
    plot        = p,
    data        = df_plot,
    chr_centers = chr_center
  )
}
