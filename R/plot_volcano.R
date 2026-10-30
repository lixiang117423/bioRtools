#' Create a volcano plot for differential expression analysis
#'
#' @description
#' Create a volcano plot to visualize the results of differential expression
#' analysis, typically from DESeq2. The plot displays the relationship between
#' fold change (x-axis) and statistical significance (y-axis).
#'
#' @param data The data frame used to create the volcano plot, typically the
#'   output from DESeq2 analysis or \code{RNASeq.call_DEGs_DESeq2}.
#' @param x Character string specifying the column name for the x-axis
#'   (fold change values). Default is "log2FoldChange".
#' @param y Character string specifying the column name for the y-axis
#'   (p-values or adjusted p-values). Default is "padj".
#' @param size Character string specifying the column name used to map the
#'   size of points. Default is "baseMean".
#' @param color Character string specifying the column used to color points.
#'   The default \code{"auto"} uses \code{"regulation"} when present and
#'   otherwise falls back to \code{"group"}. Set \code{NULL} for a single
#'   color.
#' @param xintercept_abs Numeric value for the absolute fold change threshold.
#'   Vertical lines will be drawn at +/- this value. Default is 1.
#' @param yintercept Numeric value for the significance threshold (before -log10
#'   transformation). Default is 0.05.
#' @param point_alpha Numeric value between 0 and 1 for point transparency.
#'   Default is 0.8.
#' @param line_color Character string specifying the color of threshold lines.
#'   Default is "grey50".
#' @param line_type Character string specifying the line type for thresholds.
#'   Default is "dashed".
#' @param line_size Numeric value for threshold line width. Default is 0.8.
#' @param title Character string for plot title. Default is "Volcano Plot".
#' @param subtitle Character string for plot subtitle. Default is NULL.
#' @param point_size_range Numeric vector of length 2 specifying the range
#'   for point sizes. Default is c(0.5, 4).
#'
#' @return A list containing two components:
#' \describe{
#'   \item{plot_volcano}{A ggplot object that can be further customized
#'     or directly displayed.}
#'   \item{data_summary}{A summary table showing the number of genes in
#'     each category (up-regulated, down-regulated, not significant).}
#' }
#'
#' @details
#' The volcano plot is a scatter plot that shows statistical significance
#' (y-axis) versus magnitude of change (x-axis). Points are typically colored
#' based on their significance and fold change thresholds. The y-axis values
#' are automatically transformed using -log10.
#'
#' The function automatically adds:
#' \itemize{
#'   \item Horizontal line at the significance threshold
#'   \item Vertical lines at the fold change thresholds
#'   \item Color coding based on the grouping variable
#'   \item Size mapping based on expression level
#' }
#'
#' @note
#' This function is specifically designed for RNA-seq differential expression
#' results but can be adapted for other types of data with appropriate
#' column names.
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' library(dplyr)
#' library(bioRtools)
#'
#' # Load example data
#' data(df.rnaseq.plot.volcano)
#'
#' # Basic volcano plot
#' volcano_result <- plot_volcano(data = df.rnaseq.plot.volcano)
#' volcano_result$plot_volcano
#'
#' # View summary statistics
#' volcano_result$data_summary
#'
#' # Customized volcano plot with different thresholds
#' volcano_custom <- plot_volcano(
#'   data = df.rnaseq.plot.volcano,
#'   xintercept_abs = 1.5,
#'   yintercept = 0.01,
#'   title = "Differential Expression Analysis",
#'   subtitle = "Volcano Plot with Custom Thresholds"
#' )
#' volcano_custom$plot_volcano
#'
#' # Using different column names
#' plot_volcano(
#'   data = df.rnaseq.plot.volcano,
#'   x = "log2FoldChange",
#'   y = "pvalue",  # Use raw p-values instead of adjusted
#'   color = "group",
#'   size = "baseMean"
#' )
#'
plot_volcano <- function(data,
                         x = "log2FoldChange",
                         y = "padj",
                         size = "baseMean",
                         color = "auto",
                         xintercept_abs = 1,
                         yintercept = 0.05,
                         point_alpha = 0.8,
                         line_color = "grey50",
                         line_type = "dashed",
                         line_size = 0.8,
                         title = "Volcano Plot",
                         subtitle = NULL,
                         point_size_range = c(0.5, 4)) {
  # Input validation
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }

  if (identical(color, "auto")) {
    color <- if ("regulation" %in% names(data)) {
      "regulation"
    } else if ("group" %in% names(data)) {
      "group"
    } else {
      NULL
    }
  }

  # Check required columns
  required_cols <- c(x, y, size)
  missing_cols <- required_cols[!required_cols %in% colnames(data)]
  if (!is.null(color) && !color %in% colnames(data)) {
    missing_cols <- c(missing_cols, color)
  }

  if (length(missing_cols) > 0) {
    stop("Required columns missing from data: ", paste(missing_cols, collapse = ", "))
  }

  # Validate numeric parameters
  if (!is.numeric(xintercept_abs) || length(xintercept_abs) != 1 || xintercept_abs <= 0) {
    stop("'xintercept_abs' must be a single positive numeric value")
  }

  if (!is.numeric(yintercept) || length(yintercept) != 1 || yintercept <= 0 || yintercept >= 1) {
    stop("'yintercept' must be a single numeric value between 0 and 1")
  }

  if (!is.numeric(point_alpha) || length(point_alpha) != 1 || point_alpha < 0 || point_alpha > 1) {
    stop("'point_alpha' must be a single numeric value between 0 and 1")
  }

  # Process data and filter out NA values
  data_processed <- data %>%
    dplyr::filter(
      !is.na(!!rlang::sym(x)),
      !is.na(!!rlang::sym(y)),
      !is.na(!!rlang::sym(size)),
      !!rlang::sym(y) > 0  # Ensure positive p-values for log transformation
    )
  if (!is.null(color)) {
    data_processed <- data_processed %>%
      dplyr::filter(!is.na(!!rlang::sym(color)))
  }

  # Check if data remains after filtering
  if (nrow(data_processed) == 0) {
    stop("No valid data remaining after filtering NA values")
  }

  # Create summary statistics
  if (!is.null(color) && color %in% colnames(data_processed)) {
    data_summary <- data_processed %>%
      dplyr::group_by(!!rlang::sym(color)) %>%
      dplyr::summarise(
        count = dplyr::n(),
        mean_log2FC = mean(!!rlang::sym(x), na.rm = TRUE),
        median_pvalue = median(!!rlang::sym(y), na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::arrange(!!rlang::sym(color))
  } else {
    data_summary <- data.frame(
      total_genes = nrow(data_processed),
      mean_log2FC = mean(data_processed[[x]], na.rm = TRUE),
      median_pvalue = median(data_processed[[y]], na.rm = TRUE)
    )
  }

  # Create the volcano plot
  mapping <- if (is.null(color)) {
    ggplot2::aes(
      x = !!rlang::sym(x),
      y = -log10(!!rlang::sym(y)),
      size = !!rlang::sym(size)
    )
  } else {
    ggplot2::aes(
      x = !!rlang::sym(x),
      y = -log10(!!rlang::sym(y)),
      color = !!rlang::sym(color),
      size = !!rlang::sym(size)
    )
  }

  x_label <- if (x == "log2FoldChange") {
    expression(log[2]("Fold Change"))
  } else {
    x
  }
  y_label <- if (y == "padj") {
    expression(-log[10](italic(P)[adj]))
  } else if (y == "pvalue") {
    expression(-log[10](italic(P)))
  } else {
    paste0("-log10(", y, ")")
  }
  color_label <- if (is.null(color)) {
    NULL
  } else if (color == "group") {
    "Group"
  } else {
    color
  }

  plot_volcano <- data_processed %>%
    ggplot2::ggplot(mapping) +
    # Add threshold lines
    ggplot2::geom_hline(
      yintercept = -log10(yintercept),
      color = line_color,
      linetype = line_type,
      linewidth = line_size
    ) +
    ggplot2::geom_vline(
      xintercept = c(xintercept_abs, -xintercept_abs),
      color = line_color,
      linetype = line_type,
      linewidth = line_size
    ) +
    # Add points
    ggplot2::geom_point(alpha = point_alpha) +
    # Customize scales
    ggplot2::scale_size_continuous(
      range = point_size_range,
      guide = ggplot2::guide_legend(override.aes = list(alpha = 1))
    ) +
    # Add labels
    ggplot2::labs(
      x = x_label,
      y = y_label,
      title = title,
      subtitle = subtitle,
      size = ifelse(size == "baseMean", "Base Mean", size),
      color = color_label
    ) +
    # Apply custom theme (theme_bio is internal; call directly, not via ::)
    theme_bio() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 12),
      legend.position = "right"
    )

  if (!is.null(color) && identical(color, "regulation")) {
    reg_colors <- c(
      "Up" = "#D55E00",
      "Up-regulated" = "#D55E00",
      "Up_regulated" = "#D55E00",
      "Down" = "#56B4E9",
      "Down-regulated" = "#56B4E9",
      "Down_regulated" = "#56B4E9",
      "NS" = "#999999",
      "Not significant" = "#999999"
    )
    observed <- unique(as.character(data_processed$regulation))
    reg_colors <- reg_colors[intersect(names(reg_colors), observed)]
    plot_volcano <- plot_volcano +
      ggplot2::scale_color_manual(values = reg_colors, drop = FALSE)
  } else if (!is.null(color)) {
    plot_volcano <- plot_volcano + ggsci::scale_color_d3()
  }

  # Return comprehensive results
  list(
    plot_volcano = plot_volcano,
    data_summary = data_summary
  )
}
