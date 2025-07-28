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
#' @param color Character string specifying the column name used to map the 
#'   color of points (typically grouping variable). Default is "group".
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
#'   \item{plot.volcano}{A ggplot object that can be further customized 
#'     or directly displayed.}
#'   \item{data.summary}{A summary table showing the number of genes in 
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
#' @export
#'
#' @examples
#' library(dplyr)
#' library(bioRtools)
#'
#' # Load example data
#' data(df.rnaseq.plot_volcano)
#'
#' # Basic volcano plot
#' volcano_result <- plot_volcano(data = df.rnaseq.plot_volcano)
#' volcano_result$plot.volcano
#'
#' # View summary statistics
#' volcano_result$data.summary
#'
#' # Customized volcano plot with different thresholds
#' volcano_custom <- plot_volcano(
#'   data = df.rnaseq.plot_volcano,
#'   xintercept_abs = 1.5,
#'   yintercept = 0.01,
#'   title = "Differential Expression Analysis",
#'   subtitle = "Volcano Plot with Custom Thresholds"
#' )
#' volcano_custom$plot.volcano
#'
#' # Using different column names
#' plot_volcano(
#'   data = df.rnaseq.plot_volcano,
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
                         color = "group",
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
  
  # Check required columns
  required_cols <- c(x, y, size, color)
  missing_cols <- required_cols[!required_cols %in% colnames(data)]
  
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
      !is.na(!!rlang::sym(color)),
      !!rlang::sym(y) > 0  # Ensure positive p-values for log transformation
    )
  
  # Check if data remains after filtering
  if (nrow(data_processed) == 0) {
    stop("No valid data remaining after filtering NA values")
  }
  
  # Create summary statistics
  if (color %in% colnames(data_processed)) {
    data_summary <- data_processed %>%
      dplyr::group_by(!!rlang::sym(color)) %>%
      dplyr::summarise(
        count = dplyr::n(),
        mean_log2FC = mean(!!rlang::sym(x), na.rm = TRUE),
        median_pvalue = median(!!rlang::sym(y), na.rm = TRUE),
        .groups = 'drop'
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
  plot_volcano <- data_processed %>%
    ggplot2::ggplot(ggplot2::aes(
      x = !!rlang::sym(x),
      y = -log10(!!rlang::sym(y)),
      color = !!rlang::sym(color),
      size = !!rlang::sym(size)
    )) +
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
    ggsci::scale_color_d3() +
    # Add labels
    ggplot2::labs(
      x = ifelse(x == "log2FoldChange", 
                 expression(log[2]("Fold Change")), 
                 x),
      y = ifelse(y %in% c("padj", "pvalue"),
                 expression(-log[10]("p-value")),
                 paste0("-log10(", y, ")")),
      title = title,
      subtitle = subtitle,
      size = ifelse(size == "baseMean", "Base Mean", size),
      color = ifelse(color == "group", "Group", color)
    ) +
    # Apply custom theme
    bioRtools::theme_bio() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 12),
      legend.position = "right"
    )
  
  # Return comprehensive results following project conventions
  return(list(
    plot.volcano = plot_volcano,
    data.summary = data_summary
  ))
}