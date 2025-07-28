#' Plot Manhattan plot from GWAS results or genomic data
#'
#' @description
#' A comprehensive function for creating Manhattan plots from GWAS results or other
#' genome-wide data. Supports both data frames and file input, with customizable
#' visualization options.
#'
#' @param data A data frame containing the plotting data. If \code{input_file} is 
#'   provided, this parameter can be omitted.
#' @param input_file Character string specifying the path to the input data file. 
#'   If \code{data} is provided, this parameter can be omitted. Supports common 
#'   delimited formats (csv, tsv, txt).
#' @param chr_col Character string specifying the chromosome column name. 
#'   Default is "chr".
#' @param pos_col Character string specifying the physical position column name. 
#'   Default is "pos".
#' @param val_col Character string specifying the value column name (e.g., p-values). 
#'   Default is "value".
#' @param transform_log10 Logical value indicating whether to apply -log10 
#'   transformation to the values. Set to TRUE for p-values. Default is FALSE.
#' @param title Character string for the plot title. Default is "Manhattan Plot".
#' @param ylab Character string for the Y-axis label. If NULL, will be set 
#'   automatically based on \code{transform_log10}.
#' @param colors Character vector of colors for alternating chromosome coloring. 
#'   Default is c("gray50", "steelblue").
#' @param threshold_line Numeric value for drawing a horizontal threshold line. 
#'   Default is NULL (no line).
#' @param threshold_color Character string specifying the threshold line color. 
#'   Default is "red".
#' @param point_size Numeric value for point size. Default is 1.2.
#' @param point_alpha Numeric value for point transparency (0-1). Default is 0.8.
#'
#' @return A list containing three components:
#' \describe{
#'   \item{plot.manhattan}{A ggplot object that can be further customized or 
#'     directly displayed.}
#'   \item{data.processed}{The processed data frame used for plotting, including 
#'     cumulative positions and transformed values.}
#'   \item{chromosome.centers}{A data frame containing the center positions for 
#'     each chromosome, useful for custom axis labeling.}
#' }
#'
#' @details
#' The function automatically handles chromosome sorting using natural ordering 
#' (chr1, chr2, ..., chr10, chr11, ...) and calculates cumulative positions for 
#' proper Manhattan plot visualization. The plot uses alternating colors to 
#' distinguish chromosomes and includes optional threshold lines for significance 
#' cutoffs.
#'
#' @note
#' This function requires the following packages: dplyr, ggplot2, readr, and gtools.
#' The gtools package is specifically used for natural chromosome sorting.
#'
#' @export
#'
#' @examples
#' library(dplyr)
#' library(bioRtools)
#'
#' # Create mock GWAS data
#' set.seed(123)
#' n_chr <- 10
#' n_points_per_chr <- 1000
#' mock_data <- data.frame(
#'   chromosome = rep(paste0("chr", 1:n_chr), each = n_points_per_chr),
#'   position = unlist(lapply(1:n_chr, function(x) 
#'     sample(1:1e6, n_points_per_chr, replace = TRUE))),
#'   p_value = runif(n_chr * n_points_per_chr, min = 1e-8, max = 1)
#' )
#'
#' # Basic Manhattan plot with p-value transformation
#' manhattan_result <- plot_manhattan(
#'   data = mock_data,
#'   chr_col = "chromosome",
#'   pos_col = "position", 
#'   val_col = "p_value",
#'   transform_log10 = TRUE,
#'   title = "GWAS Manhattan Plot"
#' )
#'
#' # Display the plot
#' manhattan_result$plot.manhattan
#'
#' # Add significance threshold
#' manhattan_result2 <- plot_manhattan(
#'   data = mock_data,
#'   chr_col = "chromosome",
#'   pos_col = "position",
#'   val_col = "p_value", 
#'   transform_log10 = TRUE,
#'   title = "Manhattan Plot with Threshold",
#'   threshold_line = -log10(5e-8),
#'   colors = c("#276FBF", "#183059")
#' )
#'
#' # Access processed data for further analysis
#' head(manhattan_result$data.processed)
#' manhattan_result$chromosome.centers
#'
plot_manhattan <- function(data = NULL,
                           input_file = NULL,
                           chr_col = "chr",
                           pos_col = "pos", 
                           val_col = "value",
                           transform_log10 = FALSE,
                           title = "Manhattan Plot",
                           ylab = NULL,
                           colors = c("gray50", "steelblue"),
                           threshold_line = NULL,
                           threshold_color = "red",
                           point_size = 1.2,
                           point_alpha = 0.8) {
  
  # Input validation
  if (is.null(data) && is.null(input_file)) {
    stop("Either 'data' or 'input_file' must be provided")
  }
  
  if (!is.null(data) && !is.null(input_file)) {
    warning("Both 'data' and 'input_file' provided. Using 'data' and ignoring 'input_file'")
  }
  
  # Check required packages
  required_packages <- c("dplyr", "ggplot2", "gtools", "rlang")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  
  if (length(missing_packages) > 0) {
    stop("Required packages missing: ", paste(missing_packages, collapse = ", "))
  }
  
  # Data loading and preparation
  if (!is.null(data)) {
    df <- as.data.frame(data)
  } else {
    if (!file.exists(input_file)) {
      stop("Input file not found: ", input_file)
    }
    
    tryCatch({
      df <- readr::read_delim(input_file, show_col_types = FALSE)
    }, error = function(e) {
      stop("Failed to read input file: ", e$message)
    })
  }
  
  # Validate required columns
  required_cols <- c(chr_col, pos_col, val_col)
  missing_cols <- required_cols[!required_cols %in% colnames(df)]
  
  if (length(missing_cols) > 0) {
    stop("Required columns missing from data: ", paste(missing_cols, collapse = ", "))
  }
  
  # Process data
  df_processed <- df %>%
    dplyr::select(
      chr = !!rlang::sym(chr_col),
      pos = !!rlang::sym(pos_col), 
      value = !!rlang::sym(val_col)
    ) %>%
    dplyr::mutate(
      pos = as.numeric(pos),
      value = as.numeric(value)
    ) %>%
    dplyr::filter(
      !is.na(pos),
      !is.na(value),
      pos > 0,
      value > 0  # Ensure positive values for log transformation
    )
  
  # Check if data remains after filtering
  if (nrow(df_processed) == 0) {
    stop("No valid data remaining after filtering")
  }
  
  # Apply transformation if requested
  if (transform_log10) {
    df_processed <- df_processed %>%
      dplyr::mutate(value = -log10(value))
    
    if (is.null(ylab)) {
      ylab <- expression(-log[10]("p-value"))
    }
  } else {
    if (is.null(ylab)) {
      ylab <- val_col
    }
  }
  
  # Calculate cumulative chromosome positions
  data_cum <- df_processed %>%
    dplyr::group_by(chr) %>%
    dplyr::summarise(max_pos = max(pos, na.rm = TRUE), .groups = 'drop') %>%
    dplyr::mutate(
      chr = factor(chr, levels = gtools::mixedsort(unique(chr)))
    ) %>%
    dplyr::arrange(chr) %>%
    dplyr::mutate(
      cum_pos_end = cumsum(max_pos),
      cum_pos_start = dplyr::lag(cum_pos_end, default = 0)
    ) %>%
    dplyr::select(chr, cum_pos_start) %>%
    dplyr::left_join(df_processed, ., by = "chr") %>%
    dplyr::arrange(chr, pos) %>%
    dplyr::mutate(
      pos_cum = pos + cum_pos_start,
      chr_factor = factor(chr, levels = gtools::mixedsort(unique(chr)))
    )
  
  # Calculate chromosome center positions for x-axis labels
  chromosome_centers <- data_cum %>%
    dplyr::group_by(chr_factor) %>%
    dplyr::summarise(
      center = mean(pos_cum, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    dplyr::arrange(chr_factor)
  
  # Create the Manhattan plot
  plot_manhattan <- data_cum %>%
    ggplot2::ggplot(ggplot2::aes(x = pos_cum, y = value)) +
    ggplot2::geom_point(
      ggplot2::aes(color = chr_factor),
      alpha = point_alpha,
      size = point_size
    ) +
    ggplot2::scale_color_manual(
      values = rep(colors, length.out = nlevels(data_cum$chr_factor))
    ) +
    ggplot2::scale_x_continuous(
      breaks = chromosome_centers$center,
      labels = chromosome_centers$chr_factor,
      expand = ggplot2::expansion(mult = c(0.01, 0.01))
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.1))
    ) +
    ggplot2::labs(
      x = "Chromosome",
      y = ylab,
      title = title
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = 45, 
        hjust = 1, 
        size = 9
      ),
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = 14,
        face = "bold"
      )
    )
  
  # Add threshold line if specified
  if (!is.null(threshold_line)) {
    plot_manhattan <- plot_manhattan +
      ggplot2::geom_hline(
        yintercept = threshold_line,
        color = threshold_color,
        linetype = "dashed",
        linewidth = 0.8
      )
  }
  
  # Return comprehensive results following project conventions
  return(list(
    plot.manhattan = plot_manhattan,
    data.processed = data_cum,
    chromosome.centers = chromosome_centers
  ))
}