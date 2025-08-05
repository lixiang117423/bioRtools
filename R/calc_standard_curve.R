#' Standard Curve Calculation for qPCR
#'
#' @description
#' This function calculates the standard curve and amplification efficiency
#' for qPCR primers. Based on the amplification efficiency, it determines
#' which method can be used to calculate expression levels.
#'
#' @param cq_table A data frame containing position and Cq values
#' @param concentration_table A data frame containing position and concentration values
#' @param highest_concentration The highest concentration for calculation
#' @param lowest_concentration The lowest concentration for calculation
#' @param dilution_factor Dilution factor of cDNA template (default: 4)
#' @param use_mean_cq Whether to calculate using mean Cq values (default: TRUE)
#'
#' @return A list containing:
#'   \item{table}{Data frame with regression statistics and efficiency}
#'   \item{figure}{ggplot object showing the standard curve}
#'
#' @importFrom dplyr left_join filter group_by mutate ungroup summarise
#' @importFrom stats lm
#' @importFrom broom glance
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth labs theme_bw
#' @importFrom ggpmisc stat_poly_eq
#'
#' @export
#'
#' @examples
#' # Load example data
#' cq_data_path <- system.file("examples", "calsc.cq.txt", package = "qPCRtools")
#' concentration_data_path <- system.file("examples", "calsc.info.txt", package = "qPCRtools")
#' cq_data <- read.table(cq_data_path, header = TRUE)
#' concentration_data <- read.table(concentration_data_path, header = TRUE)
#' 
#' # Calculate standard curve
#' result <- calc_standard_curve(
#'   cq_table = cq_data,
#'   concentration_table = concentration_data,
#'   lowest_concentration = 4,
#'   highest_concentration = 4096,
#'   dilution_factor = 4,
#'   use_mean_cq = TRUE
#' )
#' 
#' # View results
#' result$table
#' result$figure
#'
#' @author Xiang LI <lixiang117423@gmail.com>

# Remove global variables declaration as it's not best practice
# Instead, use proper NSE handling

calc_standard_curve <- function(cq_table,
                               concentration_table,
                               highest_concentration,
                               lowest_concentration,
                               dilution_factor = 4,
                               use_mean_cq = TRUE) {
  
  # Input validation
  validate_inputs(cq_table, concentration_table, highest_concentration, 
                 lowest_concentration, dilution_factor, use_mean_cq)
  
  # Prepare data
  processed_data <- prepare_data(cq_table, concentration_table, 
                                highest_concentration, lowest_concentration, 
                                dilution_factor, use_mean_cq)
  
  # Calculate regression statistics
  regression_results <- calculate_regression_stats(processed_data, dilution_factor, use_mean_cq)
  
  # Create plot
  plot_result <- create_standard_curve_plot(processed_data, use_mean_cq)
  
  return(list(table = regression_results, figure = plot_result))
}

#' Validate input parameters
#' @keywords internal
validate_inputs <- function(cq_table, concentration_table, highest_concentration, 
                           lowest_concentration, dilution_factor, use_mean_cq) {
  
  if (!is.data.frame(cq_table)) {
    stop("cq_table must be a data frame")
  }
  
  if (!is.data.frame(concentration_table)) {
    stop("concentration_table must be a data frame")
  }
  
  if (!all(c("Position", "Cq", "Gene") %in% names(cq_table))) {
    stop("cq_table must contain columns: Position, Cq, Gene")
  }
  
  if (!all(c("Position", "Conc") %in% names(concentration_table))) {
    stop("concentration_table must contain columns: Position, Conc")
  }
  
  if (highest_concentration <= lowest_concentration) {
    stop("highest_concentration must be greater than lowest_concentration")
  }
  
  if (dilution_factor <= 1) {
    stop("dilution_factor must be greater than 1")
  }
  
  if (!is.logical(use_mean_cq)) {
    stop("use_mean_cq must be TRUE or FALSE")
  }
}

#' Prepare and filter data for analysis
#' @keywords internal
prepare_data <- function(cq_table, concentration_table, highest_concentration, 
                        lowest_concentration, dilution_factor, use_mean_cq) {
  
  # Join tables and filter concentration range
  combined_data <- cq_table %>%
    dplyr::left_join(concentration_table, by = "Position") %>%
    dplyr::filter(.data$Conc >= lowest_concentration & .data$Conc <= highest_concentration)
  
  # Calculate mean Cq if requested and log-transform concentration
  if (use_mean_cq) {
    processed_data <- combined_data %>%
      dplyr::group_by(.data$Gene, .data$Conc) %>%
      dplyr::summarise(
        mean_cq = mean(.data$Cq, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(log_conc = log(.data$Conc, base = dilution_factor)) %>%
      dplyr::group_by(.data$Gene) %>%
      dplyr::mutate(
        max_cq = max(.data$mean_cq, na.rm = TRUE),
        min_cq = min(.data$mean_cq, na.rm = TRUE)
      ) %>%
      dplyr::ungroup()
  } else {
    processed_data <- combined_data %>%
      dplyr::mutate(log_conc = log(.data$Conc, base = dilution_factor)) %>%
      dplyr::group_by(.data$Gene) %>%
      dplyr::mutate(
        max_cq = max(.data$Cq, na.rm = TRUE),
        min_cq = min(.data$Cq, na.rm = TRUE)
      ) %>%
      dplyr::ungroup()
  }
  
  return(processed_data)
}

#' Calculate regression statistics for each gene
#' @keywords internal
calculate_regression_stats <- function(processed_data, dilution_factor, use_mean_cq) {
  
  genes <- unique(processed_data$Gene)
  results_list <- vector("list", length(genes))
  
  for (i in seq_along(genes)) {
    gene_data <- processed_data %>%
      dplyr::filter(.data$Gene == genes[i])
    
    # Fit linear model
    if (use_mean_cq) {
      fit <- stats::lm(mean_cq ~ log_conc, data = gene_data)
      max_cq <- unique(gene_data$max_cq)
      min_cq <- unique(gene_data$min_cq)
    } else {
      fit <- stats::lm(Cq ~ log_conc, data = gene_data)
      max_cq <- unique(gene_data$max_cq)
      min_cq <- unique(gene_data$min_cq)
    }
    
    # Extract model parameters
    coefficients <- stats::coef(fit)
    intercept <- round(coefficients[1], 2)
    slope <- round(coefficients[2], 2)
    
    # Get model statistics
    model_stats <- broom::glance(fit)
    r_squared <- round(model_stats$r.squared, 4)
    p_value <- round(model_stats$p.value, 5)
    
    # Calculate amplification efficiency
    efficiency <- round(dilution_factor^(-1/slope) - 1, 3)
    
    # Create result row
    results_list[[i]] <- data.frame(
      Gene = genes[i],
      Formula = sprintf("y = %.2f*log_conc + %.2f", slope, intercept),
      Slope = slope,
      Intercept = intercept,
      R2 = r_squared,
      P_value = p_value,
      Max_Cq = max_cq,
      Min_Cq = min_cq,
      Efficiency = efficiency,
      Date = as.character(Sys.Date()),
      stringsAsFactors = FALSE
    )
  }
  
  return(do.call(rbind, results_list))
}

#' Create standard curve plot
#' @keywords internal
create_standard_curve_plot <- function(processed_data, use_mean_cq) {
  
  # Determine y-axis variable
  y_var <- if (use_mean_cq) "mean_cq" else "Cq"
  y_label <- if (use_mean_cq) "Mean Cq Value" else "Cq Value"
  
  # Calculate label positions
  n_genes <- length(unique(processed_data$Gene))
  label_y_positions <- seq(0.05, 0.05 + 0.06 * n_genes, 0.06)
  
  # Create plot
  p <- ggplot2::ggplot(processed_data, ggplot2::aes(.data$log_conc, .data[[y_var]], color = .data$Gene)) +
    ggplot2::geom_point(size = 2, alpha = 0.7) +
    ggplot2::geom_smooth(method = "lm", se = TRUE, show.legend = FALSE) +
    ggpmisc::stat_poly_eq(
      ggplot2::aes(
        label = paste(
          ggplot2::after_stat(rr.label),
          ggplot2::after_stat(p.value.label),
          sep = "~~~~"
        )
      ),
      show.legend = FALSE,
      formula = y ~ x,
      parse = TRUE,
      rr.digits = 4,
      coef.digits = 3,
      label.x = 0.05,
      label.y = label_y_positions[seq_len(n_genes)]
    ) +
    ggplot2::labs(
      x = "Log Relative Concentration",
      y = y_label,
      title = "qPCR Standard Curve",
      color = "Gene"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14),
      legend.position = "bottom"
    )
  
  return(p)
}

# For backward compatibility, create an alias with the original name
CalCurve <- calc_standard_curve