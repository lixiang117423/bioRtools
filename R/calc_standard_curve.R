#' Standard Curve Calculation for qPCR
#'
#' @description
#' This function calculates the standard curve and amplification efficiency
#' for qPCR primers. Based on the amplification efficiency, it determines
#' which method can be used to calculate expression levels.
#'
#' @param cq_table A data frame containing position and Cq values.
#'   Must contain columns: Position, Cq
#' @param concentration_table A data frame containing position and concentration values.
#'   Must contain columns: Position, Conc. May also contain Gene column.
#' @param gene_name Optional character string specifying gene name. If NULL and Gene
#'   column exists in concentration_table, will use that. If NULL and no Gene column,
#'   will use "Unknown_Gene" (default: NULL)
#' @param highest_concentration The highest concentration for calculation
#' @param lowest_concentration The lowest concentration for calculation
#' @param dilution_factor Dilution factor of cDNA template (default: 4)
#' @param use_mean_cq Whether to calculate using mean Cq values (default: TRUE)
#'
#' @return A list containing:
#'   \item{table}{Data frame with regression statistics and efficiency}
#'   \item{plot}{ggplot object showing the standard curve}
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
#' cq_data_path <- system.file("extdata/qPCR", "calsc.cq.txt", package = "bioRtools")
#' concentration_data_path <- system.file("extdata/qPCR", "calsc.info.txt", package = "bioRtools")
#' cq_data <- read.table(cq_data_path, header = TRUE)
#' concentration_data <- read.table(concentration_data_path, header = TRUE)
#'
#' # Calculate standard curve (gene info in concentration_table)
#' result1 <- calc_standard_curve(
#'   cq_table = cq_data,
#'   concentration_table = concentration_data,
#'   lowest_concentration = 4,
#'   highest_concentration = 4096,
#'   dilution_factor = 4,
#'   use_mean_cq = TRUE
#' )
#'
#' # Calculate standard curve (specify gene name manually)
#' result2 <- calc_standard_curve(
#'   cq_table = cq_data,
#'   concentration_table = concentration_data,
#'   gene_name = "MyGene",
#'   lowest_concentration = 4,
#'   highest_concentration = 4096,
#'   dilution_factor = 4,
#'   use_mean_cq = TRUE
#' )
#'
#' # View results
#' result1$table
#' result1$plot
#'
#' @author Xiang LI <lixiang117423@gmail.com>

calc_standard_curve <- function(cq_table,
                                concentration_table,
                                gene_name = NULL,
                                highest_concentration,
                                lowest_concentration,
                                dilution_factor = 4,
                                use_mean_cq = TRUE) {
  # Merge data first to understand the structure
  merged_data <- merge_standard_curve_data(cq_table, concentration_table, gene_name)

  # Input validation (after merging)
  validate_inputs_fixed(merged_data, highest_concentration,
    lowest_concentration, dilution_factor, use_mean_cq)

  # Prepare data
  processed_data <- prepare_data_fixed(merged_data, highest_concentration,
    lowest_concentration, dilution_factor, use_mean_cq)

  # Calculate regression statistics
  regression_results <- calculate_regression_stats(processed_data, dilution_factor, use_mean_cq)

  # Create plot
  plot_result <- create_standard_curve_plot(processed_data, use_mean_cq)

  return(list(table = regression_results, plot = plot_result))
}

#' Merge data and handle gene information flexibly
#' @keywords internal
merge_standard_curve_data <- function(cq_table, concentration_table, gene_name) {
  # Basic merge
  merged_data <- cq_table %>%
    dplyr::left_join(concentration_table, by = "Position")

  # Handle gene information flexibly
  if (!is.null(gene_name)) {
    # Use provided gene name
    merged_data$Gene <- gene_name
  } else if ("Gene" %in% names(merged_data)) {
    # Gene column already exists (likely from concentration_table)
    # Keep as is
  } else {
    # No gene information available, create a default
    merged_data$Gene <- "Unknown_Gene"
    message("No gene information provided. Using 'Unknown_Gene' as gene name.")
  }

  return(merged_data)
}

#' Validate input parameters (after merging)
#' @keywords internal
validate_inputs_fixed <- function(merged_data, highest_concentration,
                                  lowest_concentration, dilution_factor, use_mean_cq) {

  if (nrow(merged_data) == 0) {
    stop("No data available after merging tables. Check that Position columns match.")
  }

  # Check required columns in merged data
  required_cols <- c("Position", "Cq", "Conc", "Gene")
  missing_cols <- setdiff(required_cols, names(merged_data))
  if (length(missing_cols) > 0) {
    stop(sprintf("Merged data missing required columns: %s",
      paste(missing_cols, collapse = ", ")))
  }

  # Check for data quality
  if (all(is.na(merged_data$Cq))) {
    stop("All Cq values are missing")
  }

  if (all(is.na(merged_data$Conc))) {
    stop("All concentration values are missing")
  }

  # Check concentration range
  available_conc <- unique(merged_data$Conc[!is.na(merged_data$Conc)])
  if (length(available_conc) < 2) {
    stop("At least 2 different concentration values are required for standard curve")
  }

  conc_in_range <- available_conc[available_conc >= lowest_concentration &
    available_conc <= highest_concentration]
  if (length(conc_in_range) < 2) {
    stop(sprintf("At least 2 concentration values must be within range [%s, %s]. Available: %s",
      lowest_concentration, highest_concentration, paste(available_conc, collapse = ", ")))
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
prepare_data_fixed <- function(merged_data, highest_concentration,
                               lowest_concentration, dilution_factor, use_mean_cq) {
  # Filter concentration range and remove missing values
  filtered_data <- merged_data %>%
    dplyr::filter(
      !is.na(.data$Cq),
      !is.na(.data$Conc),
      .data$Conc >= lowest_concentration & .data$Conc <= highest_concentration
    )

  # Calculate mean Cq if requested and log-transform concentration
  if (use_mean_cq) {
    processed_data <- filtered_data %>%
      dplyr::group_by(.data$Gene, .data$Conc) %>%
      dplyr::summarise(
        mean_cq = mean(.data$Cq, na.rm = TRUE),
        n_replicates = dplyr::n(),
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
    processed_data <- filtered_data %>%
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

    # Check if we have enough data points
    if (nrow(gene_data) < 2) {
      warning(sprintf("Gene %s has fewer than 2 data points. Skipping.", genes[i]))
      next
    }

    tryCatch(
      {
        # Fit linear model
        if (use_mean_cq) {
          fit <- stats::lm(mean_cq ~ log_conc, data = gene_data)
          max_cq <- unique(gene_data$max_cq)[1]
          min_cq <- unique(gene_data$min_cq)[1]
        } else {
          fit <- stats::lm(Cq ~ log_conc, data = gene_data)
          max_cq <- unique(gene_data$max_cq)[1]
          min_cq <- unique(gene_data$min_cq)[1]
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
        efficiency <- round(dilution_factor^(-1 / slope) - 1, 3)

        # Determine efficiency quality
        efficiency_quality <- dplyr::case_when(
          efficiency >= 0.90 & efficiency <= 1.10 ~ "Excellent (90-110%)",
          efficiency >= 0.80 & efficiency <= 1.20 ~ "Good (80-120%)",
          efficiency >= 0.70 & efficiency <= 1.30 ~ "Acceptable (70-130%)",
          TRUE ~ "Poor (<70% or >130%)"
        )

        # Determine which qPCR method can be used
        method_recommendation <- dplyr::case_when(
          efficiency >= 0.90 & efficiency <= 1.10 ~ "2^-ΔΔCt method (relative quantification)",
          efficiency >= 0.80 & efficiency <= 1.20 ~ "Efficiency correction method recommended",
          TRUE ~ "Standard curve method required"
        )

        # Create result row
        results_list[[i]] <- data.frame(
          Gene = genes[i],
          Formula = sprintf("y = %.2f*log_conc + %.2f", slope, intercept),
          Slope = slope,
          Intercept = intercept,
          R2 = r_squared,
          P_value = p_value,
          Max_Cq = round(max_cq, 2),
          Min_Cq = round(min_cq, 2),
          Efficiency = efficiency,
          Efficiency_Percent = paste0(round(efficiency * 100, 1), "%"),
          Efficiency_Quality = efficiency_quality,
          Method_Recommendation = method_recommendation,
          Date = as.character(Sys.Date()),
          stringsAsFactors = FALSE
        )

      },
      error = function(e) {
        warning(sprintf("Failed to calculate regression for gene %s: %s", genes[i], e$message))

        # Return a row with NA values
        results_list[[i]] <- data.frame(
          Gene = genes[i],
          Formula = "Failed",
          Slope = NA,
          Intercept = NA,
          R2 = NA,
          P_value = NA,
          Max_Cq = NA,
          Min_Cq = NA,
          Efficiency = NA,
          Efficiency_Percent = "N/A",
          Efficiency_Quality = "Failed",
          Method_Recommendation = "Analysis failed",
          Date = as.character(Sys.Date()),
          stringsAsFactors = FALSE
        )
      })
  }

  # Remove NULL entries (genes that were skipped)
  results_list <- results_list[!sapply(results_list, is.null)]

  if (length(results_list) == 0) {
    stop("No valid results could be calculated for any gene")
  }

  return(do.call(rbind, results_list))
}

#' Create standard curve plot
#' @keywords internal
create_standard_curve_plot <- function(processed_data, use_mean_cq) {

  if (nrow(processed_data) == 0) {
    warning("No data available for plotting")
    return(NULL)
  }

  # Determine y-axis variable
  y_var <- if (use_mean_cq) "mean_cq" else "Cq"
  y_label <- if (use_mean_cq) "Mean Cq Value" else "Cq Value"

  # Calculate label positions for multiple genes
  n_genes <- length(unique(processed_data$Gene))
  label_y_positions <- seq(0.05, 0.05 + 0.06 * (n_genes - 1), length.out = n_genes)

  # Create base plot
  p <- ggplot2::ggplot(processed_data, ggplot2::aes(.data$log_conc, .data[[y_var]])) +
    ggplot2::geom_point(ggplot2::aes(color = .data$Gene), size = 2, alpha = 0.7) +
    ggplot2::geom_smooth(ggplot2::aes(color = .data$Gene), method = "lm", se = TRUE, show.legend = FALSE)

  # Add regression equation and R² for each gene
  if (requireNamespace("ggpmisc", quietly = TRUE)) {
    p <- p + ggpmisc::stat_poly_eq(
      ggplot2::aes(
        color = .data$Gene,
        label = paste(
          ggplot2::after_stat(eq.label),
          ggplot2::after_stat(rr.label),
          sep = "~~~~"
        )
      ),
      show.legend = FALSE,
      formula = y ~ x,
      parse = TRUE,
      rr.digits = 4,
      coef.digits = 3,
      label.x = 0.05,
      label.y = if (n_genes == 1) 0.95 else label_y_positions
    )
  }

  # Finish plot
  p <- p +
    ggplot2::labs(
      x = sprintf("Log%d Relative Concentration", unique(round(log(4, base = 4)))),
      y = y_label,
      title = "qPCR Standard Curve Analysis",
      subtitle = paste("Points:", nrow(processed_data), "| Genes:", n_genes),
      color = "Gene"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10),
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    )

  # If only one gene, remove color legend
  if (n_genes == 1) {
    p <- p + ggplot2::theme(legend.position = "none")
  }

  return(p)
}

# For backward compatibility, create an alias with the original name
CalCurve <- function(cq.table, concentration.table, highest.concentration,
                     lowest.concentration, dilution.factor = 4, use.mean.cq = TRUE) {

  result <- calc_standard_curve(
    cq_table = cq.table,
    concentration_table = concentration.table,
    gene_name = NULL,
    highest_concentration = highest.concentration,
    lowest_concentration = lowest.concentration,
    dilution_factor = dilution.factor,
    use_mean_cq = use.mean.cq
  )

  return(result)
}

# Output column explanations:
#
# table:
# - Gene: Gene name (from concentration_table, provided gene_name, or "Unknown_Gene")
# - Formula: Linear regression equation (y = slope*log_conc + intercept)
# - Slope: Slope of the regression line
# - Intercept: Y-intercept of the regression line
# - R2: R-squared value (coefficient of determination)
# - P_value: P-value of the regression model
# - Max_Cq: Maximum Cq value in the dataset for this gene
# - Min_Cq: Minimum Cq value in the dataset for this gene
# - Efficiency: Amplification efficiency (calculated as dilution_factor^(-1/slope) - 1)
# - Efficiency_Percent: Efficiency expressed as percentage
# - Efficiency_Quality: Quality assessment of efficiency (Excellent/Good/Acceptable/Poor)
# - Method_Recommendation: Recommended qPCR quantification method based on efficiency
# - Date: Date of analysis
#
# Efficiency interpretation:
# - 90-110%: Excellent - Can use 2^-ΔΔCt method
# - 80-120%: Good - Efficiency correction recommended
# - 70-130%: Acceptable - Standard curve method preferred
# - <70% or >130%: Poor - Check primer design and PCR conditions
#
# Method recommendations:
# - 2^-ΔΔCt method: For high efficiency primers (90-110%)
# - Efficiency correction method: For good efficiency primers (80-120%)
# - Standard curve method: Required for poor efficiency primers
