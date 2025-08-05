#' Calculate Gene Expression Using Delta-Ct Method
#'
#' @description
#' This function calculates gene expression levels using the delta-Ct method,
#' which normalizes target gene expression to a reference gene.
#'
#' @param cq_table A data frame containing position, Cq values, and gene information.
#'   Must contain columns: Position, Cq, Gene
#' @param design_table A data frame containing experimental design information.
#'   Must contain columns: Position, Group, BioRep
#' @param reference_gene Character string specifying the name of the reference gene
#'   (default: "Actin")
#' @param create_plot Logical indicating whether to create a visualization plot
#'   (default: TRUE)
#'
#' @return A list containing:
#'   \item{expression_data}{Data frame with calculated expression values and statistics}
#'   \item{summary_table}{Summary statistics by group and gene}
#'   \item{plot}{ggplot object showing expression levels (if create_plot = TRUE)}
#'
#' @importFrom dplyr left_join filter group_by mutate ungroup select distinct
#'   rename summarise n
#' @importFrom ggplot2 ggplot aes geom_col geom_errorbar labs theme_bw
#'   facet_wrap position_dodge
#' @importFrom stats sd
#'
#' @export
#'
#' @examples
#' # Load example data
#' cq_data_path <- system.file("examples", "dct.cq.txt", package = "qPCRtools")
#' design_data_path <- system.file("examples", "dct.design.txt", package = "qPCRtools")
#' cq_data <- read.table(cq_data_path, sep = ",", header = TRUE)
#' design_data <- read.table(design_data_path, sep = ",", header = TRUE)
#' 
#' # Calculate expression using delta-Ct method
#' result <- calc_expression_delta_ct(
#'   cq_table = cq_data,
#'   design_table = design_data,
#'   reference_gene = "Actin"
#' )
#' 
#' # View results
#' result$summary_table
#' result$plot
#'
#' @author Xiang LI <lixiang117423@gmail.com>

calc_expression_delta_ct <- function(cq_table,
                                   design_table,
                                   reference_gene = "Actin",
                                   create_plot = TRUE) {
  
  # Input validation
  validate_delta_ct_inputs(cq_table, design_table, reference_gene)
  
  # Merge and prepare data
  merged_data <- prepare_delta_ct_data(cq_table, design_table)
  
  # Calculate reference gene means
  reference_data <- calculate_reference_means(merged_data, reference_gene)
  
  # Calculate target gene expression
  expression_data <- calculate_target_expression(merged_data, reference_data, reference_gene)
  
  # Calculate summary statistics
  summary_table <- calculate_expression_summary(expression_data)
  
  # Create plot if requested
  plot_result <- NULL
  if (create_plot) {
    plot_result <- create_expression_plot(summary_table)
  }
  
  return(list(
    expression_data = expression_data,
    summary_table = summary_table,
    plot = plot_result
  ))
}

#' Validate input parameters for delta-Ct calculation
#' @keywords internal
validate_delta_ct_inputs <- function(cq_table, design_table, reference_gene) {
  
  if (!is.data.frame(cq_table)) {
    stop("cq_table must be a data frame")
  }
  
  if (!is.data.frame(design_table)) {
    stop("design_table must be a data frame")
  }
  
  # Check required columns in cq_table
  required_cq_cols <- c("Position", "Cq", "Gene")
  missing_cq_cols <- setdiff(required_cq_cols, names(cq_table))
  if (length(missing_cq_cols) > 0) {
    stop(sprintf("cq_table missing required columns: %s", 
                paste(missing_cq_cols, collapse = ", ")))
  }
  
  # Check required columns in design_table
  required_design_cols <- c("Position", "Group", "BioRep")
  missing_design_cols <- setdiff(required_design_cols, names(design_table))
  if (length(missing_design_cols) > 0) {
    stop(sprintf("design_table missing required columns: %s", 
                paste(missing_design_cols, collapse = ", ")))
  }
  
  # Check if reference gene exists in data
  if (!reference_gene %in% cq_table$Gene) {
    stop(sprintf("Reference gene '%s' not found in cq_table", reference_gene))
  }
  
  # Check for matching positions
  common_positions <- intersect(cq_table$Position, design_table$Position)
  if (length(common_positions) == 0) {
    stop("No matching positions found between cq_table and design_table")
  }
}

#' Prepare and merge data for delta-Ct calculation
#' @keywords internal
prepare_delta_ct_data <- function(cq_table, design_table) {
  
  merged_data <- cq_table %>%
    dplyr::left_join(design_table, by = "Position") %>%
    dplyr::rename(
      position = .data$Position,
      cq = .data$Cq,
      group = .data$Group,
      gene = .data$Gene,
      bio_rep = .data$BioRep
    ) %>%
    dplyr::filter(!is.na(.data$group), !is.na(.data$bio_rep))
  
  return(merged_data)
}

#' Calculate reference gene mean Cq values
#' @keywords internal
calculate_reference_means <- function(merged_data, reference_gene) {
  
  reference_means <- merged_data %>%
    dplyr::filter(.data$gene == reference_gene) %>%
    dplyr::group_by(.data$group, .data$bio_rep) %>%
    dplyr::summarise(
      mean_ref_cq = mean(.data$cq, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(group_biorep = paste0(.data$group, "_", .data$bio_rep))
  
  return(reference_means)
}

#' Calculate target gene expression levels
#' @keywords internal
calculate_target_expression <- function(merged_data, reference_data, reference_gene) {
  
  expression_data <- merged_data %>%
    dplyr::filter(.data$gene != reference_gene) %>%
    dplyr::mutate(group_biorep = paste0(.data$group, "_", .data$bio_rep)) %>%
    dplyr::left_join(
      reference_data %>% dplyr::select(.data$group_biorep, .data$mean_ref_cq),
      by = "group_biorep"
    ) %>%
    dplyr::filter(!is.na(.data$mean_ref_cq)) %>%
    dplyr::mutate(
      delta_ct = .data$cq - .data$mean_ref_cq,
      relative_expression = 2^(-(.data$delta_ct))
    ) %>%
    dplyr::select(-.data$group_biorep)
  
  return(expression_data)
}

#' Calculate summary statistics for expression data
#' @keywords internal
calculate_expression_summary <- function(expression_data) {
  
  summary_stats <- expression_data %>%
    dplyr::group_by(.data$group, .data$gene) %>%
    dplyr::summarise(
      n_replicates = dplyr::n(),
      mean_expression = mean(.data$relative_expression, na.rm = TRUE),
      sd_expression = stats::sd(.data$relative_expression, na.rm = TRUE),
      se_expression = .data$sd_expression / sqrt(.data$n_replicates),
      mean_delta_ct = mean(.data$delta_ct, na.rm = TRUE),
      sd_delta_ct = stats::sd(.data$delta_ct, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      sd_expression = ifelse(is.na(.data$sd_expression), 0, .data$sd_expression),
      se_expression = ifelse(is.na(.data$se_expression), 0, .data$se_expression),
      sd_delta_ct = ifelse(is.na(.data$sd_delta_ct), 0, .data$sd_delta_ct)
    )
  
  return(summary_stats)
}

#' Create expression level plot
#' @keywords internal
create_expression_plot <- function(summary_table) {
  
  if (nrow(summary_table) == 0) {
    warning("No data available for plotting")
    return(NULL)
  }
  
  # Calculate error bar limits
  summary_table <- summary_table %>%
    dplyr::mutate(
      error_min = pmax(0, .data$mean_expression - .data$se_expression),
      error_max = .data$mean_expression + .data$se_expression
    )
  
  # Create plot
  plot_result <- ggplot2::ggplot(
    summary_table, 
    ggplot2::aes(x = .data$group, y = .data$mean_expression, fill = .data$group)
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.8),
      alpha = 0.7,
      color = "black",
      size = 0.3
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$error_min, ymax = .data$error_max),
      position = ggplot2::position_dodge(width = 0.8),
      width = 0.25,
      size = 0.5
    ) +
    ggplot2::facet_wrap(~ .data$gene, scales = "free_y", ncol = 3) +
    ggplot2::labs(
      x = "Treatment Group",
      y = "Relative Expression (2^-ΔCt)",
      title = "Gene Expression Analysis (Delta-Ct Method)",
      subtitle = "Error bars represent standard error",
      fill = "Group"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      strip.background = ggplot2::element_rect(fill = "lightgray"),
      legend.position = "bottom"
    )
  
  return(plot_result)
}

# For backward compatibility
CalExp2dCt <- function(cq.table, design.table, ref.gene = "Actin") {
  result <- calc_expression_delta_ct(
    cq_table = cq.table,
    design_table = design.table,
    reference_gene = ref.gene,
    create_plot = FALSE
  )
  
  # Return only expression data to match original function behavior
  return(result$expression_data)
}