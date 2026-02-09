#' Calculate Gene Expression Using Delta-Ct Method
#'
#' @description
#' This function calculates gene expression levels using the delta-Ct method,
#' which normalizes target gene expression to a reference gene.
#'
#' @param cq_table A data frame containing position, Cq values, and gene information.
#'   Should contain columns: Position, Cq, Gene
#' @param design_table A data frame containing experimental design information.
#'   Should contain columns: Position, Group, BioRep
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
#' cq_data_path <- system.file("extdata/qPCR", "dct.cq.txt", package = "bioRtools")
#' design_data_path <- system.file("extdata/qPCR", "dct.design.txt", package = "bioRtools")
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
  # Step 1: Merge and prepare data (similar to original approach)
  merged_data <- prepare_delta_ct_data(cq_table, design_table)

  # Step 2: Calculate reference gene means
  reference_data <- calculate_reference_means(merged_data, reference_gene)

  # Step 3: Calculate target gene expression
  expression_data <- calculate_target_expression(merged_data, reference_data, reference_gene)

  # Step 4: Calculate summary statistics
  summary_table <- calculate_expression_summary(expression_data)

  # Step 5: Create plot if requested
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

#' Prepare and merge data for delta-Ct calculation
#' @keywords internal
prepare_delta_ct_data <- function(cq_table, design_table) {
  # Use the same approach as original function - no strict validation
  merged_data <- cq_table %>%
    dplyr::left_join(design_table, by = "Position") %>%
    dplyr::rename(
      position = Position,
      cq = Cq,
      group = Group,
      gene = Gene,
      bio_rep = BioRep
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
    dplyr::mutate(group_biorep = paste0(.data$group, .data$bio_rep)) %>%
    dplyr::select(.data$group_biorep, .data$mean_ref_cq)

  return(reference_means)
}

#' Calculate target gene expression levels
#' @keywords internal
calculate_target_expression <- function(merged_data, reference_data, reference_gene) {

  expression_data <- merged_data %>%
    dplyr::filter(.data$gene != reference_gene) %>%
    dplyr::mutate(group_biorep = paste0(.data$group, .data$bio_rep)) %>%
    dplyr::left_join(reference_data, by = "group_biorep") %>%
    dplyr::select(-group_biorep) %>%
    dplyr::mutate(
      delta_ct = .data$cq - .data$mean_ref_cq,
      expre = 2^(.data$mean_ref_cq - .data$cq)  # Keep original variable name
    ) %>%
    dplyr::group_by(.data$group, .data$gene) %>%
    dplyr::mutate(
      n = dplyr::n(),
      mean.expre = mean(.data$expre, na.rm = TRUE),
      sd.expre = stats::sd(.data$expre, na.rm = TRUE),
      se.expre = .data$sd.expre / sqrt(.data$n)
    ) %>%
    dplyr::ungroup()

  return(expression_data)
}

#' Calculate summary statistics for expression data
#' @keywords internal
calculate_expression_summary <- function(expression_data) {

  summary_stats <- expression_data %>%
    dplyr::group_by(.data$group, .data$gene) %>%
    dplyr::summarise(
      n_replicates = dplyr::n(),
      mean_expression = mean(.data$expre, na.rm = TRUE),
      sd_expression = stats::sd(.data$expre, na.rm = TRUE),
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
      y = "Relative Expression (2^-dCt)",
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

# 为向后兼容保留原函数名和行为
CalExp2dCt <- function(cq.table, design.table, ref.gene = "Actin") {
  # 直接使用原有逻辑，保持完全兼容
  merged_data <- cq.table %>%
    dplyr::left_join(design.table, by = "Position") %>%
    dplyr::rename(
      position = Position,
      cq = Cq,
      group = Group,
      gene = Gene,
      biorep = BioRep
    )

  # 计算参考基因平均值
  reference_stats <- merged_data %>%
    dplyr::filter(.data$gene == ref.gene) %>%
    dplyr::group_by(.data$group, .data$biorep) %>%
    dplyr::mutate(
      mean.cq = mean(.data$cq),
      temp = paste0(.data$group, .data$biorep)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(.data$temp, .data$mean.cq) %>%
    dplyr::distinct_all()

  # 计算目标基因表达量
  expression_results <- merged_data %>%
    dplyr::filter(.data$gene != ref.gene) %>%
    dplyr::mutate(temp = paste0(.data$group, .data$biorep)) %>%
    dplyr::left_join(reference_stats, by = "temp") %>%
    dplyr::select(-temp) %>%
    dplyr::mutate(expre = 2^(.data$mean.cq - .data$cq)) %>%
    dplyr::group_by(.data$group, .data$gene) %>%
    dplyr::mutate(
      n = dplyr::n(),
      mean.expre = mean(.data$expre),
      sd.expre = stats::sd(.data$expre),
      se.expre = .data$sd.expre / sqrt(.data$n)
    ) %>%
    dplyr::ungroup()

  return(expression_results)
}

# Output column explanations:
#
# For calc_expression_delta_ct():
# expression_data:
# - position: Original well position from the qPCR plate
# - cq: Quantification cycle value for the target gene
# - group: Experimental treatment or condition group
# - gene: Name of the target gene being analyzed
# - bio_rep: Biological replicate identifier
# - mean_ref_cq: Mean Cq value of the reference gene for the same biological replicate
# - delta_ct: Delta Ct value (target Cq - reference Cq)
# - expre: Relative expression calculated as 2^(mean_ref_cq - cq) using 2^-ΔCt method
# - n: Number of technical replicates for each group-gene combination
# - mean.expre: Mean relative expression across technical replicates
# - sd.expre: Standard deviation of relative expression values
# - se.expre: Standard error of relative expression (sd.expre / sqrt(n))
#
# summary_table:
# - group: Experimental treatment or condition group
# - gene: Name of the target gene being analyzed
# - n_replicates: Number of technical replicates
# - mean_expression: Mean relative expression level
# - sd_expression: Standard deviation of expression values
# - se_expression: Standard error of expression values
# - mean_delta_ct: Mean delta Ct value
# - sd_delta_ct: Standard deviation of delta Ct values
#
# For CalExp2dCt() (backward compatibility):
# Same as expression_data above, maintaining exact original output format
