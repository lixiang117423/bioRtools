#' Calculate Gene Expression Using Standard Curve Method
#'
#' @description
#' This function calculates gene expression levels using standard curves,
#' with optional normalization to a reference gene and statistical comparisons
#' between treatment groups.
#'
#' @param cq_table A data frame containing position, Cq values, and gene information.
#'   Must contain columns: Position, Cq, Gene
#' @param curve_table A data frame containing standard curve parameters for each gene.
#'   Must contain columns: Gene, Slope, Intercept, min.Cq, max.Cq
#' @param design_table A data frame containing experimental design information.
#'   Must contain columns: Position, Treatment (or Group)
#' @param normalize_by_reference Logical indicating whether to normalize expression
#'   values by reference gene (default: TRUE)
#' @param reference_gene Character string specifying the reference gene name
#'   (default: "OsUBQ")
#' @param statistical_method Statistical method for group comparisons.
#'   Options: "t.test", "wilcox.test", "anova" (default: "t.test")
#' @param reference_group Character string specifying the reference group for
#'   statistical comparisons (default: "CK")
#' @param plot_type Type of plot to generate. Options: "box", "bar" (default: "box")
#' @param plot_ncol Number of columns in faceted plot (default: NULL for auto)
#'
#' @return A list containing:
#'   \item{expression_data}{Data frame with calculated expression values and statistics}
#'   \item{summary_table}{Summary statistics by group and gene}
#'   \item{statistical_results}{Statistical test results}
#'   \item{curve_warnings}{Warnings about Cq values outside curve range}
#'   \item{plot}{ggplot object showing expression levels}
#'
#' @importFrom dplyr left_join filter group_by mutate ungroup select summarise
#'   case_when n rename
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_bar geom_errorbar geom_jitter
#'   geom_text geom_hline facet_wrap labs theme element_text
#' @importFrom ggthemes theme_pander
#' @importFrom rstatix t_test wilcox_test
#' @importFrom multcomp glht mcp cld
#' @importFrom stats aov sd
#' @importFrom magrittr set_colnames
#' @importFrom tibble rownames_to_column
#'
#' @export
#'
#' @examples
#' # Load example data
#' cq_data_path <- system.file("examples", "cal.exp.curve.cq.txt", package = "qPCRtools")
#' curve_data_path <- system.file("examples", "cal.expre.curve.sdc.txt", package = "qPCRtools")
#' design_data_path <- system.file("examples", "cal.exp.curve.design.txt", package = "qPCRtools")
#' 
#' cq_data <- read.table(cq_data_path, header = TRUE)
#' curve_data <- read.table(curve_data_path, sep = "\t", header = TRUE)
#' design_data <- read.table(design_data_path, header = TRUE)
#' 
#' # Calculate expression using standard curves
#' result <- calc_expression_standard_curve(
#'   cq_table = cq_data,
#'   curve_table = curve_data,
#'   design_table = design_data,
#'   normalize_by_reference = TRUE,
#'   reference_gene = "OsUBQ",
#'   statistical_method = "t.test",
#'   reference_group = "CK",
#'   plot_type = "box"
#' )
#' 
#' # View results
#' result$summary_table
#' result$statistical_results
#' result$plot
#'
#' @author Xiang LI <lixiang117423@gmail.com>

calc_expression_standard_curve <- function(cq_table,
                                         curve_table,
                                         design_table,
                                         normalize_by_reference = TRUE,
                                         reference_gene = "OsUBQ",
                                         statistical_method = "t.test",
                                         reference_group = "CK",
                                         plot_type = "box",
                                         plot_ncol = NULL) {
  
  # Input validation
  validate_curve_inputs(cq_table, curve_table, design_table, reference_gene, 
                       reference_group, statistical_method, plot_type)
  
  # Merge data and calculate expression
  expression_data <- calculate_curve_expression(cq_table, curve_table, design_table)
  
  # Check for Cq values outside curve range
  curve_warnings <- check_curve_range(expression_data)
  
  # Apply reference gene normalization if requested
  if (normalize_by_reference) {
    expression_data <- normalize_by_reference_gene(expression_data, reference_gene)
  } else {
    expression_data <- expression_data %>%
      dplyr::select(.data$treatment, .data$gene, .data$expression) %>%
      dplyr::filter(.data$gene != reference_gene)
  }
  
  # Calculate summary statistics
  summary_table <- calculate_curve_summary(expression_data)
  
  # Perform statistical analysis
  statistical_results <- perform_curve_statistical_analysis(expression_data, reference_group, statistical_method)
  
  # Add statistical annotations
  summary_with_stats <- add_curve_statistical_annotations(summary_table, statistical_results)
  
  # Create plot
  plot_result <- create_curve_expression_plot(summary_with_stats, plot_type, plot_ncol)
  
  return(list(
    expression_data = expression_data,
    summary_table = summary_with_stats,
    statistical_results = statistical_results,
    curve_warnings = curve_warnings,
    plot = plot_result
  ))
}

#' Validate input parameters for standard curve calculation
#' @keywords internal
validate_curve_inputs <- function(cq_table, curve_table, design_table, reference_gene, 
                                 reference_group, statistical_method, plot_type) {
  
  if (!is.data.frame(cq_table)) {
    stop("cq_table must be a data frame")
  }
  
  if (!is.data.frame(curve_table)) {
    stop("curve_table must be a data frame")
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
  
  # Check required columns in curve_table
  required_curve_cols <- c("Gene", "Slope", "Intercept", "min.Cq", "max.Cq")
  missing_curve_cols <- setdiff(required_curve_cols, names(curve_table))
  if (length(missing_curve_cols) > 0) {
    stop(sprintf("curve_table missing required columns: %s", 
                paste(missing_curve_cols, collapse = ", ")))
  }
  
  # Check required columns in design_table (flexible naming)
  if ("Treatment" %in% names(design_table)) {
    design_treatment_col <- "Treatment"
  } else if ("Group" %in% names(design_table)) {
    design_treatment_col <- "Group"
  } else {
    stop("design_table must contain either 'Treatment' or 'Group' column")
  }
  
  if (!"Position" %in% names(design_table)) {
    stop("design_table must contain 'Position' column")
  }
  
  # Check if reference gene exists in both tables
  if (!reference_gene %in% cq_table$Gene) {
    stop(sprintf("Reference gene '%s' not found in cq_table", reference_gene))
  }
  
  if (!reference_gene %in% curve_table$Gene) {
    stop(sprintf("Reference gene '%s' not found in curve_table", reference_gene))
  }
  
  # Check if reference group exists
  if (!reference_group %in% design_table[[design_treatment_col]]) {
    stop(sprintf("Reference group '%s' not found in design_table", reference_group))
  }
  
  # Validate statistical method
  valid_methods <- c("t.test", "wilcox.test", "anova")
  if (!statistical_method %in% valid_methods) {
    stop(sprintf("statistical_method must be one of: %s", 
                paste(valid_methods, collapse = ", ")))
  }
  
  # Validate plot type
  valid_plot_types <- c("box", "bar")
  if (!plot_type %in% valid_plot_types) {
    stop(sprintf("plot_type must be one of: %s", 
                paste(valid_plot_types, collapse = ", ")))
  }
}

#' Calculate expression using standard curves
#' @keywords internal
calculate_curve_expression <- function(cq_table, curve_table, design_table) {
  
  # Determine treatment column name
  treatment_col <- if ("Treatment" %in% names(design_table)) "Treatment" else "Group"
  
  # Merge data and calculate expression
  merged_data <- cq_table %>%
    dplyr::left_join(design_table, by = "Position") %>%
    dplyr::left_join(curve_table, by = "Gene") %>%
    dplyr::filter(!is.na(.data$Cq), !is.na(.data$Slope), !is.na(.data$Intercept))
  
  # Standardize column names
  if (treatment_col == "Group") {
    merged_data <- merged_data %>%
      dplyr::rename(Treatment = .data$Group)
  }
  
  # Calculate expression values using standard curve equation
  expression_data <- merged_data %>%
    dplyr::mutate(
      out_of_range = dplyr::case_when(
        .data$Cq > .data$max.Cq | .data$Cq < .data$min.Cq ~ "yes",
        TRUE ~ "no"
      ),
      # Expression = 10^((Cq - Intercept) / Slope) for standard curve
      expression = 10^((.data$Cq - .data$Intercept) / .data$Slope)
    ) %>%
    dplyr::select(.data$Position, .data$Treatment, .data$Gene, .data$Cq, 
                 .data$expression, .data$out_of_range) %>%
    dplyr::rename(treatment = .data$Treatment, gene = .data$Gene)
  
  return(expression_data)
}

#' Check for Cq values outside curve range
#' @keywords internal
check_curve_range <- function(expression_data) {
  
  out_of_range_data <- expression_data %>%
    dplyr::filter(.data$out_of_range == "yes")
  
  warnings_list <- list()
  
  if (nrow(out_of_range_data) > 0) {
    warnings_list$out_of_range_positions <- out_of_range_data$Position
    warnings_list$message <- sprintf(
      "Warning: Cq values for positions %s are outside standard curve range", 
      paste(out_of_range_data$Position, collapse = ", ")
    )
    
    warning(warnings_list$message)
  }
  
  return(warnings_list)
}

#' Normalize expression by reference gene
#' @keywords internal
normalize_by_reference_gene <- function(expression_data, reference_gene) {
  
  # Calculate reference gene means by treatment
  reference_means <- expression_data %>%
    dplyr::filter(.data$gene == reference_gene) %>%
    dplyr::group_by(.data$treatment) %>%
    dplyr::summarise(mean_reference_expression = mean(.data$expression, na.rm = TRUE), 
                    .groups = "drop")
  
  # Normalize target genes by reference gene
  normalized_data <- expression_data %>%
    dplyr::filter(.data$gene != reference_gene) %>%
    dplyr::left_join(reference_means, by = "treatment") %>%
    dplyr::mutate(
      normalized_expression = .data$expression / .data$mean_reference_expression
    ) %>%
    dplyr::select(.data$treatment, .data$gene, expression = .data$normalized_expression)
  
  return(normalized_data)
}

#' Calculate summary statistics
#' @keywords internal
calculate_curve_summary <- function(expression_data) {
  
  summary_stats <- expression_data %>%
    dplyr::group_by(.data$treatment, .data$gene) %>%
    dplyr::summarise(
      n_replicates = dplyr::n(),
      mean_expression = mean(.data$expression, na.rm = TRUE),
      sd_expression = stats::sd(.data$expression, na.rm = TRUE),
      se_expression = .data$sd_expression / sqrt(.data$n_replicates),
      min_expression = min(.data$expression, na.rm = TRUE),
      max_expression = max(.data$expression, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      sd_expression = ifelse(is.na(.data$sd_expression), 0, .data$sd_expression),
      se_expression = ifelse(is.na(.data$se_expression), 0, .data$se_expression)
    )
  
  return(summary_stats)
}

#' Perform statistical analysis for standard curve data
#' @keywords internal
perform_curve_statistical_analysis <- function(expression_data, reference_group, statistical_method) {
  
  if (statistical_method == "t.test") {
    stat_results <- expression_data %>%
      dplyr::group_by(.data$gene) %>%
      rstatix::t_test(expression ~ treatment, ref.group = reference_group) %>%
      dplyr::ungroup() %>%
      dplyr::select(.data$gene, .data$group2, .data$p) %>%
      dplyr::rename(treatment = .data$group2)
    
  } else if (statistical_method == "wilcox.test") {
    stat_results <- expression_data %>%
      dplyr::group_by(.data$gene) %>%
      rstatix::wilcox_test(expression ~ treatment, ref.group = reference_group) %>%
      dplyr::ungroup() %>%
      dplyr::select(.data$gene, .data$group2, .data$p) %>%
      dplyr::rename(treatment = .data$group2)
    
  } else { # anova with Tukey post-hoc
    genes <- unique(expression_data$gene)
    results_list <- vector("list", length(genes))
    
    for (i in seq_along(genes)) {
      gene_data <- expression_data %>%
        dplyr::filter(.data$gene == genes[i]) %>%
        dplyr::mutate(treatment = factor(.data$treatment))
      
      if (length(unique(gene_data$treatment)) > 1) {
        fit <- stats::aov(expression ~ treatment, data = gene_data)
        tukey_test <- multcomp::glht(fit, linfct = multcomp::mcp(treatment = "Tukey"))
        letters <- multcomp::cld(tukey_test, level = 0.95, decreasing = TRUE)
        
        results_list[[i]] <- data.frame(
          gene = genes[i],
          treatment = names(letters$mcletters$Letters),
          significance_letter = as.character(letters$mcletters$Letters),
          stringsAsFactors = FALSE
        )
      }
    }
    
    stat_results <- do.call(rbind, results_list)
  }
  
  return(stat_results)
}

#' Add statistical annotations to summary table
#' @keywords internal
add_curve_statistical_annotations <- function(summary_table, statistical_results) {
  
  if ("significance_letter" %in% names(statistical_results)) {
    # For ANOVA results
    annotated_summary <- summary_table %>%
      dplyr::left_join(statistical_results, by = c("treatment", "gene")) %>%
      dplyr::rename(significance = .data$significance_letter)
    
  } else {
    # For t-test and wilcox test results
    stat_with_significance <- statistical_results %>%
      dplyr::mutate(
        significance = dplyr::case_when(
          .data$p < 0.001 ~ "***",
          .data$p < 0.01 ~ "**",
          .data$p < 0.05 ~ "*",
          TRUE ~ "NS"
        )
      )
    
    # Add reference group with no significance
    ref_groups <- summary_table %>%
      dplyr::anti_join(stat_with_significance, by = c("treatment", "gene")) %>%
      dplyr::select(.data$treatment, .data$gene) %>%
      dplyr::mutate(significance = "ref")
    
    all_significance <- rbind(
      stat_with_significance %>% dplyr::select(.data$treatment, .data$gene, .data$significance),
      ref_groups
    )
    
    annotated_summary <- summary_table %>%
      dplyr::left_join(all_significance, by = c("treatment", "gene"))
  }
  
  return(annotated_summary)
}

#' Create standard curve expression plot
#' @keywords internal
create_curve_expression_plot <- function(summary_data, plot_type, plot_ncol) {
  
  if (plot_type == "box") {
    p <- ggplot2::ggplot(summary_data, ggplot2::aes(.data$treatment, .data$mean_expression, 
                                                    fill = .data$treatment)) +
      ggplot2::geom_col(alpha = 0.7, width = 0.6) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data$mean_expression - .data$se_expression,
                    ymax = .data$mean_expression + .data$se_expression),
        width = 0.2
      ) +
      ggplot2::geom_text(
        ggplot2::aes(label = .data$significance),
        vjust = -0.5,
        size = 4
      ) +
      ggplot2::facet_wrap(~ .data$gene, scales = "free_y", ncol = plot_ncol) +
      ggplot2::labs(
        x = "Treatment",
        y = "Relative Expression",
        title = "Gene Expression Analysis (Standard Curve Method)"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.position = "none",
        strip.text.x = ggplot2::element_text(face = "italic"),
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )
    
  } else { # bar plot
    p <- ggplot2::ggplot(summary_data, ggplot2::aes(.data$treatment, .data$mean_expression, 
                                                    fill = .data$treatment)) +
      ggplot2::geom_bar(stat = "identity", width = 0.6, alpha = 0.7) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data$mean_expression - .data$sd_expression,
                    ymax = .data$mean_expression + .data$sd_expression),
        width = 0.2
      ) +
      ggplot2::geom_text(
        ggplot2::aes(y = .data$max_expression * 1.08, 
                    label = .data$significance),
        size = 4,
        color = "red"
      ) +
      ggplot2::facet_wrap(~ .data$gene, scales = "free_y", ncol = plot_ncol) +
      ggplot2::labs(
        x = "Treatment",
        y = "Relative Expression",
        title = "Gene Expression Analysis (Standard Curve Method)"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.position = "none",
        strip.text.x = ggplot2::element_text(face = "italic"),
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )
  }
  
  return(p)
}

# For backward compatibility
CalExpCurve <- function(cq.table, curve.table, design.table, 
                       correction = TRUE, ref.gene = "OsUBQ", 
                       stat.method = "t.test", ref.group = "CK", 
                       fig.type = "box", fig.ncol = NULL) {
  
  result <- calc_expression_standard_curve(
    cq_table = cq.table,
    curve_table = curve.table,
    design_table = design.table,
    normalize_by_reference = correction,
    reference_gene = ref.gene,
    statistical_method = stat.method,
    reference_group = ref.group,
    plot_type = fig.type,
    plot_ncol = fig.ncol
  )
  
  # Return in original format for compatibility
  return(list(
    table = result$summary_table,
    figure = result$plot
  ))
}