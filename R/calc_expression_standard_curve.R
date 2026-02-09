#' Calculate Gene Expression Using Standard Curve Method
#'
#' @description
#' This function calculates gene expression levels using standard curves,
#' with optional normalization to a reference gene and statistical comparisons
#' between treatment groups.
#'
#' @param cq_table A data frame containing position and Cq values.
#'   Must contain columns: Position, Cq
#' @param curve_table A data frame containing standard curve parameters for each gene.
#'   Must contain columns: Gene, Slope, Intercept, min.Cq, max.Cq
#' @param design_table A data frame containing experimental design information.
#'   Must contain columns: Position, Treatment (or Group), Gene
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
#' cq_data_path <- system.file("extdata/qPCR", "cal.exp.curve.cq.txt", package = "bioRtools")
#' curve_data_path <- system.file("extdata/qPCR", "cal.expre.curve.sdc.txt", package = "bioRtools")
#' design_data_path <- system.file("extdata/qPCR", "cal.exp.curve.design.txt", package = "bioRtools")
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
  # Merge data and calculate expression first
  expression_data <- calculate_curve_expression(cq_table, curve_table, design_table)

  # Validate merged data (after merging, not before)
  validate_curve_inputs_fixed(expression_data, curve_table, reference_gene,
    reference_group, statistical_method, plot_type)

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
  summary_with_stats <- add_curve_statistical_annotations(summary_table, statistical_results, reference_group)

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

#' Validate input parameters for standard curve calculation (after merging)
#' @keywords internal
validate_curve_inputs_fixed <- function(merged_data, curve_table, reference_gene,
                                        reference_group, statistical_method, plot_type) {

  if (nrow(merged_data) == 0) {
    stop("No data available after merging tables. Check that Position and Gene columns match between tables.")
  }

  # Check required columns in merged data
  required_merged_cols <- c("Position", "treatment", "gene", "Cq", "expression")
  missing_merged_cols <- setdiff(required_merged_cols, names(merged_data))
  if (length(missing_merged_cols) > 0) {
    stop(sprintf("Merged data missing required columns: %s",
      paste(missing_merged_cols, collapse = ", ")))
  }

  # Check if reference gene exists in merged data
  available_genes <- unique(merged_data$gene)
  if (!reference_gene %in% available_genes) {
    stop(sprintf("Reference gene '%s' not found. Available genes: %s",
      reference_gene, paste(available_genes, collapse = ", ")))
  }

  # Check if reference gene has curve parameters
  if (!reference_gene %in% curve_table$Gene) {
    stop(sprintf("Reference gene '%s' not found in curve_table", reference_gene))
  }

  # Check if reference group exists
  available_groups <- unique(merged_data$treatment)
  if (!reference_group %in% available_groups) {
    stop(sprintf("Reference group '%s' not found. Available groups: %s",
      reference_group, paste(available_groups, collapse = ", ")))
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

  tryCatch(
    {
      if (statistical_method == "t.test") {
        stat_results <- expression_data %>%
          dplyr::group_by(.data$gene) %>%
          rstatix::t_test(expression ~ treatment, ref.group = reference_group) %>%
          dplyr::ungroup() %>%
          dplyr::select(.data$gene, .data$group2, .data$p, .data$statistic, .data$df) %>%
          dplyr::mutate(
            p_value = round(.data$p, 6),
            t_statistic = round(.data$statistic, 4)
          ) %>%
          dplyr::rename(treatment = .data$group2) %>%
          dplyr::select(.data$gene, .data$treatment, .data$p_value, .data$t_statistic, .data$df)

      } else if (statistical_method == "wilcox.test") {
        stat_results <- expression_data %>%
          dplyr::group_by(.data$gene) %>%
          rstatix::wilcox_test(expression ~ treatment, ref.group = reference_group) %>%
          dplyr::ungroup() %>%
          dplyr::select(.data$gene, .data$group2, .data$p, .data$statistic) %>%
          dplyr::mutate(
            p_value = round(.data$p, 6),
            w_statistic = round(.data$statistic, 4)
          ) %>%
          dplyr::rename(treatment = .data$group2) %>%
          dplyr::select(.data$gene, .data$treatment, .data$p_value, .data$w_statistic)

      } else { # anova with Tukey post-hoc
        genes <- unique(expression_data$gene)
        results_list <- vector("list", length(genes))
        anova_results <- vector("list", length(genes))

        for (i in seq_along(genes)) {
          gene_data <- expression_data %>%
            dplyr::filter(.data$gene == genes[i]) %>%
            dplyr::mutate(treatment = factor(.data$treatment))

          if (length(unique(gene_data$treatment)) > 1) {
            # ANOVA F-test
            fit <- stats::aov(expression ~ treatment, data = gene_data)
            anova_summary <- summary(fit)
            f_value <- anova_summary[[1]]["treatment", "F value"]
            p_value <- anova_summary[[1]]["treatment", "Pr(>F)"]
            df1 <- anova_summary[[1]]["treatment", "Df"]
            df2 <- anova_summary[[1]]["Residuals", "Df"]

            anova_results[[i]] <- data.frame(
              gene = genes[i],
              f_statistic = round(f_value, 4),
              p_value_anova = round(p_value, 6),
              df1 = df1,
              df2 = df2,
              stringsAsFactors = FALSE
            )

            # Post-hoc Tukey test
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

        stat_results <- list(
          letters = do.call(rbind, results_list),
          anova_stats = do.call(rbind, anova_results)
        )
      }

      return(stat_results)

    },
    error = function(e) {
      warning("Statistical analysis failed: ", e$message)
      return(data.frame(gene = character(0), treatment = character(0), p_value = numeric(0)))
    })
}

#' Add statistical annotations to summary table
#' @keywords internal
add_curve_statistical_annotations <- function(summary_table, statistical_results, reference_group) {

  if (is.list(statistical_results) && "letters" %in% names(statistical_results)) {
    # For ANOVA results
    letters_data <- statistical_results$letters %>%
      # Make reference group letter empty
      dplyr::mutate(significance_letter = ifelse(.data$treatment == reference_group, "", .data$significance_letter))

    anova_stats <- statistical_results$anova_stats

    annotated_summary <- summary_table %>%
      dplyr::left_join(letters_data, by = c("treatment", "gene")) %>%
      dplyr::left_join(anova_stats, by = "gene") %>%
      dplyr::rename(significance = .data$significance_letter)

  } else if (nrow(statistical_results) > 0) {
    # For t-test and wilcox test results
    stat_with_significance <- statistical_results %>%
      dplyr::mutate(
        significance = dplyr::case_when(
          .data$p_value < 0.001 ~ "***",
          .data$p_value < 0.01 ~ "**",
          .data$p_value < 0.05 ~ "*",
          TRUE ~ "NS"
        )
      )

    # Add empty significance for reference group
    ref_groups <- summary_table %>%
      dplyr::filter(.data$treatment == reference_group) %>%
      dplyr::select(.data$treatment, .data$gene) %>%
      dplyr::mutate(
        significance = "",
        p_value = NA_real_
      )

    # Add statistical columns if they exist
    if ("t_statistic" %in% names(stat_with_significance)) {
      ref_groups$t_statistic <- NA_real_
      ref_groups$df <- NA_real_
    }
    if ("w_statistic" %in% names(stat_with_significance)) {
      ref_groups$w_statistic <- NA_real_
    }

    all_significance <- rbind(
      stat_with_significance %>% dplyr::select(.data$treatment, .data$gene, .data$significance, dplyr::everything()),
      ref_groups
    )

    annotated_summary <- summary_table %>%
      dplyr::left_join(all_significance, by = c("treatment", "gene"))

  } else {
    # No statistical results
    annotated_summary <- summary_table %>%
      dplyr::mutate(
        significance = ifelse(.data$treatment == reference_group, "", "NS"),
        p_value = NA_real_
      )
  }

  return(annotated_summary)
}

#' Create standard curve expression plot
#' @keywords internal
create_curve_expression_plot <- function(summary_data, plot_type, plot_ncol) {

  if (nrow(summary_data) == 0) {
    warning("No data available for plotting")
    return(NULL)
  }

  # Ensure significance column exists
  if (!"significance" %in% names(summary_data)) {
    summary_data$significance <- "NS"
  }

  if (plot_type == "box") {
    p <- ggplot2::ggplot(summary_data, ggplot2::aes(.data$treatment, .data$mean_expression,
      fill = .data$treatment)) +
      ggplot2::geom_col(alpha = 0.7, width = 0.6) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = pmax(0, .data$mean_expression - .data$se_expression),
          ymax = .data$mean_expression + .data$se_expression),
        width = 0.2
      ) +
      ggplot2::geom_text(
        ggplot2::aes(y = .data$mean_expression + .data$se_expression + max(.data$mean_expression, na.rm = TRUE) * 0.05,
          label = .data$significance),
        vjust = 0,
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
        ggplot2::aes(ymin = pmax(0, .data$mean_expression - .data$sd_expression),
          ymax = .data$mean_expression + .data$sd_expression),
        width = 0.2
      ) +
      ggplot2::geom_text(
        ggplot2::aes(y = pmax(.data$mean_expression + .data$sd_expression, .data$max_expression) * 1.08,
          label = .data$significance),
        size = 4,
        color = "black"
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

  tryCatch(
    {
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

    },
    error = function(e) {
      warning("Standard curve analysis failed: ", e$message)
      return(list(
        table = data.frame(),
        figure = NULL
      ))
    })
}

# Output column explanations:
#
# expression_data:
# - treatment: Experimental treatment or condition group
# - gene: Name of the target gene being analyzed
# - expression: Calculated expression value (normalized if reference gene used)
#
# summary_table:
# - treatment: Experimental treatment or condition group
# - gene: Name of the target gene being analyzed
# - n_replicates: Number of replicates per group
# - mean_expression: Mean expression level
# - sd_expression: Standard deviation of expression values
# - se_expression: Standard error of expression values
# - min_expression: Minimum expression value
# - max_expression: Maximum expression value
# - significance: Statistical significance annotation (*, **, ***, NS, letters; empty for reference group)
#
# Additional statistical columns (depending on method):
# For t.test:
# - p_value: P-value from t-test (NA for reference group)
# - t_statistic: T-statistic value (NA for reference group)
# - df: Degrees of freedom (NA for reference group)
#
# For wilcox.test:
# - p_value: P-value from Wilcoxon test (NA for reference group)
# - w_statistic: W-statistic value (NA for reference group)
#
# For anova:
# - f_statistic: F-statistic from ANOVA
# - p_value_anova: P-value from ANOVA F-test
# - df1, df2: Degrees of freedom
#
# curve_warnings:
# - out_of_range_positions: Positions with Cq values outside curve range
# - message: Warning message about out-of-range values
#
# Expression calculation:
# 1. Standard curve: expression = 10^((Cq - Intercept) / Slope)
# 2. Reference normalization: normalized_expression = target_expression / reference_expression
# 3. Quality control: warnings for Cq values outside calibration range
