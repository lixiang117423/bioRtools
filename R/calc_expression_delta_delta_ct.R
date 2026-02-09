#' Calculate Gene Expression Using Delta-Delta Ct Method
#'
#' @description
#' This function calculates relative gene expression using the delta-delta Ct
#' method (2^-ddCt), which normalizes target gene expression to a reference
#' gene and compares to a reference group.
#'
#' @param cq_table A data frame containing position and Cq values.
#'   Must contain columns: Position, Cq
#' @param design_table A data frame containing experimental design information.
#'   Must contain columns: Position, Group, Gene, BioRep
#' @param reference_gene Character string specifying the reference gene name
#'   (default: "OsUBQ")
#' @param reference_group Character string specifying the reference group name
#'   (default: "CK")
#' @param statistical_method Statistical method for group comparisons.
#'   Options: "t.test", "wilcox.test", "anova" (default: "t.test")
#' @param remove_outliers Logical indicating whether to remove outliers using
#'   IQR method (default: TRUE)
#' @param plot_type Type of plot to generate. Options: "box", "bar" (default: "box")
#' @param plot_ncol Number of columns in faceted plot (default: NULL for auto)
#'
#' @return A list containing:
#'   \item{expression_data}{Data frame with calculated expression values and statistics}
#'   \item{summary_table}{Summary statistics by group and gene}
#'   \item{statistical_results}{Statistical test results}
#'   \item{plot}{ggplot object showing expression levels}
#'
#' @importFrom dplyr left_join filter group_by mutate ungroup select rename
#'   case_when add_row n
#' @importFrom tidyr pivot_wider
#' @importFrom magrittr set_names
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_bar geom_errorbar geom_jitter
#'   geom_text geom_hline facet_wrap labs theme element_text
#' @importFrom ggthemes theme_pander
#' @importFrom rstatix t_test wilcox_test
#' @importFrom multcomp glht mcp cld
#' @importFrom stats aov sd
#'
#' @export
#'
#' @examples
#' library(bioRtools)
#' # Load example data
#' cq_data_path <- system.file("extdata/qPCR", "ddct.cq.txt", package = "bioRtools")
#' design_data_path <- system.file("extdata/qPCR", "ddct.design.txt", package = "bioRtools")
#' cq_data <- read.table(cq_data_path, header = TRUE)
#' design_data <- read.table(design_data_path, header = TRUE)
#'
#' # Calculate expression using delta-delta Ct method
#' result <- calc_expression_delta_delta_ct(
#'   cq_table = cq_data,
#'   design_table = design_data,
#'   reference_gene = "OsUBQ",
#'   reference_group = "CK",
#'   statistical_method = "t.test",
#'   remove_outliers = TRUE,
#'   plot_type = "bar"
#' )
#'
#' # View results
#' result$summary_table
#' result$statistical_results
#' result$plot
#'
#' @author Xiang LI <lixiang117423@gmail.com>

calc_expression_delta_delta_ct <- function(cq_table,
                                           design_table,
                                           reference_gene = "OsUBQ",
                                           reference_group = "CK",
                                           statistical_method = "t.test",
                                           remove_outliers = TRUE,
                                           plot_type = "box",
                                           plot_ncol = NULL) {
  # Prepare data first (merge tables before validation)
  merged_data <- prepare_ddct_data(cq_table, design_table)

  # Basic validation (after merging)
  validate_ddct_inputs_fixed(merged_data, reference_gene, reference_group,
    statistical_method, plot_type)

  # Calculate delta-delta Ct values
  expression_data <- calculate_ddct_expression(merged_data, reference_gene, reference_group)

  # Remove outliers if requested
  if (remove_outliers) {
    expression_data <- remove_expression_outliers(expression_data)
  }

  # Calculate summary statistics
  summary_table <- calculate_ddct_summary(expression_data)

  # Perform statistical analysis
  statistical_results <- perform_statistical_analysis(expression_data, reference_group, statistical_method)

  # Add statistical annotations to summary table
  summary_with_stats <- add_statistical_annotations(summary_table, statistical_results)

  # Create plot
  plot_result <- create_ddct_plot(summary_with_stats, plot_type, plot_ncol)

  return(list(
    expression_data = expression_data,
    summary_table = summary_with_stats,
    statistical_results = statistical_results,
    plot = plot_result
  ))
}

#' Prepare and merge data for delta-delta Ct calculation
#' @keywords internal
prepare_ddct_data <- function(cq_table, design_table) {

  merged_data <- cq_table %>%
    dplyr::left_join(design_table, by = "Position") %>%
    dplyr::rename(
      position = Position,
      cq = Cq,
      group = Group,
      gene = Gene,
      bio_rep = BioRep
    ) %>%
    dplyr::filter(!is.na(.data$group), !is.na(.data$bio_rep), !is.na(.data$cq))

  return(merged_data)
}

#' Validate input parameters for delta-delta Ct calculation (after merging)
#' @keywords internal
validate_ddct_inputs_fixed <- function(merged_data, reference_gene,
                                       reference_group, statistical_method, plot_type) {
  # Check if reference gene exists in merged data
  if (!reference_gene %in% merged_data$gene) {
    available_genes <- unique(merged_data$gene)
    stop(sprintf("Reference gene '%s' not found. Available genes: %s",
      reference_gene, paste(available_genes, collapse = ", ")))
  }

  # Check if reference group exists in merged data
  if (!reference_group %in% merged_data$group) {
    available_groups <- unique(merged_data$group)
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

#' Calculate delta-delta Ct expression values
#' @keywords internal
calculate_ddct_expression <- function(merged_data, reference_gene, reference_group) {
  # Get target genes (all except reference gene)
  target_genes <- setdiff(unique(merged_data$gene), reference_gene)

  # Calculate reference values (reference gene in reference group)
  ref_gene_ref_group_cq <- merged_data %>%
    dplyr::filter(.data$gene == reference_gene, .data$group == reference_group) %>%
    dplyr::pull(.data$cq) %>%
    mean(na.rm = TRUE)

  expression_results <- vector("list", length(target_genes))

  for (i in seq_along(target_genes)) {
    target_gene <- target_genes[i]

    # Calculate dCt for reference group (target - reference in reference group)
    target_gene_ref_group_cq <- merged_data %>%
      dplyr::filter(.data$gene == target_gene, .data$group == reference_group) %>%
      dplyr::pull(.data$cq) %>%
      mean(na.rm = TRUE)

    delta_ct_ref <- target_gene_ref_group_cq - ref_gene_ref_group_cq

    # Calculate expression for all groups
    gene_data <- merged_data %>%
      dplyr::filter(.data$gene %in% c(target_gene, reference_gene)) %>%
      dplyr::group_by(.data$group, .data$bio_rep, .data$gene) %>%
      dplyr::summarise(mean_cq = mean(.data$cq, na.rm = TRUE), .groups = "drop") %>%
      tidyr::pivot_wider(
        id_cols = c("group", "bio_rep"),
        names_from = "gene",
        values_from = "mean_cq"
      )

    # Ensure column order and names
    if (reference_gene %in% names(gene_data) && target_gene %in% names(gene_data)) {
      gene_data <- gene_data %>%
        dplyr::rename(
          target_cq = !!rlang::sym(target_gene),
          reference_cq = !!rlang::sym(reference_gene)
        ) %>%
        dplyr::mutate(
          delta_ct = .data$target_cq - .data$reference_cq,
          delta_delta_ct = .data$delta_ct - delta_ct_ref,
          relative_expression = 2^(-.data$delta_delta_ct),
          gene = target_gene
        ) %>%
        dplyr::select(.data$group, .data$gene, .data$bio_rep, .data$relative_expression,
          .data$delta_ct, .data$delta_delta_ct)

      expression_results[[i]] <- gene_data
    }
  }

  # Combine all results
  all_expression_data <- do.call(rbind, expression_results)

  return(all_expression_data)
}

#' Remove outliers using IQR method
#' @keywords internal
remove_expression_outliers <- function(expression_data) {

  find_outliers <- function(x) {
    q1 <- stats::quantile(x, 0.25, na.rm = TRUE)
    q3 <- stats::quantile(x, 0.75, na.rm = TRUE)
    iqr <- q3 - q1
    lower_bound <- q1 - 1.5 * iqr
    upper_bound <- q3 + 1.5 * iqr

    return(x < lower_bound | x > upper_bound)
  }

  cleaned_data <- expression_data %>%
    dplyr::group_by(.data$group, .data$gene) %>%
    dplyr::filter(!find_outliers(.data$relative_expression)) %>%
    dplyr::ungroup()

  return(cleaned_data)
}

#' Calculate summary statistics
#' @keywords internal
calculate_ddct_summary <- function(expression_data) {

  summary_stats <- expression_data %>%
    dplyr::group_by(.data$group, .data$gene) %>%
    dplyr::summarise(
      n_replicates = dplyr::n(),
      mean_expression = mean(.data$relative_expression, na.rm = TRUE),
      sd_expression = stats::sd(.data$relative_expression, na.rm = TRUE),
      se_expression = .data$sd_expression / sqrt(.data$n_replicates),
      mean_delta_ct = mean(.data$delta_ct, na.rm = TRUE),
      mean_delta_delta_ct = mean(.data$delta_delta_ct, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      sd_expression = ifelse(is.na(.data$sd_expression), 0, .data$sd_expression),
      se_expression = ifelse(is.na(.data$se_expression), 0, .data$se_expression)
    )

  return(summary_stats)
}

#' Perform statistical analysis
#' @keywords internal
perform_statistical_analysis <- function(expression_data, reference_group, statistical_method) {

  tryCatch(
    {
      if (statistical_method == "t.test") {
        stat_results <- expression_data %>%
          dplyr::group_by(.data$gene) %>%
          rstatix::t_test(relative_expression ~ group, ref.group = reference_group) %>%
          dplyr::ungroup() %>%
          dplyr::select(.data$gene, .data$group2, .data$p) %>%
          dplyr::rename(group = .data$group2)

      } else if (statistical_method == "wilcox.test") {
        stat_results <- expression_data %>%
          dplyr::group_by(.data$gene) %>%
          rstatix::wilcox_test(relative_expression ~ group, ref.group = reference_group) %>%
          dplyr::ungroup() %>%
          dplyr::select(.data$gene, .data$group2, .data$p) %>%
          dplyr::rename(group = .data$group2)

      } else { # anova with Tukey post-hoc
        genes <- unique(expression_data$gene)
        results_list <- vector("list", length(genes))

        for (i in seq_along(genes)) {
          gene_data <- expression_data %>%
            dplyr::filter(.data$gene == genes[i]) %>%
            dplyr::mutate(group = factor(.data$group))

          if (length(unique(gene_data$group)) > 1) {
            fit <- stats::aov(relative_expression ~ group, data = gene_data)
            tukey_test <- multcomp::glht(fit, linfct = multcomp::mcp(group = "Tukey"))
            letters <- multcomp::cld(tukey_test, level = 0.95, decreasing = TRUE)

            results_list[[i]] <- data.frame(
              gene = genes[i],
              group = names(letters$mcletters$Letters),
              significance_letter = as.character(letters$mcletters$Letters),
              stringsAsFactors = FALSE
            )
          }
        }

        stat_results <- do.call(rbind, results_list)
      }

      return(stat_results)

    },
    error = function(e) {
      warning("Statistical analysis failed: ", e$message)
      return(data.frame(gene = character(0), group = character(0), p = numeric(0)))
    })
}

#' Add statistical annotations to summary table
#' @keywords internal
add_statistical_annotations <- function(summary_table, statistical_results) {

  if (nrow(statistical_results) == 0) {
    # If statistical analysis failed, add NS for all
    summary_table$significance <- "NS"
    return(summary_table)
  }

  if ("significance_letter" %in% names(statistical_results)) {
    # For ANOVA results
    annotated_summary <- summary_table %>%
      dplyr::left_join(statistical_results, by = c("group", "gene")) %>%
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
      dplyr::anti_join(stat_with_significance, by = c("group", "gene")) %>%
      dplyr::select(.data$group, .data$gene) %>%
      dplyr::mutate(significance = "")

    all_significance <- rbind(
      stat_with_significance %>% dplyr::select(.data$group, .data$gene, .data$significance),
      ref_groups
    )

    annotated_summary <- summary_table %>%
      dplyr::left_join(all_significance, by = c("group", "gene"))
  }

  # Ensure significance column exists and has default values
  if (!"significance" %in% names(annotated_summary)) {
    annotated_summary$significance <- "NS"
  }

  annotated_summary$significance[is.na(annotated_summary$significance)] <- "NS"

  return(annotated_summary)
}

#' Create delta-delta Ct plot
#' @keywords internal
create_ddct_plot <- function(summary_data, plot_type, plot_ncol) {

  if (nrow(summary_data) == 0) {
    warning("No data available for plotting")
    return(NULL)
  }

  # Prepare data for plotting
  plot_data <- summary_data %>%
    dplyr::rename(
      treatment = .data$group,
      expression = .data$mean_expression,
      sd_expr = .data$sd_expression,
      se_expr = .data$se_expression
    )

  if (plot_type == "bar") {
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(.data$treatment, .data$expression,
      fill = .data$treatment)) +
      ggplot2::geom_col(alpha = 0.7, width = 0.6) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = pmax(0, .data$expression - .data$se_expr),
          ymax = .data$expression + .data$se_expr),
        width = 0.2
      ) +
      ggplot2::geom_text(
        ggplot2::aes(y = .data$expression + .data$se_expr + max(.data$expression) * 0.05,
          label = .data$significance),
        vjust = 0,
        size = 4
      ) +
      ggplot2::facet_wrap(~ .data$gene, scales = "free_y", ncol = plot_ncol) +
      ggplot2::labs(
        x = "Treatment",
        y = "Relative Expression (2^-ddCt)",
        title = "Gene Expression Analysis (Delta-Delta Ct Method)"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.position = "none",
        strip.text.x = ggplot2::element_text(face = "italic"),
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )

  } else { # box plot style with bars
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(.data$treatment, .data$expression,
      fill = .data$treatment)) +
      ggplot2::geom_col(alpha = 0.7, width = 0.6) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = pmax(0, .data$expression - .data$se_expr),
          ymax = .data$expression + .data$se_expr),
        width = 0.2
      ) +
      ggplot2::geom_text(
        ggplot2::aes(y = .data$expression + .data$se_expr + max(.data$expression) * 0.05,
          label = .data$significance),
        vjust = 0,
        size = 4
      ) +
      ggplot2::facet_wrap(~ .data$gene, scales = "free_y", ncol = plot_ncol) +
      ggplot2::labs(
        x = "Treatment",
        y = "Relative Expression (2^-ddCt)",
        title = "Gene Expression Analysis (Delta-Delta Ct Method)"
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

# For backward compatibility - simplified version that works with original data structure
CalExp2ddCt <- function(cq.table, design.table, ref.gene = "OsUBQ",
                        ref.group = "CK", stat.method = "t.test",
                        remove.outliers = TRUE, fig.type = "box",
                        fig.ncol = NULL) {

  tryCatch(
    {
      result <- calc_expression_delta_delta_ct(
        cq_table = cq.table,
        design_table = design.table,
        reference_gene = ref.gene,
        reference_group = ref.group,
        statistical_method = stat.method,
        remove_outliers = remove.outliers,
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
      warning("Analysis failed: ", e$message)
      return(list(
        table = data.frame(),
        figure = NULL
      ))
    })
}

# Output column explanations:
#
# expression_data:
# - group: Experimental treatment or condition group
# - gene: Name of the target gene being analyzed (excludes reference gene)
# - bio_rep: Biological replicate identifier
# - relative_expression: Relative expression calculated using 2^-ddCt method
# - delta_ct: Delta Ct value (target gene Cq - reference gene Cq)
# - delta_delta_ct: Delta-delta Ct value (delta_ct - delta_ct_reference_group)
#
# summary_table:
# - group: Experimental treatment or condition group
# - gene: Name of the target gene being analyzed
# - n_replicates: Number of biological replicates
# - mean_expression: Mean relative expression level (2^-ddCt)
# - sd_expression: Standard deviation of relative expression values
# - se_expression: Standard error of relative expression values
# - mean_delta_ct: Mean delta Ct value
# - mean_delta_delta_ct: Mean delta-delta Ct value
# - significance: Statistical significance annotation (*, **, ***, NS)
#
# statistical_results:
# - gene: Target gene name
# - group: Treatment group (compared to reference group)
# - p: P-value from statistical test
# - significance: Significance level or letters for multiple comparisons
