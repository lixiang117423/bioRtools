#' Calculate Gene Expression Using qPCR with Efficiency Correction
#'
#' @description
#' Calculate expression using qPCR with efficiency correction and automatic
#' reference gene selection using GeNorm algorithm if reference genes are not specified.
#' This function preserves the original algorithm logic from the author.
#'
#' @param cq_table A data frame containing qPCR data.
#'   Must contain columns: Position, Cq
#' @param design_table A data frame containing experimental design information.
#'   Must contain columns: Position, Group, Gene, BioRep, TechRep, Eff
#' @param reference_gene Character vector of reference gene names.
#'   If NULL, the most stable genes will be selected automatically
#' @param reference_group Character string specifying the reference group for
#'   statistical comparisons (default: "CK")
#' @param statistical_method Statistical method for group comparisons.
#'   Options: "t.test", "wilcox.test", "anova" (default: "t.test")
#' @param plot_type Type of plot to generate. Options: "box", "bar" (default: "box")
#' @param plot_ncol Number of columns in faceted plot (default: NULL for auto)
#'
#' @return A list containing:
#'   \item{table}{Data frame with calculated expression values and statistics}
#'   \item{figure}{ggplot object showing expression levels}
#'
#' @importFrom dplyr left_join filter group_by mutate ungroup select summarise
#'   case_when n rename
#' @importFrom tidyr pivot_wider
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_bar geom_errorbar
#'   geom_text facet_wrap labs theme element_text
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
#' library(bioRtools)
#' # Load example data
#' cq_data_path <- system.file("extdata/qPCR", "cal.expre.rqpcr.cq.txt", package = "bioRtools")
#' design_data_path <- system.file("extdata/qPCR", "cal.expre.rqpcr.design.txt", package = "bioRtools")
#'
#' cq_data <- read.table(cq_data_path, header = TRUE)
#' design_data <- read.table(design_data_path, header = TRUE)
#'
#' # Calculate expression
#' result <- calc_expression_qpcr_efficiency(
#'   cq_table = cq_data,
#'   design_table = design_data,
#'   reference_gene = NULL,  # Automatic selection
#'   reference_group = "CK",
#'   statistical_method = "t.test",
#'   plot_type = "box"
#' )
#'
#' # View results
#' result$table
#' result$figure
#'
#' @author Xiang LI <lixiang117423@gmail.com>

calc_expression_qpcr_efficiency <- function(cq_table,
                                            design_table,
                                            reference_gene = NULL,
                                            reference_group = "CK",
                                            statistical_method = "t.test",
                                            plot_type = "box",
                                            plot_ncol = NULL) {
  # Merge data first - preserving original algorithm logic
  df <- merge_qpcr_data(cq_table, design_table)

  # Validate merged data (after merging, not before)
  validate_qpcr_inputs_fixed(df, statistical_method, plot_type)

  # Calculate expression using original algorithm
  df_expression <- calculate_qpcr_expression(df)

  # Find reference gene using original GeNorm algorithm if not provided
  if (is.null(reference_gene)) {
    reference_gene <- find_reference_genes_original(df_expression)
  }

  # Validate reference gene exists in merged data
  available_genes <- unique(df$gene)
  if (!all(reference_gene %in% available_genes)) {
    missing_genes <- setdiff(reference_gene, available_genes)
    stop(sprintf("Reference gene(s) '%s' not found. Available genes: %s",
      paste(missing_genes, collapse = ", "),
      paste(available_genes, collapse = ", ")))
  }

  # Calculate normalization factors using original algorithm
  normalization_factors <- calculate_normalization_factors_original(df_expression, reference_gene)

  # Calculate corrected expression using original algorithm
  results_all <- calculate_corrected_expression_original(df_expression, reference_gene, normalization_factors)

  # Perform statistical analysis using original method
  results_with_stats <- perform_statistical_analysis_original(results_all, reference_group, statistical_method)

  # Create plot using original approach
  plot_result <- create_plot_original(results_with_stats, plot_type, plot_ncol)

  # Clean up results to include statistical information
  final_results <- results_with_stats %>%
    dplyr::select(-.data$temp, -.data$efficiency) %>%
    dplyr::select(.data$group, .data$gene, .data$bio_rep, .data$expression_value,
      .data$mean_expression, .data$sd_expression, .data$se_expression, .data$significance,
      dplyr::everything()) %>%
    dplyr::mutate(
      across(c(.data$expression_value, .data$mean_expression, .data$sd_expression, .data$se_expression), ~ round(.x, 4)),
      across(matches("p_value|statistic"), ~ round(.x, 6))
    )

  list(table = final_results, figure = plot_result)
}

#' Merge data exactly as in original
#' @keywords internal
merge_qpcr_data <- function(cq_table, design_table) {

  df <- cq_table %>%
    dplyr::left_join(design_table, by = "Position") %>%
    dplyr::rename(
      position = Position,
      cq = Cq,
      group = Group,
      gene = Gene,
      bio_rep = BioRep,
      tech_rep = TechRep,
      efficiency = Eff
    )

  df
}

#' Validate merged data (after merging, not separate tables)
#' @keywords internal
validate_qpcr_inputs_fixed <- function(merged_data, statistical_method, plot_type) {

  if (!is.data.frame(merged_data)) {
    stop("merged_data must be a data frame")
  }

  # Check required columns in merged data
  required_cols <- c("position", "cq", "group", "gene", "bio_rep", "tech_rep", "efficiency")
  missing_cols <- setdiff(required_cols, names(merged_data))
  if (length(missing_cols) > 0) {
    stop(sprintf("Merged data missing required columns: %s",
      paste(missing_cols, collapse = ", ")))
  }

  # Check for basic data issues
  if (nrow(merged_data) == 0) {
    stop("No data available after merging tables")
  }

  if (all(is.na(merged_data$cq))) {
    stop("All Cq values are missing")
  }

  if (all(is.na(merged_data$efficiency))) {
    stop("All efficiency values are missing")
  }

  valid_methods <- c("t.test", "wilcox.test", "anova")
  if (!statistical_method %in% valid_methods) {
    stop(sprintf("statistical_method must be one of: %s",
      paste(valid_methods, collapse = ", ")))
  }

  valid_plot_types <- c("box", "bar")
  if (!plot_type %in% valid_plot_types) {
    stop(sprintf("plot_type must be one of: %s",
      paste(valid_plot_types, collapse = ", ")))
  }
}

#' Calculate expression using original algorithm exactly
#' @keywords internal
calculate_qpcr_expression <- function(df) {

  df_expression <- df %>%
    dplyr::group_by(.data$bio_rep, .data$group, .data$gene) %>%
    dplyr::mutate(
      mean_cq = mean(.data$cq, na.rm = TRUE),
      sd_cq = stats::sd(.data$cq, na.rm = TRUE),
      sd_cq = ifelse(is.na(.data$sd_cq), 0, .data$sd_cq)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.data$bio_rep, .data$gene) %>%
    dplyr::mutate(
      min_mean_cq = min(.data$mean_cq),
      corrected_cq = .data$efficiency^(.data$min_mean_cq - .data$mean_cq),
      sd_corrected_cq = .data$sd_cq * .data$corrected_cq * log(.data$efficiency)
    ) %>%
    dplyr::ungroup()

  df_expression
}

#' Find reference genes using original GeNorm implementation
#' @keywords internal
find_reference_genes_original <- function(df_expression) {
  # Handle case with insufficient genes for GeNorm
  all_genes <- unique(df_expression$gene)
  if (length(all_genes) < 3) {
    # If less than 3 genes, return the first one or two genes
    return(all_genes[1:min(2, length(all_genes))])
  }

  tryCatch(
    {
      # Reproduce original GeNorm algorithm exactly
      df_ref <- df_expression %>%
        dplyr::select(.data$group, .data$gene, .data$cq, .data$bio_rep, .data$tech_rep) %>%
        dplyr::mutate(
          treatment = paste0(.data$group, .data$bio_rep, .data$tech_rep)
        ) %>%
        tidyr::pivot_wider(names_from = .data$gene, values_from = .data$cq)

      df_temp <- df_ref[, 5:ncol(df_ref)] %>% as.data.frame()

      # Remove columns with all NA values
      df_temp <- df_temp[, colSums(is.na(df_temp)) != nrow(df_temp)]

      if (ncol(df_temp) < 2) {
        return(colnames(df_temp))
      }

      n <- ncol(df_temp)
      M <- numeric(n)

      for (j in 1:n) {
        A <- log2(df_temp[, j] / df_temp[, -j, drop = FALSE])
        if (n > 2) {
          M[j] <- mean(apply(A, 2, stats::sd, na.rm = TRUE), na.rm = TRUE)
        } else {
          M[j] <- stats::sd(A, na.rm = TRUE)
        }
      }

      if (is.data.frame(df_temp)) {
        names(M) <- names(df_temp)
      } else {
        names(M) <- colnames(df_temp)
      }

      gene_symbol <- colnames(df_temp)
      n <- ncol(df_temp)
      num_ref <- min(2, n)  # Ensure we don't ask for more genes than available

      V <- numeric(max(0, n - num_ref))
      if (length(V) > 0) {
        names(V) <- paste(((n - 1):num_ref), "/", (n:(num_ref + 1)), sep = "")
      }

      mean_m <- numeric(n - num_ref + 1)
      names(mean_m) <- as.character(n:num_ref)
      R <- character(n)
      names(R) <- as.character(c(rep(1, num_ref), (num_ref + 1):length(R)))

      # Original geometric mean function
      geometric_mean <- function(x) {
        x <- x[!is.na(x)]
        if (length(x) == 0) return(1)
        if (any(x <= 0)) {
          x[x <= 0] <- min(x[x > 0], na.rm = TRUE) / 100  # Replace with small positive
        }
        return(prod(x)^(1 / length(x)))
      }

      # Original gene stability function
      gene_stable <- function(data, na.rm = TRUE) {
        if (!is.data.frame(data) & !is.matrix(data)) {
          stop("'data' has to of class matrix or data.frame")
        }
        n <- ncol(data)
        if (n == 1) return(setNames(0, colnames(data)[1]))  # Return single gene with 0 stability
        M <- numeric(n)
        for (j in 1:n) {
          A <- log2(data[, j] / data[, -j, drop = FALSE])
          if (n > 2) {
            M[j] <- mean(apply(A, 2, stats::sd, na.rm = na.rm), na.rm = TRUE)
          } else {
            M[j] <- stats::sd(A, na.rm = na.rm)
          }
        }
        if (is.data.frame(data)) {
          names(M) <- names(data)
        } else {
          names(M) <- colnames(data)
        }
        return(M)
      }

      # Original selection algorithm
      for (i in n:num_ref) {
        M <- gene_stable(df_temp, na.rm = TRUE)
        if (all(is.na(M))) break

        ind <- which.max(M)
        mean_m[n - i + 1] <- mean(M, na.rm = TRUE)
        if (i == num_ref) {
          R[1:num_ref] <- gene_symbol
        } else {
          R[i] <- gene_symbol[ind]
        }
        if (i > 2 && ncol(df_temp) > 1) {
          nf_old <- apply(df_temp, 1, geometric_mean)
          nf_new <- apply(df_temp[, -ind, drop = FALSE], 1, geometric_mean)
          V[n - i + 1] <- stats::sd(log2(nf_new / nf_old), na.rm = TRUE)
        }
        df_temp <- df_temp[, -ind, drop = FALSE]
        gene_symbol <- gene_symbol[-ind]
      }

      ref_gene <- as.character(R[1:num_ref])
      ref_gene <- ref_gene[!is.na(ref_gene) & ref_gene != ""]
      return(ref_gene)

    },
    error = function(e) {
      warning("GeNorm analysis failed: ", e$message, ". Using first available genes.")
      return(all_genes[1:min(2, length(all_genes))])
    })
}

#' Calculate normalization factors using original algorithm
#' @keywords internal
calculate_normalization_factors_original <- function(df_expression, reference_gene) {
  # Original geometric mean function
  geometric_mean <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(1)
    if (any(x <= 0)) {
      x[x <= 0] <- min(x[x > 0], na.rm = TRUE) / 100  # Replace with small positive
    }
    return(prod(x)^(1 / length(x)))
  }

  df_factor <- df_expression %>%
    dplyr::filter(.data$gene %in% reference_gene) %>%
    dplyr::mutate(temp_2 = paste0(.data$group, .data$bio_rep, .data$gene))

  df_factor <- df_factor[!duplicated(df_factor$temp_2), ] %>% as.data.frame()

  factor_df <- data.frame()

  for (i in unique(df_factor$bio_rep)) {
    df_temp <- df_factor %>% dplyr::filter(.data$bio_rep == i)
    for (j in unique(df_temp$group)) {
      df_temp_2 <- df_temp %>%
        dplyr::filter(.data$group == j) %>%
        dplyr::select(.data$group, .data$corrected_cq)

      if (nrow(df_temp_2) > 0) {
        fac <- data.frame(group = j, bio_rep = i, factor = geometric_mean(df_temp_2$corrected_cq))
        factor_df <- rbind(factor_df, fac)
      }
    }
  }

  if (nrow(factor_df) == 0) {
    stop("Unable to calculate normalization factors")
  }

  factor_df <- factor_df %>%
    dplyr::mutate(temp_2 = paste0(.data$group, .data$bio_rep)) %>%
    dplyr::select(.data$temp_2, .data$factor)

  df_factor <- df_factor %>%
    dplyr::mutate(temp_2 = paste0(.data$group, .data$bio_rep)) %>%
    merge(factor_df, by = "temp_2") %>%
    dplyr::mutate(sd_factor = (.data$sd_corrected_cq / (length(reference_gene) * (.data$corrected_cq)))^2) %>%
    dplyr::group_by(.data$bio_rep, .data$group) %>%
    dplyr::mutate(sd_factor = sqrt(sum(.data$sd_factor, na.rm = TRUE)) * .data$factor)

  df_factor
}

#' Calculate corrected expression using original algorithm
#' @keywords internal
calculate_corrected_expression_original <- function(df_expression, reference_gene, df_factor) {

  df_goi <- df_expression %>%
    dplyr::filter(!.data$gene %in% reference_gene) %>%
    dplyr::mutate(temp_2 = paste0(.data$group, .data$bio_rep)) %>%
    merge(df_factor[, c("temp_2", "factor", "sd_factor")], by = "temp_2") %>%
    dplyr::mutate(
      expression = .data$corrected_cq / .data$factor,
      sd_value = .data$expression * sqrt((.data$sd_corrected_cq / .data$corrected_cq)^2 + (.data$sd_factor / .data$factor)^2),
      se_value = .data$sd_value / sqrt(2)
    )

  # Calculate mean expression
  res_all <- df_goi %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.data$gene, .data$group) %>%
    dplyr::mutate(
      mean_expression = mean(unique(.data$expression), na.rm = TRUE),
      sd_expression = stats::sd(unique(.data$expression), na.rm = TRUE),
      se_expression = .data$sd_expression / sqrt(length(unique(.data$bio_rep)))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.data$gene) %>%
    dplyr::mutate(
      min_expression = min(.data$mean_expression, na.rm = TRUE),
      mean_expression = .data$mean_expression / .data$min_expression,
      sd_expression = .data$sd_expression / .data$min_expression,
      se_expression = .data$se_expression / .data$min_expression
    ) %>%
    dplyr::select(
      .data$group, .data$gene, .data$efficiency, .data$expression, .data$bio_rep,
      .data$mean_expression, .data$sd_expression, .data$se_expression
    ) %>%
    dplyr::rename(
      expression_value = .data$expression
    ) %>%
    dplyr::mutate(temp = paste0(.data$group, .data$gene, .data$bio_rep)) %>%
    dplyr::filter(!duplicated(.data$temp)) %>%
    dplyr::select(-.data$temp) %>%
    dplyr::mutate(temp = paste0(.data$gene, .data$group))

  res_all
}

#' Perform statistical analysis using original method
#' @keywords internal
perform_statistical_analysis_original <- function(res_all, reference_group, statistical_method) {

  tryCatch(
    {
      if (statistical_method == "t.test") {
        # Perform t-test for all non-reference groups
        stat_results <- res_all %>%
          dplyr::group_by(.data$gene) %>%
          rstatix::t_test(.data$expression_value ~ .data$group, ref.group = reference_group) %>%
          dplyr::ungroup() %>%
          dplyr::select(.data$gene, .data$group2, .data$p, .data$statistic, .data$df) %>%
          dplyr::mutate(
            p_value = round(.data$p, 6),
            t_statistic = round(.data$statistic, 4),
            significance = dplyr::case_when(
              .data$p < 0.001 ~ "***",
              .data$p > 0.001 & .data$p < 0.01 ~ "**",
              .data$p > 0.01 & .data$p < 0.05 ~ "*",
              TRUE ~ "NS"
            )
          ) %>%
          dplyr::rename(group = .data$group2) %>%
          dplyr::mutate(temp = paste0(.data$gene, .data$group)) %>%
          dplyr::select(.data$temp, .data$p_value, .data$t_statistic, .data$df, .data$significance)

        # Add empty statistics for reference group (no statistical annotation)
        ref_stat <- res_all %>%
          dplyr::filter(.data$group == reference_group) %>%
          dplyr::mutate(temp = paste0(.data$gene, .data$group)) %>%
          dplyr::select(.data$temp) %>%
          dplyr::distinct() %>%
          dplyr::mutate(
            p_value = NA_real_,
            t_statistic = NA_real_,
            df = NA_real_,
            significance = ""
          )

        # Combine all significance results
        all_stat <- rbind(stat_results, ref_stat)

        res_all <- res_all %>%
          dplyr::left_join(all_stat, by = "temp")

      } else if (statistical_method == "wilcox.test") {
        # Perform wilcox test for all non-reference groups
        stat_results <- res_all %>%
          dplyr::group_by(.data$gene) %>%
          rstatix::wilcox_test(.data$expression_value ~ .data$group, ref.group = reference_group) %>%
          dplyr::ungroup() %>%
          dplyr::select(.data$gene, .data$group2, .data$p, .data$statistic) %>%
          dplyr::mutate(
            p_value = round(.data$p, 6),
            w_statistic = round(.data$statistic, 4),
            significance = dplyr::case_when(
              .data$p < 0.001 ~ "***",
              .data$p > 0.001 & .data$p < 0.01 ~ "**",
              .data$p > 0.01 & .data$p < 0.05 ~ "*",
              TRUE ~ "NS"
            )
          ) %>%
          dplyr::rename(group = .data$group2) %>%
          dplyr::mutate(temp = paste0(.data$gene, .data$group)) %>%
          dplyr::select(.data$temp, .data$p_value, .data$w_statistic, .data$significance)

        # Add empty statistics for reference group (no statistical annotation)
        ref_stat <- res_all %>%
          dplyr::filter(.data$group == reference_group) %>%
          dplyr::mutate(temp = paste0(.data$gene, .data$group)) %>%
          dplyr::select(.data$temp) %>%
          dplyr::distinct() %>%
          dplyr::mutate(
            p_value = NA_real_,
            w_statistic = NA_real_,
            significance = ""
          )

        # Combine all significance results
        all_stat <- rbind(stat_results, ref_stat)

        res_all <- res_all %>%
          dplyr::left_join(all_stat, by = "temp")

      } else {
        # ANOVA with post-hoc Tukey test
        anova_results <- NULL
        df_stat <- NULL

        for (i in unique(res_all$gene)) {
          df_sub <- res_all %>%
            dplyr::filter(.data$gene == i) %>%
            dplyr::mutate(group = factor(.data$group))

          if (length(unique(df_sub$group)) > 1) {
            # ANOVA F-test
            fit <- stats::aov(.data$expression_value ~ .data$group, data = df_sub)
            anova_summary <- summary(fit)
            f_value <- anova_summary[[1]]["group", "F value"]
            p_value <- anova_summary[[1]]["group", "Pr(>F)"]
            df1 <- anova_summary[[1]]["group", "Df"]
            df2 <- anova_summary[[1]]["Residuals", "Df"]

            # Store ANOVA results for each gene
            anova_gene_result <- data.frame(
              gene = i,
              f_statistic = round(f_value, 4),
              p_value_anova = round(p_value, 6),
              df1 = df1,
              df2 = df2,
              stringsAsFactors = FALSE
            )
            anova_results <- rbind(anova_results, anova_gene_result)

            # Post-hoc Tukey test for letter annotations
            tuk <- multcomp::glht(fit, linfct = multcomp::mcp(group = "Tukey"))
            letters_result <- multcomp::cld(tuk, level = 0.95, decreasing = TRUE)

            letters_df <- letters_result[["mcletters"]][["Letters"]] %>%
              as.data.frame() %>%
              dplyr::mutate(gene = i) %>%
              tibble::rownames_to_column(var = "group") %>%
              magrittr::set_colnames(c("group", "significance", "gene")) %>%
              dplyr::select(.data$group, .data$gene, .data$significance) %>%
              dplyr::mutate(temp = paste0(.data$gene, .data$group)) %>%
              dplyr::select(.data$temp, .data$significance)

            # For ANOVA, make reference group letter empty
            letters_df <- letters_df %>%
              dplyr::mutate(significance = ifelse(grepl(paste0(reference_group, "$"), .data$temp), "", .data$significance))

            df_stat <- rbind(df_stat, letters_df)
          }
        }

        if (!is.null(df_stat) && !is.null(anova_results)) {
          # Add ANOVA statistics to each group
          res_all <- res_all %>%
            dplyr::left_join(anova_results, by = "gene") %>%
            dplyr::mutate(temp = paste0(.data$gene, .data$group)) %>%
            dplyr::left_join(df_stat, by = "temp")
        } else {
          res_all$significance <- ""
          res_all$f_statistic <- NA_real_
          res_all$p_value_anova <- NA_real_
          res_all$df1 <- NA_real_
          res_all$df2 <- NA_real_
        }
      }

      # Ensure reference group has empty significance
      res_all <- res_all %>%
        dplyr::mutate(significance = ifelse(.data$group == reference_group, "", .data$significance))

      return(res_all)

    },
    error = function(e) {
      warning("Statistical analysis failed: ", e$message)
      res_all$significance <- ""
      res_all$p_value <- NA_real_
      return(res_all)
    })
}

#' Create plot using original approach
#' @keywords internal
create_plot_original <- function(res_all, plot_type, plot_ncol) {

  if (nrow(res_all) == 0) {
    warning("No data available for plotting")
    return(NULL)
  }

  df_plot <- res_all %>%
    dplyr::rename(
      treatment = .data$group,
      expression = .data$expression_value,
      mean_expression = .data$mean_expression,
      sd_expression = .data$sd_expression,
      se_expression = .data$se_expression
    ) %>%
    dplyr::group_by(.data$gene, .data$treatment) %>%
    dplyr::mutate(n = dplyr::n()) %>%
    dplyr::ungroup()

  # Ensure significance column exists
  if (!"significance" %in% names(df_plot)) {
    df_plot$significance <- "NS"
  }

  if (plot_type == "box") {
    p <- df_plot %>%
      ggplot2::ggplot(ggplot2::aes(.data$treatment, .data$expression, fill = .data$treatment)) +
      ggplot2::geom_boxplot(width = 0.6) +
      ggplot2::facet_wrap(. ~ .data$gene, scales = "free_y", ncol = plot_ncol) +
      ggplot2::geom_text(ggplot2::aes(.data$treatment, min(.data$expression, na.rm = TRUE), label = .data$significance),
        check_overlap = TRUE, size = 3, color = "black"
      ) +
      ggthemes::theme_pander() +
      ggplot2::labs(y = "Relative expression") +
      ggplot2::theme(
        legend.position = "none",
        strip.text.x = ggplot2::element_text(face = "italic")
      )
  } else if (plot_type == "bar") {
    p <- df_plot %>%
      dplyr::group_by(.data$gene) %>%
      dplyr::mutate(max_temp = max(.data$mean_expression, na.rm = TRUE)) %>%
      dplyr::ungroup() %>%
      ggplot2::ggplot(ggplot2::aes(.data$treatment, .data$mean_expression / .data$n, fill = .data$treatment)) +
      ggplot2::geom_bar(stat = "identity", width = 0.6) +
      ggplot2::geom_errorbar(ggplot2::aes(.data$treatment,
        ymin = pmax(0, .data$mean_expression - .data$sd_expression),
        ymax = .data$mean_expression + .data$sd_expression
      ),
      width = 0.2
      ) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = .data$max_temp * 1.15), color = NA) +
      ggplot2::facet_wrap(. ~ .data$gene, scales = "free_y", ncol = plot_ncol) +
      ggplot2::geom_text(ggplot2::aes(.data$treatment, (.data$mean_expression + .data$sd_expression) * 1.08, label = .data$significance),
        check_overlap = TRUE, size = 4, color = "black"
      ) +
      ggthemes::theme_pander() +
      ggplot2::labs(y = "Relative expression") +
      ggplot2::theme(
        legend.position = "none",
        strip.text.x = ggplot2::element_text(face = "italic")
      )
  }

  p
}

# For backward compatibility - preserving exact original interface
CalExpRqPCR <- function(cq.table,
                        design.table,
                        ref.gene = NULL,
                        ref.group = "CK",
                        stat.method = "t.test",
                        fig.type = "box",
                        fig.ncol = NULL) {

  tryCatch(
    {
      result <- calc_expression_qpcr_efficiency(
        cq_table = cq.table,
        design_table = design.table,
        reference_gene = ref.gene,
        reference_group = ref.group,
        statistical_method = stat.method,
        plot_type = fig.type,
        plot_ncol = fig.ncol
      )

      return(result)

    },
    error = function(e) {
      warning("qPCR efficiency analysis failed: ", e$message)
      return(list(
        table = data.frame(),
        figure = NULL
      ))
    })
}


