#' Calculate Gene Expression Using qPCR with Efficiency Correction
#'
#' @description
#' Calculate expression using qPCR with efficiency correction and automatic
#' reference gene selection using GeNorm algorithm if reference genes are not specified.
#' This function preserves the original algorithm logic from the author.
#'
#' @param cq_table A data frame containing qPCR data.
#'   Must contain columns: Position, Cq, Gene
#' @param design_table A data frame containing experimental design information.
#'   Must contain columns: Position, Group, BioRep, TechRep, Eff
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
#' @importFrom tidyr spread
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
#' # Load example data
#' cq_data_path <- system.file("examples", "cal.expre.rqpcr.cq.txt", package = "qPCRtools")
#' design_data_path <- system.file("examples", "cal.expre.rqpcr.design.txt", package = "qPCRtools")
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
  
  # Input validation
  validate_qpcr_inputs(cq_table, design_table, statistical_method, plot_type)
  
  # Merge data - preserving original algorithm logic
  df <- merge_qpcr_data(cq_table, design_table)
  
  # Calculate expression using original algorithm
  df_expression <- calculate_qpcr_expression(df)
  
  # Find reference gene using original GeNorm algorithm if not provided
  if (is.null(reference_gene)) {
    reference_gene <- find_reference_genes_original(df_expression)
  }
  
  # Calculate normalization factors using original algorithm
  normalization_factors <- calculate_normalization_factors_original(df_expression, reference_gene)
  
  # Calculate corrected expression using original algorithm
  results_all <- calculate_corrected_expression_original(df_expression, reference_gene, normalization_factors)
  
  # Perform statistical analysis using original method
  results_with_stats <- perform_statistical_analysis_original(results_all, reference_group, statistical_method)
  
  # Create plot using original approach
  plot_result <- create_plot_original(results_with_stats, plot_type, plot_ncol)
  
  # Clean up results to match original format
  final_results <- results_with_stats %>%
    dplyr::select(-.data$temp, -.data$eff) %>%
    dplyr::select(1, 4, 2, 3, 5:7)
  
  return(list(table = final_results, figure = plot_result))
}

#' Validate input parameters - keeping original requirements
#' @keywords internal
validate_qpcr_inputs <- function(cq_table, design_table, statistical_method, plot_type) {
  
  if (!is.data.frame(cq_table)) {
    stop("cq_table must be a data frame")
  }
  
  if (!is.data.frame(design_table)) {
    stop("design_table must be a data frame")
  }
  
  # Check required columns exactly as original
  required_cq_cols <- c("Position", "Cq", "Gene")
  missing_cq_cols <- setdiff(required_cq_cols, names(cq_table))
  if (length(missing_cq_cols) > 0) {
    stop(sprintf("cq_table missing required columns: %s", 
                paste(missing_cq_cols, collapse = ", ")))
  }
  
  required_design_cols <- c("Position", "Group", "BioRep", "TechRep", "Eff")
  missing_design_cols <- setdiff(required_design_cols, names(design_table))
  if (length(missing_design_cols) > 0) {
    stop(sprintf("design_table missing required columns: %s", 
                paste(missing_design_cols, collapse = ", ")))
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

#' Merge data exactly as in original
#' @keywords internal
merge_qpcr_data <- function(cq_table, design_table) {
  
  df <- cq_table %>%
    dplyr::left_join(design_table, by = "Position") %>%
    dplyr::rename(
      position = .data$Position,
      cq = .data$Cq,
      group = .data$Group,
      gene = .data$Gene,
      biorep = .data$BioRep,
      techrep = .data$TechRep,
      eff = .data$Eff
    )
  
  return(df)
}

#' Calculate expression using original algorithm exactly
#' @keywords internal
calculate_qpcr_expression <- function(df) {
  
  df_expression <- df %>%
    dplyr::group_by(.data$biorep, .data$group, .data$gene) %>%
    dplyr::mutate(
      mean.cq = mean(.data$cq, na.rm = TRUE),
      sd.cq = stats::sd(.data$cq, na.rm = TRUE),
      sd.cq = ifelse(is.na(.data$sd.cq), 0, .data$sd.cq)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.data$biorep, .data$gene) %>%
    dplyr::mutate(
      min.mean.cq = min(.data$mean.cq),
      QCq = .data$eff^(.data$min.mean.cq - .data$mean.cq),
      SD_QCq = .data$sd.cq * .data$QCq * log(.data$eff)
    ) %>%
    dplyr::ungroup()
  
  return(df_expression)
}

#' Find reference genes using original GeNorm implementation
#' @keywords internal
find_reference_genes_original <- function(df_expression) {
  
  # Reproduce original GeNorm algorithm exactly
  df_ref <- df_expression %>%
    dplyr::select(.data$group, .data$gene, .data$cq, .data$biorep, .data$techrep) %>%
    dplyr::mutate(
      Treatment = paste0(.data$group, .data$biorep, .data$techrep)
    ) %>%
    tidyr::spread(key = .data$gene, value = .data$cq)
  
  df_temp <- df_ref[, 5:ncol(df_ref)] %>% as.data.frame()
  
  n <- length(unique(df_expression$gene))
  M <- numeric(n)
  
  for (j in 1:n) {
    A <- log2(df_temp[, j] / df_temp[, -j])
    if (n > 2) {
      M[j] <- mean(apply(A, 2, stats::sd, na.rm = TRUE))
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
  num_ref <- 2
  
  V <- numeric(n - num_ref)
  names(V) <- paste(((n - 1):num_ref), "/", (n:(num_ref + 1)), sep = "")
  mean_M <- numeric(n - num_ref + 1)
  names(mean_M) <- as.character(n:num_ref)
  R <- character(n)
  names(R) <- as.character(c(rep(1, num_ref), (num_ref + 1):length(R)))
  
  # Original geometric mean function
  geometric_mean <- function(x) {
    x <- x[!is.na(x)]
    if (any(x < 0)) {
      stop("'x' contains negative value(s)")
    } else {
      return(prod(x)^(1 / length(x)))
    }
  }
  
  # Original gene stability function
  gene_stable <- function(data, na.rm = TRUE) {
    if (!is.data.frame(data) & !is.matrix(data)) {
      stop("'data' has to of class matrix or data.frame")
    }
    n <- ncol(data)
    if (n == 1) stop("you need at least two genes for this computation")
    M <- numeric(n)
    for (j in 1:n) {
      A <- log2(data[, j] / data[, -j])
      if (n > 2) {
        M[j] <- mean(apply(A, 2, stats::sd, na.rm = na.rm))
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
    ind <- which.max(M)
    mean_M[n - i + 1] <- mean(M)
    if (i == num_ref) {
      R[1:num_ref] <- gene_symbol
    } else {
      R[i] <- gene_symbol[ind]
    }
    if (i > 2) {
      NF_old <- apply(df_temp, 1, geometric_mean)
      NF_new <- apply(df_temp[, -ind], 1, geometric_mean)
      V[n - i + 1] <- stats::sd(log2(NF_new / NF_old), na.rm = TRUE)
    }
    df_temp <- df_temp[, -ind]
    gene_symbol <- gene_symbol[-ind]
  }
  
  ref_gene <- as.character(R[1:num_ref])
  return(ref_gene)
}

#' Calculate normalization factors using original algorithm
#' @keywords internal
calculate_normalization_factors_original <- function(df_expression, reference_gene) {
  
  # Original geometric mean function
  geometric_mean <- function(x) {
    x <- x[!is.na(x)]
    if (any(x < 0)) {
      stop("'x' contains negative value(s)")
    } else {
      return(prod(x)^(1 / length(x)))
    }
  }
  
  df_factor <- df_expression %>%
    dplyr::filter(.data$gene %in% reference_gene) %>%
    dplyr::mutate(temp_2 = paste0(.data$group, .data$biorep, .data$gene))
  
  df_factor <- df_factor[!duplicated(df_factor$temp_2), ] %>% as.data.frame()
  
  factor_df <- data.frame()
  
  for (i in unique(df_factor$biorep)) {
    df_temp <- df_factor %>% dplyr::filter(.data$biorep == i)
    for (j in unique(df_temp$group)) {
      df_temp_2 <- df_temp %>%
        dplyr::filter(.data$group == j) %>%
        dplyr::select(.data$group, .data$QCq)
      fac <- data.frame(group = j, biorep = i, factor = geometric_mean(df_temp_2$QCq))
      factor_df <- rbind(factor_df, fac)
    }
  }
  
  factor_df <- factor_df %>%
    dplyr::mutate(temp_2 = paste0(.data$group, .data$biorep)) %>%
    dplyr::select(.data$temp_2, .data$factor)
  
  df_factor <- df_factor %>%
    dplyr::mutate(temp_2 = paste0(.data$group, .data$biorep)) %>%
    merge(factor_df, by = "temp_2") %>%
    dplyr::mutate(SD.factor = (.data$SD_QCq / (length(reference_gene) * (.data$QCq)))^2) %>%
    dplyr::group_by(.data$biorep, .data$group) %>%
    dplyr::mutate(SD.factor = sqrt(sum(.data$SD.factor)) * .data$factor)
  
  return(df_factor)
}

#' Calculate corrected expression using original algorithm
#' @keywords internal
calculate_corrected_expression_original <- function(df_expression, reference_gene, df_factor) {
  
  df_goi <- df_expression %>%
    dplyr::filter(!.data$gene %in% reference_gene) %>%
    dplyr::mutate(temp_2 = paste0(.data$group, .data$biorep)) %>%
    merge(df_factor[, c("temp_2", "factor", "SD.factor")], by = "temp_2") %>%
    dplyr::mutate(
      expression = .data$QCq / .data$factor,
      SD_1 = .data$expression * sqrt((.data$SD_QCq / .data$QCq)^2 + (.data$SD.factor / .data$factor)^2),
      SE_1 = .data$SD_1 / sqrt(2)
    )
  
  # Calculate mean expression exactly as original
  res_all <- df_goi %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.data$gene, .data$group) %>%
    dplyr::mutate(
      mean.expression = mean(unique(.data$expression), na.rm = TRUE),
      sd.expression = stats::sd(unique(.data$expression), na.rm = TRUE),
      se.expression = .data$sd.expression / sqrt(length(unique(.data$biorep)))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.data$gene) %>%
    dplyr::mutate(
      min.expression = min(.data$mean.expression),
      mean.expression = .data$mean.expression / .data$min.expression,
      sd.expression = .data$sd.expression / .data$min.expression,
      se.expression = .data$se.expression / .data$min.expression
    ) %>%
    dplyr::select(
      .data$group, .data$gene, .data$eff, .data$expression, .data$biorep,
      .data$mean.expression, .data$sd.expression, .data$se.expression
    ) %>%
    dplyr::rename(
      Expre4Stat = .data$expression,
      Expression = .data$mean.expression,
      SD = .data$sd.expression,
      SE = .data$se.expression
    ) %>%
    dplyr::mutate(temp = paste0(.data$group, .data$gene, .data$biorep)) %>%
    dplyr::filter(!duplicated(.data$temp)) %>%
    dplyr::select(-.data$temp) %>%
    dplyr::mutate(temp = paste0(.data$gene, .data$group))
  
  return(res_all)
}

#' Perform statistical analysis using original method
#' @keywords internal
perform_statistical_analysis_original <- function(res_all, reference_group, statistical_method) {
  
  if (statistical_method == "t.test") {
    res_all %>%
      dplyr::group_by(.data$gene) %>%
      rstatix::t_test(.data$Expre4Stat ~ .data$group, ref.group = reference_group) %>%
      dplyr::ungroup() %>%
      dplyr::select(.data$gene, .data$group2, .data$p) %>%
      dplyr::mutate(signif = dplyr::case_when(
        .data$p < 0.001 ~ "***",
        .data$p > 0.001 & .data$p < 0.01 ~ "**",
        .data$p > 0.01 & .data$p < 0.05 ~ "*",
        TRUE ~ "NS"
      )) %>%
      dplyr::add_row(group2 = reference_group, p = NA, signif = NA) %>%
      dplyr::rename(group = .data$group2) %>%
      dplyr::mutate(temp = paste0(.data$gene, .data$group)) %>%
      dplyr::select(.data$temp, .data$signif) -> df_stat
    
    res_all <- res_all %>%
      dplyr::left_join(df_stat, by = "temp")
    
  } else if (statistical_method == "wilcox.test") {
    res_all %>%
      dplyr::group_by(.data$gene) %>%
      rstatix::wilcox_test(.data$Expre4Stat ~ .data$group, ref.group = reference_group) %>%
      dplyr::ungroup() %>%
      dplyr::select(.data$gene, .data$group2, .data$p) %>%
      dplyr::mutate(signif = dplyr::case_when(
        .data$p < 0.001 ~ "***",
        .data$p > 0.001 & .data$p < 0.01 ~ "**",
        .data$p > 0.01 & .data$p < 0.05 ~ "*",
        TRUE ~ "NS"
      )) %>%
      dplyr::add_row(group2 = reference_group, p = NA, signif = NA) %>%
      dplyr::rename(group = .data$group2) %>%
      dplyr::mutate(temp = paste0(.data$gene, .data$group)) %>%
      dplyr::select(.data$temp, .data$signif) -> df_stat
    
    res_all <- res_all %>%
      dplyr::left_join(df_stat, by = "temp")
    
  } else {
    df_stat <- NULL
    for (i in unique(res_all$gene)) {
      df_sub <- res_all %>%
        dplyr::filter(.data$gene == i) %>%
        dplyr::mutate(group = factor(.data$group))
      fit <- stats::aov(.data$Expre4Stat ~ .data$group, data = df_sub)
      tuk <- multcomp::glht(fit, linfct = multcomp::mcp(group = "Tukey"))
      multcomp::cld(tuk, level = 0.95, ddecreasing = TRUE)[["mcletters"]][["Letters"]] %>%
        as.data.frame() %>%
        dplyr::mutate(gene = i) %>%
        tibble::rownames_to_column(var = "group") %>%
        magrittr::set_colnames(c("group", "signif", "gene")) %>%
        dplyr::select(.data$group, .data$gene, .data$signif) %>%
        dplyr::mutate(temp = paste0(.data$group, .data$gene)) %>%
        dplyr::select(.data$temp, .data$signif) %>%
        rbind(df_stat) -> df_stat
    }
    res_all %>%
      dplyr::mutate(temp = paste0(.data$group, .data$gene)) %>%
      dplyr::left_join(df_stat, by = "temp") -> res_all
  }
  
  return(res_all)
}

#' Create plot using original approach
#' @keywords internal
create_plot_original <- function(res_all, plot_type, plot_ncol) {
  
  df_plot <- res_all %>%
    dplyr::rename(
      Treatment = .data$group,
      gene = .data$gene,
      expre = .data$Expre4Stat,
      mean.expre = .data$Expression,
      sd.expre = .data$SD,
      se.expre = .data$SE
    ) %>%
    dplyr::group_by(.data$gene, .data$Treatment) %>%
    dplyr::mutate(n = dplyr::n()) %>%
    dplyr::ungroup()
  
  if (plot_type == "box") {
    p <- df_plot %>%
      ggplot2::ggplot(ggplot2::aes(.data$Treatment, .data$expre, fill = .data$Treatment)) +
      ggplot2::geom_boxplot(width = 0.6) +
      ggplot2::facet_wrap(. ~ .data$gene, scales = "free_y", ncol = plot_ncol) +
      ggplot2::geom_text(ggplot2::aes(.data$Treatment, min(.data$expre), label = .data$signif),
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
      dplyr::mutate(max.temp = max(.data$mean.expre)) %>%
      dplyr::ungroup() %>%
      ggplot2::ggplot(ggplot2::aes(.data$Treatment, .data$mean.expre / .data$n, fill = .data$Treatment)) +
      ggplot2::geom_bar(stat = "identity", width = 0.6) +
      ggplot2::geom_errorbar(ggplot2::aes(.data$Treatment,
        ymin = .data$mean.expre - .data$sd.expre,
        ymax = .data$mean.expre + .data$sd.expre
      ),
      width = 0.2
      ) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = .data$max.temp * 1.15), color = NA) +
      ggplot2::facet_wrap(. ~ .data$gene, scales = "free_y", ncol = plot_ncol) +
      ggplot2::geom_text(ggplot2::aes(.data$Treatment, (.data$mean.expre + .data$sd.expre) * 1.08, label = .data$signif),
        check_overlap = TRUE, size = 4, color = "black"
      ) +
      ggthemes::theme_pander() +
      ggplot2::labs(y = "Relative expression") +
      ggplot2::theme(
        legend.position = "none",
        strip.text.x = ggplot2::element_text(face = "italic")
      )
  }
  
  return(p)
}

# For backward compatibility - preserving exact original interface
CalExpRqPCR <- function(cq.table,
                        design.table,
                        ref.gene = NULL,
                        ref.group = "CK",
                        stat.method = "t.test",
                        fig.type = "box",
                        fig.ncol = NULL) {
  
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
}