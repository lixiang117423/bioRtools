#' Weighted Gene Co-expression Network Analysis (WGCNA)
#'
#' A comprehensive function to perform WGCNA including data preprocessing,
#' network construction, module identification, and module-trait association analysis.
#'
#' @param expression_data Either a file path to expression data (tab-delimited) or
#'   a data frame containing expression data. Required columns: sample ID, gene ID,
#'   and expression values.
#' @param sample_info Either a file path to sample information (Excel format) or
#'   a data frame containing sample metadata. Must include sample IDs and trait columns.
#' @param output_dir Character. Directory path where results and plots will be saved.
#' @param trait_columns Character vector. Column names in sample_info to use for
#'   module-trait association analysis.
#' @param gene_id_col Character. Column name containing gene IDs in expression_data.
#'   Default: "gene_id".
#' @param sample_col Character. Column name containing sample IDs in expression_data.
#'   Default: "sample".
#' @param expression_col Character. Column name containing expression values in
#'   expression_data. Default: "FPKM".
#' @param sample_id_col Character. Column name containing sample IDs in sample_info.
#'   Default: "treatment".
#' @param min_gene_sd Numeric. Minimum standard deviation threshold for gene filtering.
#'   Genes with SD <= this value will be excluded. Default: 0.1.
#' @param power_vector Integer vector. Powers to test for soft threshold selection.
#'   Default: 1:30.
#' @param min_module_size Integer. Minimum number of genes required for a module.
#'   Default: 50.
#' @param max_block_size Integer. Maximum block size for module detection. If NULL,
#'   uses all genes in one block. Default: NULL.
#' @param tom_type Character. Type of Topological Overlap Matrix. Options: "signed",
#'   "unsigned", "signed hybrid". Default: "signed".
#' @param correlation_method Character. Correlation method for module-trait association.
#'   Options: "pearson", "spearman", "kendall". Default: "pearson".
#' @param save_network Logical. Whether to save the network object as RData file.
#'   Default: TRUE.
#' @param plot_width Numeric. Width of saved plots in inches. Default: 12.
#' @param plot_height Numeric. Height of saved plots in inches. Default: 6.
#' @param plot_dpi Numeric. DPI resolution for saved plots. Default: 500.
#' @param verbose Logical. Whether to print progress messages. Default: TRUE.
#'
#' @return A list containing:
#' \describe{
#'   \item{network}{The WGCNA network object from blockwiseModules}
#'   \item{module_colors}{Character vector of module color assignments for each gene}
#'   \item{module_trait_cor}{Matrix of module-trait correlations}
#'   \item{module_trait_pvalue}{Matrix of module-trait correlation p-values}
#'   \item{soft_power}{Selected soft threshold power value}
#'   \item{gene_module_mapping}{Data frame mapping genes to modules with colors}
#'   \item{module_eigengenes}{Matrix of module eigengenes}
#'   \item{module_summary}{Data frame summarizing module sizes}
#'   \item{expression_matrix}{Processed expression matrix used for analysis}
#'   \item{trait_matrix}{Processed trait matrix used for correlation analysis}
#' }
#'
#' @details
#' This function performs a complete WGCNA pipeline:
#' \enumerate{
#'   \item Data loading and preprocessing with log2(x+1) transformation
#'   \item Quality control and gene filtering based on standard deviation
#'   \item Soft threshold power selection for scale-free topology
#'   \item Network construction using blockwiseModules
#'   \item Module identification and gene assignment
#'   \item Module-trait association analysis
#'   \item Generation of diagnostic plots and result files
#' }
#'
#' The function automatically creates the output directory if it doesn't exist and
#' generates multiple output files including correlation matrices, gene-module
#' mappings, and visualization plots.
#'
#' @note
#' Required packages: tidyverse, WGCNA, patchwork, readxl, pheatmap.
#' Gene IDs containing ":" will be split and the second part will be used.
#' Expression values are log2(x+1) transformed before analysis.
#'
#' @examples
#' library(bioRtools)
#' \dontrun{
#' # Using file paths
#' result1 <- run_wgcna_analysis(
#'   expression_data = "expression_data.txt",
#'   sample_info = "sample_info.xlsx",
#'   output_dir = "wgcna_results",
#'   trait_columns = c("time", "treatment")
#' )
#'
#' # Using data frames with custom column names
#' expr_df <- read.delim("expression.txt")
#' sample_df <- read_excel("samples.xlsx")
#'
#' result2 <- run_wgcna_analysis(
#'   expression_data = expr_df,
#'   sample_info = sample_df,
#'   output_dir = "wgcna_results",
#'   trait_columns = c("timepoint", "genotype"),
#'   gene_id_col = "transcript_id",
#'   expression_col = "TPM",
#'   sample_id_col = "sample_name"
#' )
#'
#' # Access results
#' print(result1$module_trait_cor)
#' head(result1$gene_module_mapping)
#' }
#'
#' @seealso
#' \code{\link[WGCNA]{blockwiseModules}}, \code{\link[WGCNA]{pickSoftThreshold}}
#'
#' @author Your Name
#' @export
run_wgcna_analysis <- function(
  expression_data,
  sample_info,
  output_dir,
  trait_columns,
  gene_id_col = "gene_id",
  sample_col = "sample",
  expression_col = "FPKM",
  sample_id_col = "treatment",
  min_gene_sd = 0.1,
  power_vector = 1:30,
  min_module_size = 50,
  max_block_size = NULL,
  tom_type = "signed",
  correlation_method = "pearson",
  save_network = TRUE,
  plot_width = 12,
  plot_height = 6,
  plot_dpi = 500,
  verbose = TRUE
) {
  # Input validation
  if (missing(expression_data)) {
    stop("'expression_data' is required")
  }
  if (missing(sample_info)) {
    stop("'sample_info' is required")
  }
  if (missing(output_dir)) {
    stop("'output_dir' is required")
  }
  if (missing(trait_columns)) {
    stop("'trait_columns' is required")
  }

  # Validate parameters
  if (!is.character(output_dir) || length(output_dir) != 1) {
    stop("'output_dir' must be a character string")
  }
  if (!is.character(trait_columns) || length(trait_columns) == 0) {
    stop("'trait_columns' must be a non-empty character vector")
  }
  if (!is.numeric(min_gene_sd) || min_gene_sd < 0) {
    stop("'min_gene_sd' must be a non-negative number")
  }
  if (!is.numeric(power_vector) || any(power_vector <= 0)) {
    stop("'power_vector' must contain positive numbers")
  }
  if (!is.numeric(min_module_size) || min_module_size <= 0) {
    stop("'min_module_size' must be a positive integer")
  }
  if (!tom_type %in% c("signed", "unsigned", "signed hybrid")) {
    stop("'tom_type' must be one of: 'signed', 'unsigned', 'signed hybrid'")
  }
  if (!correlation_method %in% c("pearson", "spearman", "kendall")) {
    stop("'correlation_method' must be one of: 'pearson', 'spearman', 'kendall'")
  }

  # Load required packages
  required_packages <- c("tidyverse", "WGCNA", "patchwork", "readxl", "pheatmap")
  missing_packages <- setdiff(required_packages, rownames(installed.packages()))

  if (length(missing_packages) > 0) {
    stop("Required packages missing: ", paste(missing_packages, collapse = ", "),
      "\nPlease install them first.")
  }

  # Load libraries quietly
  suppressPackageStartupMessages({
    library(tidyverse)
    library(WGCNA)
    library(patchwork)
    library(readxl)
    library(pheatmap)
  })

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    if (verbose) cat("Created output directory:", output_dir, "\n")
  }

  if (verbose) {
    cat("=== WGCNA Analysis Started ===\n")
    cat("Output directory:", output_dir, "\n")
    cat("Trait columns:", paste(trait_columns, collapse = ", "), "\n\n")
  }

  # 1. Data loading and preprocessing
  if (verbose) cat("1. Loading and preprocessing expression data...\n")

  # Load expression data
  if (is.character(expression_data)) {
    if (!file.exists(expression_data)) {
      stop("Expression file not found: ", expression_data)
    }
    df <- readr::read_delim(expression_data, delim = "\t", show_col_types = FALSE)
    if (verbose) cat("   Loaded expression data from file:", expression_data, "\n")
  } else if (is.data.frame(expression_data)) {
    df <- expression_data
    if (verbose) cat("   Using provided expression data frame\n")
  } else {
    stop("'expression_data' must be either a file path (character) or a data frame")
  }

  # Validate required columns in expression data
  required_expr_cols <- c(sample_col, gene_id_col, expression_col)
  missing_expr_cols <- setdiff(required_expr_cols, colnames(df))
  if (length(missing_expr_cols) > 0) {
    stop("Missing columns in expression_data: ", paste(missing_expr_cols, collapse = ", "))
  }

  # Process expression data
  df <- df %>%
    dplyr::select(all_of(c(sample_col, gene_id_col, expression_col))) %>%
    dplyr::rename(sample = !!sym(sample_col),
      gene_id = !!sym(gene_id_col),
      expression = !!sym(expression_col)) %>%
    dplyr::mutate(
      gene_id = ifelse(grepl(":", gene_id),
        sapply(strsplit(gene_id, ":"), "[", 2),
        gene_id)
    ) %>%
    dplyr::distinct_all() %>%
    dplyr::group_by(sample, gene_id) %>%
    dplyr::summarise(expression = sum(expression, na.rm = TRUE), .groups = "drop") %>%
    dplyr::filter(!is.na(expression), is.finite(expression), expression >= 0)

  if (verbose) {
    cat("   Processed expression data:", nrow(df), "entries for",
      length(unique(df$gene_id)), "genes and",
      length(unique(df$sample)), "samples\n")
  }

  # 2. Generate expression distribution plots
  if (verbose) cat("\n2. Generating expression distribution plots...\n")

  # Pre-transformation boxplot
  p1 <- df %>%
    ggplot(aes(x = sample, y = expression)) +
    geom_boxplot(fill = "lightblue", alpha = 0.7) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
    labs(title = "Expression Distribution (Before Log2 Transformation)",
      x = "Sample", y = "Expression Value")

  ggsave(p1, file = file.path(output_dir, "01_expression_distribution_before_log.png"),
    width = plot_width, height = plot_height, dpi = plot_dpi)

  # Post-transformation boxplot
  p2 <- df %>%
    ggplot(aes(x = sample, y = log2(expression + 1))) +
    geom_boxplot(fill = "lightcoral", alpha = 0.7) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
    labs(title = "Expression Distribution (After Log2 Transformation)",
      x = "Sample", y = "Log2(Expression + 1)")

  ggsave(p2, file = file.path(output_dir, "02_expression_distribution_after_log.png"),
    width = plot_width, height = plot_height, dpi = plot_dpi)

  # 3. Data filtering and quality control
  if (verbose) cat("\n3. Performing quality control and gene filtering...\n")

  # Apply log2 transformation and filter genes
  df_filtered <- df %>%
    dplyr::mutate(expression = log2(expression + 1)) %>%
    dplyr::group_by(gene_id) %>%
    dplyr::mutate(
      gene_sd = sd(expression, na.rm = TRUE),
      gene_mean = mean(expression, na.rm = TRUE)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(
      !is.na(gene_sd),
      is.finite(gene_sd),
      gene_sd > min_gene_sd,
      gene_sd > 0
    ) %>%
    dplyr::select(-gene_sd, -gene_mean)

  if (verbose) {
    cat("   After filtering: retained", length(unique(df_filtered$gene_id)), "genes\n")
  }

  # Convert to matrix format
  df_matrix <- df_filtered %>%
    tidyr::pivot_wider(names_from = gene_id, values_from = expression, values_fill = 0) %>%
    tibble::column_to_rownames(var = "sample")

  # Check data quality
  quality_check <- WGCNA::goodSamplesGenes(df_matrix, verbose = 0)
  if (verbose) {
    cat("   Quality check - Good samples:", sum(quality_check$goodSamples),
      "/ Good genes:", sum(quality_check$goodGenes), "\n")
  }

  # Remove bad samples/genes if any
  if (!quality_check$allOK) {
    if (sum(!quality_check$goodSamples) > 0) {
      df_matrix <- df_matrix[quality_check$goodSamples, ]
      if (verbose) cat("   Removed", sum(!quality_check$goodSamples), "problematic samples\n")
    }
    if (sum(!quality_check$goodGenes) > 0) {
      df_matrix <- df_matrix[, quality_check$goodGenes]
      if (verbose) cat("   Removed", sum(!quality_check$goodGenes), "problematic genes\n")
    }
  }

  # 4. Soft threshold selection
  if (verbose) cat("\n4. Selecting optimal soft threshold...\n")

  soft_threshold <- WGCNA::pickSoftThreshold(
    df_matrix,
    powerVector = power_vector,
    verbose = 0,
    blockSize = min(ncol(df_matrix), 5000)
  )

  selected_power <- soft_threshold$powerEstimate
  if (is.na(selected_power) || is.null(selected_power)) {
    # Find the first power where R^2 > 0.8 or use 6 as fallback
    fit_indices <- soft_threshold$fitIndices
    r_squared <- -sign(fit_indices$slope) * fit_indices$SFT.R.sq
    power_candidates <- which(r_squared > 0.8)

    if (length(power_candidates) > 0) {
      selected_power <- fit_indices$Power[power_candidates[1]]
    } else {
      selected_power <- 6
    }

    if (verbose) {
      cat("   Warning: Could not determine optimal power automatically.\n")
      cat("   Using power =", selected_power, "\n")
    }
  } else {
    if (verbose) cat("   Selected soft threshold power:", selected_power, "\n")
  }

  # Plot soft threshold selection
  soft_plot_data <- soft_threshold$fitIndices %>%
    as.data.frame() %>%
    dplyr::mutate(
      scale_free_r2 = -sign(slope) * SFT.R.sq,
      mean_connectivity = mean.k.
    )

  p3 <- soft_plot_data %>%
    ggplot(aes(x = Power, y = scale_free_r2)) +
    geom_point(size = 2, color = "red") +
    geom_line(alpha = 0.6) +
    geom_hline(yintercept = 0.85, color = "red", linetype = "dashed", alpha = 0.7) +
    geom_vline(xintercept = selected_power, color = "blue", linetype = "dashed", alpha = 0.7) +
    labs(x = "Soft Threshold (power)",
      y = expression("Scale Free Topology Model Fit, signed R"^2),
      title = "Scale Independence") +
    theme_minimal() +
    ylim(0, 1)

  p4 <- soft_plot_data %>%
    ggplot(aes(x = Power, y = mean_connectivity)) +
    geom_point(size = 2, color = "red") +
    geom_line(alpha = 0.6) +
    geom_vline(xintercept = selected_power, color = "blue", linetype = "dashed", alpha = 0.7) +
    labs(x = "Soft Threshold (power)",
      y = "Mean Connectivity",
      title = "Mean Connectivity") +
    theme_minimal()

  p_combined <- (p3 | p4) +
    patchwork::plot_annotation(
      title = "Soft Threshold Selection",
      tag_levels = "A"
    )

  ggsave(p_combined, file = file.path(output_dir, "03_soft_threshold_selection.png"),
    width = 10, height = 6, dpi = plot_dpi)

  # 5. Network construction
  if (verbose) cat("\n5. Constructing co-expression network...\n")

  # Set max block size
  if (is.null(max_block_size)) {
    max_block_size <- ncol(df_matrix)
  }

  # Temporarily override cor function to avoid conflicts
  cor_backup <- cor
  cor <- WGCNA::cor

  tryCatch({
    # Build network
    network <- WGCNA::blockwiseModules(
      df_matrix,
      power = selected_power,
      TOMType = tom_type,
      maxBlockSize = max_block_size,
      minModuleSize = min_module_size,
      mergeCutHeight = 0.25,
      numericLabels = TRUE,
      pamRespectsDendro = FALSE,
      saveTOMs = FALSE,
      verbose = 0
    )

    if (verbose) {
      cat("   Network construction completed\n")
      cat("   Number of modules identified:", length(unique(network$colors)) - 1, "\n")
    }

  }, error = function(e) {
    stop("Error in network construction: ", e$message)
  }, finally = {
    # Restore original cor function
    cor <- cor_backup
  })

  # Save network object
  if (save_network) {
    save(network, file = file.path(output_dir, "04_wgcna_network.RData"))
    if (verbose) cat("   Network object saved\n")
  }

  # Generate dynamic tree cut dendrogram
  if (verbose) cat("   Generating dynamic tree cut dendrogram...\n")

  tryCatch(
    {
      if (!is.null(network$dendrograms) && length(network$dendrograms) > 0) {
        # Create the dendrogram plot
        dendrogram_file <- file.path(output_dir, "04b_dynamic_tree_cut_dendrogram.png")

        png(filename = dendrogram_file,
          width = max(1000, ncol(df_matrix) * 2),
          height = 600,
          res = 150)

        # Plot dendrogram with module colors
        module_colors_plot <- WGCNA::labels2colors(network$colors)
        if (!is.null(network$blockGenes) && length(network$blockGenes) > 0) {
          colors_for_plot <- module_colors_plot[network$blockGenes[[1]]]
        } else {
          colors_for_plot <- module_colors_plot
        }

        WGCNA::plotDendroAndColors(
          dendro = network$dendrograms[[1]],
          colors = colors_for_plot,
          groupLabels = c("Dynamic Tree Cut"),
          dendroLabels = FALSE,
          hang = 0.03,
          addGuide = TRUE,
          guideHang = 0.05,
          main = "Gene Clustering and Module Assignment"
        )

        dev.off()

        if (verbose) cat("   Dendrogram saved to:", basename(dendrogram_file), "\n")

      } else {
        if (verbose) cat("   Warning: No dendrogram available for plotting\n")
      }
    },
    error = function(e) {
      if (verbose) cat("   Warning: Could not generate dendrogram plot:", e$message, "\n")
    })

  # 6. Module identification and gene assignment
  if (verbose) cat("\n6. Processing modules and gene assignments...\n")

  module_colors <- WGCNA::labels2colors(network$colors)

  # Create gene-module mapping
  gene_module_mapping <- data.frame(
    gene_id = names(network$colors),
    module_number = network$colors,
    module_color = module_colors,
    stringsAsFactors = FALSE
  )

  # Save gene-module mapping
  write.table(gene_module_mapping,
    file = file.path(output_dir, "05_gene_module_mapping.txt"),
    sep = "\t", row.names = FALSE, quote = FALSE)

  # Module size summary
  module_summary <- gene_module_mapping %>%
    dplyr::group_by(module_color, module_number) %>%
    dplyr::summarise(gene_count = n(), .groups = "drop") %>%
    dplyr::arrange(desc(gene_count))

  write.table(module_summary,
    file = file.path(output_dir, "06_module_size_summary.txt"),
    sep = "\t", row.names = FALSE, quote = FALSE)

  if (verbose) {
    cat("   Identified", nrow(module_summary), "modules (including grey)\n")
    cat("   Module sizes range from", min(module_summary$gene_count),
      "to", max(module_summary$gene_count), "genes\n")
  }

  # 7. Module-trait association analysis
  if (verbose) cat("\n7. Performing module-trait association analysis...\n")

  # Load sample information
  if (is.character(sample_info)) {
    if (!file.exists(sample_info)) {
      stop("Sample info file not found: ", sample_info)
    }
    sample_info_df <- readxl::read_excel(sample_info)
    if (verbose) cat("   Loaded sample info from file:", sample_info, "\n")
  } else if (is.data.frame(sample_info)) {
    sample_info_df <- sample_info
    if (verbose) cat("   Using provided sample info data frame\n")
  } else {
    stop("'sample_info' must be either a file path (character) or a data frame")
  }

  # Validate required columns in sample info
  required_sample_cols <- c(sample_id_col, trait_columns)
  missing_sample_cols <- setdiff(required_sample_cols, colnames(sample_info_df))
  if (length(missing_sample_cols) > 0) {
    stop("Missing columns in sample_info: ", paste(missing_sample_cols, collapse = ", "))
  }

  # Process sample information
  sample_info_processed <- sample_info_df %>%
    dplyr::select(all_of(c(sample_id_col, trait_columns))) %>%
    dplyr::rename(sample = !!sym(sample_id_col)) %>%
    dplyr::distinct_all() %>%
    dplyr::filter(sample %in% rownames(df_matrix))

  if (nrow(sample_info_processed) == 0) {
    stop("No matching samples between expression data and sample info")
  }

  # Create design matrix for traits
  trait_formula <- paste("~ 0 +", paste(trait_columns, collapse = " + "))

  tryCatch(
    {
      trait_matrix <- model.matrix(as.formula(trait_formula), data = sample_info_processed) %>%
        as.data.frame()
      rownames(trait_matrix) <- sample_info_processed$sample
    },
    error = function(e) {
      stop("Error creating trait matrix. Check trait columns for missing values or factor levels: ", e$message)
    })

  # Calculate module eigengenes
  module_eigengenes <- WGCNA::moduleEigengenes(df_matrix, module_colors, verbose = 0)$eigengenes %>%
    WGCNA::orderMEs()

  # Ensure samples match between eigengenes and traits
  common_samples <- intersect(rownames(module_eigengenes), rownames(trait_matrix))
  if (length(common_samples) == 0) {
    stop("No common samples between expression matrix and trait matrix")
  }

  module_eigengenes <- module_eigengenes[common_samples, , drop = FALSE]
  trait_matrix <- trait_matrix[common_samples, , drop = FALSE]

  # Calculate module-trait correlations
  module_trait_cor <- cor(module_eigengenes, trait_matrix,
    use = "pairwise.complete.obs", method = correlation_method)
  module_trait_pvalue <- WGCNA::corPvalueStudent(module_trait_cor, length(common_samples))

  # Save correlation results
  write.table(cbind(Module = rownames(module_trait_cor), as.data.frame(module_trait_cor)),
    file = file.path(output_dir, "07_module_trait_correlations.txt"),
    sep = "\t", row.names = FALSE, quote = FALSE)

  write.table(cbind(Module = rownames(module_trait_pvalue), as.data.frame(module_trait_pvalue)),
    file = file.path(output_dir, "08_module_trait_pvalues.txt"),
    sep = "\t", row.names = FALSE, quote = FALSE)

  if (verbose) {
    cat("   Module-trait correlation analysis completed\n")
    cat("   Analyzed", nrow(module_trait_cor), "modules against", ncol(module_trait_cor), "traits\n")
  }

  # 8. Generate module-trait heatmap
  if (verbose) cat("\n8. Generating module-trait association heatmap...\n")

  # Prepare data for plotting
  cor_plot_data <- module_trait_cor %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "module") %>%
    tidyr::pivot_longer(cols = -module, names_to = "trait", values_to = "correlation")

  pvalue_plot_data <- module_trait_pvalue %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "module") %>%
    tidyr::pivot_longer(cols = -module, names_to = "trait", values_to = "pvalue")

  plot_data <- cor_plot_data %>%
    dplyr::left_join(pvalue_plot_data, by = c("module", "trait")) %>%
    dplyr::mutate(
      significance = dplyr::case_when(
        pvalue < 0.001 ~ "***",
        pvalue < 0.01 ~ "**",
        pvalue < 0.05 ~ "*",
        TRUE ~ ""
      ),
      correlation_text = sprintf("%.2f", correlation)
    )

  # Create heatmap
  p5 <- plot_data %>%
    ggplot(aes(x = module, y = trait, fill = correlation)) +
    geom_tile(color = "white", size = 0.1) +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0,
      name = "Correlation",
      limits = c(-1, 1)
    ) +
    geom_text(aes(label = significance), size = 4, vjust = 0.7, fontface = "bold") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(title = "Module-Trait Associations",
      caption = "* p < 0.05, ** p < 0.01, *** p < 0.001") +
    coord_fixed()

  ggsave(p5, file = file.path(output_dir, "09_module_trait_heatmap.png"),
    width = plot_width, height = plot_height, dpi = plot_dpi)

  # 9. Summary and completion
  if (verbose) {
    cat("\n=== WGCNA Analysis Completed Successfully! ===\n")
    cat("Results saved to:", output_dir, "\n")
    cat("\nGenerated files:\n")
    cat("- 01_expression_distribution_before_log.png\n")
    cat("- 02_expression_distribution_after_log.png\n")
    cat("- 03_soft_threshold_selection.png\n")
    if (save_network) cat("- 04_wgcna_network.RData\n")
    cat("- 04b_dynamic_tree_cut_dendrogram.png\n")
    cat("- 05_gene_module_mapping.txt\n")
    cat("- 06_module_size_summary.txt\n")
    cat("- 07_module_trait_correlations.txt\n")
    cat("- 08_module_trait_pvalues.txt\n")
    cat("- 09_module_trait_heatmap.png\n")

    # Summary statistics
    sig_correlations <- sum(module_trait_pvalue < 0.05, na.rm = TRUE)
    cat("\nSummary:\n")
    cat("- Total modules:", nrow(module_summary), "\n")
    cat("- Total genes:", length(unique(gene_module_mapping$gene_id)), "\n")
    cat("- Significant module-trait correlations (p < 0.05):", sig_correlations, "\n")
  }

  # Return comprehensive results
  results <- list(
    network = network,
    module_colors = module_colors,
    module_trait_cor = module_trait_cor,
    module_trait_pvalue = module_trait_pvalue,
    soft_power = selected_power,
    gene_module_mapping = gene_module_mapping,
    module_eigengenes = module_eigengenes,
    module_summary = module_summary,
    expression_matrix = df_matrix,
    trait_matrix = trait_matrix
  )

  # Add metadata
  attr(results, "analysis_info") <- list(
    timestamp = Sys.time(),
    parameters = list(
      min_gene_sd = min_gene_sd,
      min_module_size = min_module_size,
      tom_type = tom_type,
      correlation_method = correlation_method,
      soft_power = selected_power
    ),
    sample_counts = list(
      total_samples = nrow(df_matrix),
      total_genes = ncol(df_matrix),
      total_modules = nrow(module_summary)
    )
  )

  class(results) <- c("wgcna_results", "list")

  return(results)
}

#' Print method for WGCNA results
#' @param x A wgcna_results object
#' @param ... Additional arguments (not used)
#' @author Xiang LI \email{lixiang117423@@foxmail.com}
#' @export
print.wgcna_results <- function(x, ...) {
  cat("WGCNA Analysis Results\n")
  cat("======================\n\n")

  info <- attr(x, "analysis_info")
  if (!is.null(info)) {
    cat("Analysis timestamp:", format(info$timestamp), "\n")
    cat("Soft threshold power:", info$parameters$soft_power, "\n")
    cat("Minimum module size:", info$parameters$min_module_size, "\n")
    cat("TOM type:", info$parameters$tom_type, "\n")
    cat("Correlation method:", info$parameters$correlation_method, "\n\n")

    cat("Data summary:\n")
    cat("- Samples:", info$sample_counts$total_samples, "\n")
    cat("- Genes:", info$sample_counts$total_genes, "\n")
    cat("- Modules:", info$sample_counts$total_modules, "\n\n")
  }

  cat("Available components:\n")
  cat("- network: WGCNA network object\n")
  cat("- module_colors: Module color assignments\n")
  cat("- module_trait_cor: Module-trait correlations\n")
  cat("- module_trait_pvalue: Module-trait p-values\n")
  cat("- gene_module_mapping: Gene-to-module assignments\n")
  cat("- module_eigengenes: Module eigengenes\n")
  cat("- expression_matrix: Processed expression data\n")
  cat("- trait_matrix: Processed trait data\n\n")

  sig_cors <- sum(x$module_trait_pvalue < 0.05, na.rm = TRUE)
  cat("Significant correlations (p < 0.05):", sig_cors, "\n")
}

#' Summary method for WGCNA results
#' @param object A wgcna_results object
#' @param ... Additional arguments (not used)
#' @author Xiang LI \email{lixiang117423@@foxmail.com}
#' @export
summary.wgcna_results <- function(object, ...) {
  cat("WGCNA Analysis Summary\n")
  cat("======================\n\n")

  # Module summary
  cat("Module Summary:\n")
  print(object$module_summary)
  cat("\n")

  # Top correlations
  cor_data <- as.data.frame(object$module_trait_cor) %>%
    tibble::rownames_to_column("Module") %>%
    tidyr::pivot_longer(-Module, names_to = "Trait", values_to = "Correlation") %>%
    dplyr::arrange(desc(abs(Correlation))) %>%
    dplyr::slice_head(n = 10)

  cat("Top 10 Module-Trait Correlations:\n")
  print(cor_data)
}
