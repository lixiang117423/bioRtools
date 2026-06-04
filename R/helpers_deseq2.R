#' DESeq2 Analysis Workflow Helper Functions
#'
#' Internal helper functions for DESeq2-based differential analysis workflows.
#' These functions encapsulate common DESeq2 steps used in find_degs_deseq2 and
#' find_dams_deseq2, reducing code duplication and improving maintainability.
#'
#' @keywords internal
#' @noRd

# ============================================================================
# DESeq2 Data Preparation
# ============================================================================

#' Create DESeqDataSet object
#'
#' @param count_data Count matrix (genes/features × samples)
#' @param sample_metadata Sample metadata data frame with rownames matching colnames of count_data
#' @param design_formula Design formula for DESeq2 model
#'
#' @return DESeqDataSet object
#' @keywords internal
#' @noRd
create_deseq2_dataset <- function(count_data, sample_metadata, design_formula) {
  tryCatch(
    {
      DESeq2::DESeqDataSetFromMatrix(
        countData = count_data,
        colData = sample_metadata,
        design = design_formula
      )
    },
    error = function(e) {
      stop(err_deseq2_workflow("DESeqDataSet creation", e$message))
    }
  )
}

#' Run DESeq2 normalization and dispersion estimation
#'
#' @param dds DESeqDataSet object
#'
#' @return Analyzed DESeqDataSet object
#' @keywords internal
#' @noRd
run_deseq2_analysis <- function(dds) {
  tryCatch(
    {
      DESeq2::DESeq(dds, quiet = TRUE)
    },
    error = function(e) {
      stop(err_deseq2_workflow("DESeq2 fitting", e$message))
    }
  )
}

# ============================================================================
# DESeq2 Results Extraction
# ============================================================================

#' Extract and optionally shrink log2 fold changes
#'
#' @param dds_analyzed Analyzed DESeqDataSet object
#' @param alpha Alpha level for results (default: 0.1)
#' @param shrink_lfc Apply LFC shrinkage (default: TRUE)
#' @param independent_filtering Apply independent filtering (default: TRUE)
#'
#' @return Data frame of results
#' @keywords internal
#' @noRd
extract_deseq2_results <- function(dds_analyzed, alpha = 0.1, shrink_lfc = TRUE,
                                   independent_filtering = TRUE) {
  results_df <- tryCatch(
    {
      if (shrink_lfc) {
        message("Applying LFC shrinkage with apeglm method...")
        results_raw <- DESeq2::lfcShrink(
          dds_analyzed,
          coef = DESeq2::resultsNames(dds_analyzed)[
            length(DESeq2::resultsNames(dds_analyzed))
          ],
          type = "apeglm",
          quiet = TRUE
        )
        message("LFC shrinkage successful.")
      } else {
        results_raw <- DESeq2::results(
          dds_analyzed,
          alpha = alpha,
          independentFiltering = independent_filtering
        )
      }

      as.data.frame(results_raw)
    },
    error = function(e) {
      warn_data_processing(
        paste("LFC shrinkage failed, using unshrunken results:", e$message)
      )

      # Fallback: extract without shrinkage
      results_raw <- DESeq2::results(
        dds_analyzed,
        alpha = alpha,
        independentFiltering = independent_filtering
      )
      as.data.frame(results_raw)
    }
  )

  results_df
}

# ============================================================================
# Results Classification and Processing
# ============================================================================

#' Classify regulation direction and add annotation columns
#'
#' @param results_df Results data frame from DESeq2
#' @param fc_threshold Log2 fold change threshold (default: 1)
#' @param padj_threshold Adjusted p-value threshold (default: 0.05)
#' @param regulation_labels Named character vector for regulation classifications
#'   (default: list(up = "Up-regulated", down = "Down-regulated", ns = "Not significant"))
#' @param feature_id_col Column name for feature IDs (default: auto-detect from rownames)
#'
#' @return Annotated results data frame with regulation classifications
#' @keywords internal
#' @noRd
classify_deseq2_results <- function(results_df, fc_threshold = 1, padj_threshold = 0.05,
                                    regulation_labels = list(
                                      up = "Up-regulated",
                                      down = "Down-regulated",
                                      ns = "Not significant"
                                    ),
                                    feature_id_col = NULL) {
  # Auto-generate feature ID column if not provided
  if (is.null(feature_id_col)) {
    feature_id_col <- if (!is.null(rownames(results_df))) "feature_id" else "feature"
  }

  results_processed <- results_df %>%
    tibble::rownames_to_column(var = feature_id_col) %>%
    dplyr::mutate(
      # Main regulation classification
      regulation = dplyr::case_when(
        log2FoldChange > !!fc_threshold & padj < !!padj_threshold ~
          !!regulation_labels$up,
        log2FoldChange < -!!fc_threshold & padj < !!padj_threshold ~
          !!regulation_labels$down,
        TRUE ~ !!regulation_labels$ns
      ),
      # Fold change conversion
      fold_change = 2^log2FoldChange,
      abs_log2fc = abs(log2FoldChange),
      # Handle missing adjusted p-values
      padj = ifelse(is.na(padj), 1, padj),
      # Detailed significance levels
      significance_level = dplyr::case_when(
        padj >= 0.05 ~ "NS (p≥0.05)",
        padj >= 0.01 ~ "* (p<0.05)",
        padj >= 0.001 ~ "** (p<0.01)",
        padj >= 0.0001 ~ "*** (p<0.001)",
        TRUE ~ "**** (p<0.0001)"
      )
    ) %>%
    # Sort by significance and effect size
    dplyr::arrange(padj, desc(abs_log2fc))

  results_processed
}

# ============================================================================
# Results Summary and Statistics
# ============================================================================

#' Calculate summary statistics from DESeq2 results
#'
#' @param results_processed Processed results data frame with 'regulation' column
#' @param results_raw Raw results data frame before processing
#' @param count_matrix Original count matrix
#' @param design_formula Design formula used
#' @param regulation_labels Regulation label names (for counting)
#'
#' @return List with summary statistics
#' @keywords internal
#' @noRd
calculate_deseq2_summary <- function(results_processed, results_raw, count_matrix,
                                     design_formula, regulation_labels = list(
                                       up = "Up-regulated",
                                       down = "Down-regulated"
                                     )) {
  n_up <- sum(results_processed$regulation == regulation_labels$up, na.rm = TRUE)
  n_down <- sum(results_processed$regulation == regulation_labels$down, na.rm = TRUE)
  n_total_de <- n_up + n_down
  n_tested <- sum(!is.na(results_processed$padj))

  list(
    n_features_input = nrow(count_matrix),
    n_features_tested = n_tested,
    n_samples = ncol(count_matrix),
    n_up = n_up,
    n_down = n_down,
    n_total_de = n_total_de,
    design_formula = deparse(design_formula),
    mean_basemean = mean(results_processed$baseMean, na.rm = TRUE),
    median_basemean = median(results_processed$baseMean, na.rm = TRUE),
    pct_significant = (n_total_de / n_tested) * 100
  )
}

#' Attach summary statistics as attributes
#'
#' @param results_df Results data frame
#' @param summary_stats List of summary statistics
#'
#' @return Results data frame with attributes
#' @keywords internal
#' @noRd
attach_deseq2_metadata <- function(results_df, summary_stats) {
  for (name in names(summary_stats)) {
    attr(results_df, paste0("deseq2_", name)) <- summary_stats[[name]]
  }

  results_df
}

# ============================================================================
# Quality Control and Diagnostics
# ============================================================================

#' Check library size consistency
#'
#' @param count_matrix Count matrix
#' @param warn_ratio Maximum acceptable ratio between largest and smallest library size
#'   (default: 3)
#'
#' @return Invisibly returns TRUE; issues warning if ratio exceeds threshold
#' @keywords internal
#' @noRd
check_library_size_ratio <- function(count_matrix, warn_ratio = 3) {
  lib_sizes <- colSums(count_matrix)
  lib_size_ratio <- max(lib_sizes) / min(lib_sizes)

  if (lib_size_ratio > warn_ratio) {
    warn_data_processing(
      sprintf(
        "Large variation in library sizes (%.1fx). Results may be influenced by sequencing depth differences.",
        lib_size_ratio
      )
    )
  }

  if (min(lib_sizes) < 1000) {
    warn_data_processing(
      "Very small library sizes detected (<1000 counts). Results may be unreliable."
    )
  }

  invisible(TRUE)
}

#' Identify low-abundance features for potential filtering
#'
#' @param count_matrix Count matrix
#' @param min_count_threshold Minimum count per feature (default: 10)
#'
#' @return Character vector of feature names to consider removing
#' @keywords internal
#' @noRd
identify_low_abundance_features <- function(count_matrix, min_count_threshold = 10) {
  feature_sums <- rowSums(count_matrix)
  low_abundance_idx <- which(feature_sums < min_count_threshold)

  if (length(low_abundance_idx) > 0) {
    rownames(count_matrix)[low_abundance_idx]
  } else {
    character(0)
  }
}

# ============================================================================
# Convenience Wrappers
# ============================================================================

#' Complete DESeq2 workflow from count matrix to results
#'
#' Combines all DESeq2 steps into a single function for common use cases.
#'
#' @param count_data Count matrix
#' @param sample_metadata Sample metadata
#' @param design_formula Design formula
#' @param fc_threshold Log2 fold change threshold
#' @param padj_threshold Adjusted p-value threshold
#' @param shrink_lfc Apply LFC shrinkage
#' @param independent_filtering Use independent filtering
#' @param verbose Print progress messages
#'
#' @return List with:
#'   - results: Processed results data frame
#'   - summary: Summary statistics
#'   - dds: Analyzed DESeqDataSet object (for additional inspection)
#'
#' @keywords internal
#' @noRd
run_complete_deseq2_workflow <- function(count_data, sample_metadata, design_formula,
                                         fc_threshold = 1, padj_threshold = 0.05,
                                         shrink_lfc = TRUE, independent_filtering = TRUE,
                                         verbose = TRUE) {
  if (verbose) cat("Creating DESeq2 dataset...\n")
  dds <- create_deseq2_dataset(count_data, sample_metadata, design_formula)

  if (verbose) cat("Running DESeq2 analysis...\n")
  dds_analyzed <- run_deseq2_analysis(dds)

  if (verbose) cat("Extracting results...\n")
  results_raw <- extract_deseq2_results(dds_analyzed,
    shrink_lfc = shrink_lfc,
    independent_filtering = independent_filtering
  )

  if (verbose) cat("Classifying results...\n")
  results_processed <- classify_deseq2_results(
    results_raw,
    fc_threshold = fc_threshold,
    padj_threshold = padj_threshold
  )

  if (verbose) cat("Calculating summary statistics...\n")
  summary_stats <- calculate_deseq2_summary(
    results_processed, results_raw, count_data,
    design_formula
  )

  results_with_metadata <- attach_deseq2_metadata(results_processed, summary_stats)

  if (verbose) {
    cat(sprintf("Analysis complete: %d up, %d down, %d significant total\n",
      summary_stats$n_up, summary_stats$n_down, summary_stats$n_total_de))
  }

  list(
    results = results_with_metadata,
    summary = summary_stats,
    dds = dds_analyzed
  )
}
