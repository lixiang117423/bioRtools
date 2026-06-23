#' Identify Differentially Expressed Genes Using edgeR
#'
#' This function performs differential gene expression analysis using edgeR's
#' negative binomial models. It identifies genes that are significantly
#' up-regulated or down-regulated between experimental conditions, using TMM
#' normalization and empirical Bayes dispersion estimation. edgeR is
#' well-suited for experiments with small sample sizes and complex designs.
#'
#' @param data Integer count matrix where rows represent genes and columns
#'   represent samples. Row names should be gene identifiers (e.g., gene symbols,
#'   Ensembl IDs) and column names should match sample identifiers in the sample
#'   metadata. All values must be non-negative integers representing raw counts
#'   (not normalized, FPKM, or TPM values)
#' @param sample Data frame containing sample metadata where rows represent
#'   samples and columns contain experimental factors and covariates. Row names
#'   should exactly match column names in the count matrix. Must contain variables
#'   referenced in the design formula
#' @param formula Design formula specifying the experimental design for edgeR.
#'   Default is \code{~group}, assuming a "group" column exists in sample metadata.
#'   Common examples:
#'   \itemize{
#'     \item \code{~condition}: Simple two-group comparison
#'     \item \code{~condition + batch}: Control for batch effects
#'     \item \code{~condition + sex + age}: Multiple covariates
#'   }
#' @param log2_fold_change Minimum absolute log2 fold change threshold for
#'   significance classification (default: 1, equivalent to 2-fold change).
#'   Genes with |log2FC| >= this value AND fdr < fdr threshold are considered
#'   differentially expressed
#' @param fdr Adjusted p-value (FDR) threshold for statistical significance
#'   (default: 0.05). Uses Benjamini-Hochberg correction
#' @param normalize_method Normalization method for library size correction.
#'   One of "TMM" (default), "TMMwsp", "RLE", "upperquartile", "none".
#'   TMM is recommended for most RNA-seq experiments
#' @param filter_min_cpm Minimum counts per million (CPM) threshold for
#'   filtering low-expressed genes (default: 1). Genes must have CPM >= this
#'   value in at least \code{filter_min_samples} samples to be retained
#' @param filter_min_samples Minimum number of samples in which a gene must
#'   meet the CPM threshold (default: NULL, automatically set to the size of
#'   the smallest group). Controls the stringency of low-expression filtering
#' @param trend_method Method for fitting the mean-dispersion trend.
#'   One of "loess" (default), "locfit", "movingave", "locfit.mixed", "none".
#'   Default "loess" uses local regression. "locfit.mixed" is recommended for
#'   small datasets
#' @param robust Logical indicating whether to use robust estimation for
#'   the quasi-likelihood model (default: TRUE). Improves performance with
#'   outliers and unequal sample sizes
#' @param remove_problematic_genes Logical indicating whether to automatically
#'   remove genes with non-finite values (Inf, -Inf, NaN) (default: TRUE).
#'   If TRUE, problematic genes are removed with a warning. If FALSE, the function
#'   will stop with an error message
#' @param pairwise Logical indicating whether to perform pairwise comparisons
#'   between groups (default: TRUE). If TRUE, all pairs of groups are compared.
#'   If FALSE, uses the full design formula
#' @param ref_group Optional character string specifying the reference group
#'   for pairwise comparisons. If NULL (default), all pairwise comparisons
#'   are made. If specified, only comparisons against this reference group
#'   are performed
#'
#' @return A data frame containing differential expression results. In pairwise
#'   mode, each row represents a gene-comparison combination. Columns include:
#'   \itemize{
#'     \item \code{gene}: Gene identifiers
#'     \item \code{log2FoldChange}: Log2 fold change between conditions
#'     \item \code{logCPM}: Average log2 counts per million
#'     \item \code{pvalue}: Raw p-value from the statistical test
#'     \item \code{padj}: Benjamini-Hochberg adjusted p-value (FDR)
#'     \item \code{regulation}: Gene regulation classification:
#'       \itemize{
#'         \item "Up": log2FC > threshold & fdr < threshold
#'         \item "Down": log2FC < -threshold & fdr < threshold
#'         \item "NS": Does not meet both criteria
#'       }
#'     \item \code{fold_change}: Actual fold change (2^log2FoldChange)
#'     \item \code{abs_log2fc}: Absolute log2 fold change for ranking
#'     \item \code{comparison}: Comparison label (pairwise mode only)
#'     \item \code{group}: Test group name (pairwise mode only)
#'     \item \code{ref_group}: Reference group name (pairwise mode only)
#'     \item \code{group_mean}, \code{group_sd}, \code{group_n}:
#'       Descriptive stats for test group (pairwise mode)
#'     \item \code{ref_mean}, \code{ref_sd}, \code{ref_n}:
#'       Descriptive stats for reference group (pairwise mode)
#'     \item \code{test_method}: Statistical method used
#'     \item Original count columns from input data
#'   }
#'
#' @details
#' edgeR workflow implemented in this function:
#' \enumerate{
#'   \item \strong{Create DGEList}: Combines count matrix with group information
#'   \item \strong{Filter low-expression genes}: Removes genes below CPM threshold
#'   \item \strong{TMM normalization}: Corrects for library size and composition bias
#'   \item \strong{Dispersion estimation}: Estimates common, trended, and tagwise dispersions
#'   \item \strong{Model fitting}: Fits negative binomial GLM or uses exact test
#'   \item \strong{Statistical testing}: Tests for differential expression
#'   \item \strong{Multiple testing correction}: Applies Benjamini-Hochberg FDR
#' }
#'
#' Key differences from DESeq2:
#' \itemize{
#'   \item edgeR uses TMM normalization vs DESeq2's median-of-ratios
#'   \item edgeR's quasi-likelihood approach provides more conservative results
#'   \item edgeR is generally faster for large datasets
#'   \item edgeR handles small sample sizes well with empirical Bayes shrinkage
#' }
#'
#' @note
#' \itemize{
#'   \item Input must be raw counts (integers), not normalized data
#'   \item Low-expressed genes are filtered by CPM before analysis
#'   \item edgeR handles library size normalization internally via TMM
#'   \item For simple two-group comparisons, exact test is used; for complex
#'     designs, GLM quasi-likelihood F-test is used
#'   \item Results are sorted by adjusted p-value and absolute fold change
#' }
#'
#' @references
#' Robinson MD, McCarthy DJ, Smyth GK (2010). edgeR: a Bioconductor package
#' for differential expression analysis of digital gene expression data.
#' Bioinformatics, 26(1), 139-140.
#'
#' McCarthy DJ, Chen Y, Smyth GK (2012). Differential expression analysis of
#' multifactor RNA-Seq experiments with respect to biological variation.
#' Nucleic Acids Research, 40(10), 4288-4297.
#'
#' Chen Y, Lun ATL, Smyth GK (2016). From reads to genes to pathways:
#' differential expression analysis of RNA-Seq experiments using Rsubread and
#' the edgeR quasi-likelihood pipeline. F1000Research, 5, 1438.
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{find_degs_deseq2}} for DESeq2-based differential expression
#'   \item \code{\link{volcano_plot}} for visualization of DE results
#'   \item \code{\link{enrich_go}} for functional enrichment of DE genes
#' }
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' library(bioRtools)
#'
#' # Example 1: Basic usage with simulated data
#' set.seed(123)
#' n_genes <- 1000
#' n_samples <- 12
#'
#' mock_counts <- matrix(
#'   rnbinom(n_genes * n_samples, size = 10, mu = 50),
#'   nrow = n_genes,
#'   ncol = n_samples,
#'   dimnames = list(
#'     paste0("Gene_", 1:n_genes),
#'     paste0("Sample_", 1:n_samples)
#'   )
#' )
#'
#' mock_sample <- data.frame(
#'   row.names = colnames(mock_counts),
#'   group = rep(c("Control", "Treatment"), each = 6)
#' )
#'
#' \dontrun{
#' de_results <- find_degs_edger(
#'   data = mock_counts,
#'   sample = mock_sample
#' )
#'
#' print(table(de_results$regulation))
#' }
#'
#' # Example 2: Using DESeq2-style data
#' \dontrun{
#' data(df.rnaseq.gene)
#' data(df.rnaseq.sample)
#'
#' de_results <- find_degs_edger(
#'   data = df.rnaseq.gene,
#'   sample = df.rnaseq.sample,
#'   log2_fold_change = 1,
#'   fdr = 0.05
#' )
#'
#' print(head(de_results))
#' }
#'
#' # Example 3: With reference group and custom thresholds
#' \dontrun{
#' de_ref <- find_degs_edger(
#'   data = df.rnaseq.gene,
#'   sample = df.rnaseq.sample,
#'   ref_group = "Control",
#'   log2_fold_change = 1.5,
#'   fdr = 0.01
#' )
#' }
find_degs_edger <- function(data, sample, formula = ~group, log2_fold_change = 1, fdr = 0.05,
                            normalize_method = "TMM", filter_min_cpm = 1,
                            filter_min_samples = NULL, trend_method = "loess",
                            robust = TRUE, remove_problematic_genes = TRUE,
                            pairwise = TRUE, ref_group = NULL) {

  # --- Input validation -------------------------------------------------------

  if (!is.matrix(data) && !is.data.frame(data)) {
    stop(err_invalid_input("data", "a matrix or data frame"))
  }

  if (!is.data.frame(sample)) {
    stop(err_invalid_input("sample", "a data frame"))
  }

  if (nrow(data) == 0 || ncol(data) == 0) {
    stop(err_missing_required("data"))
  }

  if (nrow(sample) == 0) {
    stop(err_missing_required("sample"))
  }

  if (!is.logical(remove_problematic_genes) || length(remove_problematic_genes) != 1) {
    stop(err_invalid_input("remove_problematic_genes", "a single logical value (TRUE or FALSE)"))
  }

  if (any(data < 0, na.rm = TRUE)) {
    min_val <- min(data, na.rm = TRUE)
    stop(err_negative_values("data", min_val))
  }

  # --- Prepare data matrix ----------------------------------------------------

  data <- as.data.frame(data)
  numeric_cols <- sapply(data, is.numeric)
  if (!all(numeric_cols)) {
    if (sum(!numeric_cols) == 1) {
      rownames(data) <- data[[which(!numeric_cols)]]
      data <- data[, numeric_cols, drop = FALSE]
    } else {
      data <- data[, numeric_cols, drop = FALSE]
    }
  }
  data_matrix <- as.matrix(data)
  n_genes_original <- nrow(data_matrix)

  # Save original expression values for output
  expr_df <- as.data.frame(data_matrix)
  expr_df$gene <- rownames(expr_df)
  rownames(expr_df) <- NULL

  # Round non-integer values
  if (!is.integer(data_matrix)) {
    if (any(data_matrix != round(data_matrix), na.rm = TRUE)) {
      warning("Non-integer values detected. Rounding to integers. edgeR requires raw counts; results from rounded normalized values may be less reliable.")
    }
    data_matrix[!is.finite(data_matrix)] <- 0
    data_matrix <- round(data_matrix)
  }

  has_nonfinite <- any(!is.finite(data_matrix), na.rm = TRUE)

  if (has_nonfinite) {
    problem_genes_idx <- which(apply(data_matrix, 1, function(x) any(!is.finite(x))))
    problem_samples_idx <- which(apply(data_matrix, 2, function(x) any(!is.finite(x))))

    gene_names <- if (!is.null(rownames(data_matrix))) rownames(data_matrix) else paste0("Row_", 1:nrow(data_matrix))
    sample_names_v <- if (!is.null(colnames(data_matrix))) colnames(data_matrix) else paste0("Col_", 1:ncol(data_matrix))

    problematic_gene_names <- gene_names[problem_genes_idx]
    problematic_sample_names <- sample_names_v[problem_samples_idx]

    if (remove_problematic_genes) {
      warn_genes_removed(length(problematic_gene_names), "non-finite values (Inf, -Inf, NaN)")
      data_matrix <- data_matrix[-problem_genes_idx, , drop = FALSE]
      if (nrow(data_matrix) == 0) {
        stop("All genes removed after filtering. Check your data - edgeR requires raw integer counts, not normalized values (FPKM/TPM/RPKM).")
      }
    } else {
      stop(err_infinite_values(problem_genes_idx, problem_samples_idx))
    }
  }

  storage.mode(data_matrix) <- "integer"

  # --- Validate sample-column correspondence ----------------------------------

  if (ncol(data_matrix) != nrow(sample)) {
    stop("Number of columns in 'data' must equal number of rows in 'sample'")
  }

  data_samples <- colnames(data_matrix)
  sample_names <- rownames(sample)

  if (is.null(sample_names) || identical(sample_names, as.character(seq_len(nrow(sample))))) {
    for (col in names(sample)) {
      if (is.character(sample[[col]]) && all(data_samples %in% sample[[col]])) {
        rownames(sample) <- sample[[col]]
        sample_names <- sample[[col]]
        break
      }
    }
  }

  if (is.null(data_samples) || is.null(sample_names)) {
    stop(err_row_col_mismatch("data", "sample"))
  }

  if (!all(data_samples %in% sample_names)) {
    stop(err_row_col_mismatch("data", "sample metadata"))
  }

  if (!all(sample_names %in% data_samples)) {
    extra_samples <- setdiff(sample_names, data_samples)
    warn_data_processing(paste("Extra samples in metadata (ignored):", paste(extra_samples, collapse = ", ")))
  }

  sample_aligned <- sample[data_samples, , drop = FALSE]

  # --- Validate formula -------------------------------------------------------

  if (!inherits(formula, "formula")) {
    stop(err_invalid_input("formula", "a formula object (e.g., ~condition)"))
  }

  validate_formula_vars(formula, sample_aligned)

  formula_vars <- all.vars(formula)
  if (length(formula_vars) == 0) {
    stop(err_invalid_input("formula", "at least one variable"))
  }

  for (var in formula_vars) {
    if (any(is.na(sample_aligned[[var]]))) {
      stop(err_missing_values(sprintf("design variable '%s'", var), sum(is.na(sample_aligned[[var]]))))
    }
  }

  # --- Validate parameters ----------------------------------------------------

  if (!is.numeric(log2_fold_change) || length(log2_fold_change) != 1 || log2_fold_change < 0) {
    stop("'log2_fold_change' must be a single non-negative number")
  }

  if (!is.numeric(fdr) || length(fdr) != 1 || fdr <= 0 || fdr > 1) {
    stop("'fdr' must be a single number between 0 and 1")
  }

  valid_norm_methods <- c("TMM", "TMMwsp", "RLE", "upperquartile", "none")
  if (!normalize_method %in% valid_norm_methods) {
    stop(err_invalid_input("normalize_method", paste("one of:", paste(valid_norm_methods, collapse = ", "))))
  }

  if (!is.numeric(filter_min_cpm) || length(filter_min_cpm) != 1 || filter_min_cpm < 0) {
    stop("'filter_min_cpm' must be a single non-negative number")
  }

  valid_trend_methods <- c("loess", "locfit", "movingave", "locfit.mixed", "none")
  if (!trend_method %in% valid_trend_methods) {
    stop(err_invalid_input("trend_method", paste("one of:", paste(valid_trend_methods, collapse = ", "))))
  }

  if (!is.logical(robust) || length(robust) != 1) {
    stop("'robust' must be a single logical value")
  }

  # --- Check experimental design ----------------------------------------------

  main_factor <- formula_vars[length(formula_vars)]
  if (main_factor %in% names(sample_aligned)) {
    group_counts <- table(sample_aligned[[main_factor]])
    n_groups <- length(group_counts)

    if (n_groups < 2) {
      stop("Need at least 2 groups for differential expression analysis")
    }

    if (any(group_counts < 2)) {
      stop("Each group must have at least 2 biological replicates")
    }

    if (any(group_counts < 3)) {
      warning("Groups with fewer than 3 replicates may give unreliable results. 6+ replicates recommended.")
    }

    min_replicates <- min(group_counts)
    if (min_replicates < 6) {
      warning(paste("Minimum replicates per group:", min_replicates,
        ". Consider 6+ replicates for robust differential expression."))
    }
  }

  # Check library sizes
  lib_sizes <- colSums(data_matrix)
  lib_size_ratio <- max(lib_sizes) / min(lib_sizes)

  if (lib_size_ratio > 10) {
    warning(paste("Large library size differences detected (max/min ratio:",
      round(lib_size_ratio, 1), "). Consider checking data quality."))
  }

  if (min(lib_sizes) < 1000) {
    warning("Very small library sizes detected (<1000 counts). Results may be unreliable.")
  }

  # --- Set filter_min_samples -------------------------------------------------

  if (is.null(filter_min_samples)) {
    filter_min_samples <- min(group_counts)
  }

  # --- Pairwise mode ----------------------------------------------------------

  if (pairwise) {
    groups <- unique(as.character(sample_aligned[[main_factor]]))
    if (length(groups) < 2) stop("Need at least 2 groups for pairwise analysis")

    if (!is.null(ref_group)) {
      if (!ref_group %in% groups) stop("'ref_group' must be one of: ", paste(groups, collapse = ", "))
      others <- setdiff(groups, ref_group)
      pairs <- lapply(others, function(g) c(ref_group, g))
    } else {
      pairs <- utils::combn(groups, 2, simplify = FALSE)
    }

    run_one_pair <- function(pair) {
      idx <- sample_aligned[[main_factor]] %in% pair
      d_sub <- data_matrix[, idx, drop = FALSE]
      s_sub <- sample_aligned[idx, , drop = FALSE]
      s_sub[[main_factor]] <- factor(as.character(s_sub[[main_factor]]), levels = pair)

      if (any(table(s_sub[[main_factor]]) < 2)) return(NULL)

      comp_label <- paste(pair[2], "vs", pair[1])
      group_factor <- s_sub[[main_factor]]

      # Create DGEList
      dge <- tryCatch(
        edgeR::DGEList(counts = d_sub, group = group_factor),
        error = function(e) NULL
      )
      if (is.null(dge)) return(NULL)

      # Filter low-expression genes
      keep <- edgeR::filterByExpr(dge, group = group_factor, min.count = filter_min_cpm)
      if (sum(keep) == 0) return(NULL)
      dge <- dge[keep, , keep.lib.sizes = FALSE]

      # Normalize
      dge <- tryCatch(
        edgeR::calcNormFactors(dge, method = normalize_method),
        error = function(e) {
          warning(paste("Normalization failed for", comp_label, ":", e$message))
          dge
        }
      )

      # Estimate dispersion
      dge <- tryCatch(
        edgeR::estimateDisp(dge, design = model.matrix(~group_factor), robust = robust),
        error = function(e) NULL
      )
      if (is.null(dge)) return(NULL)

      # Use exact test for two-group comparison
      res_raw <- tryCatch(
        edgeR::exactTest(dge, pair = as.character(pair)),
        error = function(e) NULL
      )
      if (is.null(res_raw)) return(NULL)

      res_table <- edgeR::topTags(res_raw, n = Inf, sort.by = "none")
      res <- as.data.frame(res_table)

      res <- tibble::rownames_to_column(res, "gene")
      colnames(res)[colnames(res) == "logFC"] <- "log2FoldChange"
      colnames(res)[colnames(res) == "PValue"] <- "pvalue"
      colnames(res)[colnames(res) == "FDR"] <- "padj"

      res$padj <- ifelse(is.na(res$padj), 1, res$padj)
      res$regulation <- ifelse(res$log2FoldChange > log2_fold_change & res$padj < fdr, "Up",
                        ifelse(res$log2FoldChange < -log2_fold_change & res$padj < fdr, "Down", "NS"))
      res$fold_change <- 2^res$log2FoldChange
      res$abs_log2fc <- abs(res$log2FoldChange)
      res$comparison <- comp_label
      res$group <- pair[2]
      res$ref_group <- pair[1]

      # Group descriptive statistics (CPM-normalized)
      cpm_vals <- edgeR::cpm(dge, normalized.lib.sizes = TRUE)
      idx_grp <- group_factor == pair[2]
      idx_ref <- group_factor == pair[1]
      n_grp <- sum(idx_grp)
      n_ref <- sum(idx_ref)

      grp_means <- rowMeans(cpm_vals[idx_grp, , drop = FALSE], na.rm = TRUE)
      grp_sds <- apply(cpm_vals[idx_grp, , drop = FALSE], 1, sd, na.rm = TRUE)
      ref_means <- rowMeans(cpm_vals[idx_ref, , drop = FALSE], na.rm = TRUE)
      ref_sds <- apply(cpm_vals[idx_ref, , drop = FALSE], 1, sd, na.rm = TRUE)

      # Match gene names back
      gene_idx <- match(res$gene, rownames(cpm_vals))
      res$group_mean <- round(grp_means[gene_idx], 2)
      res$group_sd <- round(grp_sds[gene_idx], 2)
      res$group_n <- n_grp
      res$ref_mean <- round(ref_means[gene_idx], 2)
      res$ref_sd <- round(ref_sds[gene_idx], 2)
      res$ref_n <- n_ref
      res$test_method <- "edgeR-exactTest"

      res
    }

    pair_results <- lapply(pairs, run_one_pair)
    pair_results <- pair_results[!sapply(pair_results, is.null)]

    if (length(pair_results) == 0) {
      message("No pairwise comparisons produced results")
      return(data.frame())
    }

    combined <- do.call(rbind, pair_results)
    rownames(combined) <- NULL

    combined <- dplyr::left_join(combined, expr_df, by = "gene")
    combined <- combined[order(combined$comparison, combined$padj, -combined$abs_log2fc), ]

    if (interactive()) {
      cat("edgeR Pairwise Differential Expression Summary:\n")
      cat("==================================================\n")
      cat("Comparisons:", length(pair_results), "\n")
      for (comp in unique(combined$comparison)) {
        sub <- combined[combined$comparison == comp, ]
        n_up <- sum(sub$regulation == "Up")
        n_down <- sum(sub$regulation == "Down")
        cat(sprintf("  %s: %d up, %d down, %d total DE\n", comp, n_up, n_down, n_up + n_down))
      }
      cat("\n")
    }

    return(combined)
  }

  # --- Single comparison (GLM approach) ---------------------------------------

  group_factor <- sample_aligned[[main_factor]]
  design_mat <- model.matrix(formula, data = sample_aligned)

  # Create DGEList
  dge <- tryCatch(
    edgeR::DGEList(counts = data_matrix, group = group_factor),
    error = function(e) {
      stop(paste("Error creating DGEList:", e$message))
    }
  )

  # Filter low-expression genes
  keep <- edgeR::filterByExpr(dge, design = design_mat, min.count = filter_min_cpm)
  n_filtered <- sum(!keep)
  if (n_filtered > 0) {
    message(paste("Filtered", n_filtered, "low-expression genes (CPM <", filter_min_cpm,
      "in fewer than", filter_min_samples, "samples)"))
  }

  if (sum(keep) == 0) {
    stop("All genes filtered out. Consider lowering 'filter_min_cpm' or checking data quality.")
  }

  dge <- dge[keep, , keep.lib.sizes = FALSE]

  # Normalize
  dge <- tryCatch(
    edgeR::calcNormFactors(dge, method = normalize_method),
    error = function(e) {
      warning(paste("Normalization failed:", e$message, "Proceeding without normalization factors."))
      dge
    }
  )

  # Estimate dispersion
  dge <- tryCatch(
    edgeR::estimateDisp(dge, design = design_mat, robust = robust, trend.method = trend_method),
    error = function(e) {
      stop(paste("Error estimating dispersions:", e$message))
    }
  )

  # Fit GLM
  fit <- tryCatch(
    edgeR::glmQLFit(dge, design = design_mat, robust = robust),
    error = function(e) {
      stop(paste("Error fitting GLM:", e$message))
    }
  )

  # Test for DE genes (last coefficient = main effect)
  coef_to_test <- ncol(design_mat)
  res_raw <- tryCatch(
    edgeR::glmQLFTest(fit, coef = coef_to_test),
    error = function(e) {
      # Fallback: try likelihood ratio test
      fit_lr <- edgeR::glmFit(dge, design = design_mat)
      edgeR::glmLRT(fit_lr, coef = coef_to_test)
    }
  )

  res_table <- edgeR::topTags(res_raw, n = Inf, sort.by = "none")
  results_df <- as.data.frame(res_table)

  # Rename columns to match DESeq2 output
  colnames(results_df)[colnames(results_df) == "logFC"] <- "log2FoldChange"
  colnames(results_df)[colnames(results_df) == "PValue"] <- "pvalue"
  colnames(results_df)[colnames(results_df) == "FDR"] <- "padj"

  # Process and enhance results
  results_processed <- results_df %>%
    tibble::rownames_to_column(var = "gene") %>%
    dplyr::mutate(
      regulation = dplyr::case_when(
        log2FoldChange > !!log2_fold_change & padj < !!fdr ~ "Up-regulated",
        log2FoldChange < -!!log2_fold_change & padj < !!fdr ~ "Down-regulated",
        TRUE ~ "Not significant"
      ),
      fold_change = 2^log2FoldChange,
      abs_log2fc = abs(log2FoldChange),
      padj = ifelse(is.na(padj), 1, padj),
      significance_level = dplyr::case_when(
        padj >= 0.05 ~ "NS (p≥0.05)",
        padj >= 0.01 ~ "* (p<0.05)",
        padj >= 0.001 ~ "** (p<0.01)",
        padj >= 0.0001 ~ "*** (p<0.001)",
        TRUE ~ "**** (p<0.0001)"
      )
    ) %>%
    dplyr::left_join(expr_df, by = "gene") %>%
    dplyr::arrange(padj, desc(abs_log2fc))

  # Calculate summary statistics
  n_upregulated <- sum(results_processed$regulation == "Up-regulated", na.rm = TRUE)
  n_downregulated <- sum(results_processed$regulation == "Down-regulated", na.rm = TRUE)
  n_total_de <- n_upregulated + n_downregulated
  n_tested <- sum(!is.na(results_processed$padj))

  # Add metadata
  attr(results_processed, "n_genes_input") <- nrow(data_matrix)
  attr(results_processed, "n_genes_original") <- n_genes_original
  attr(results_processed, "n_genes_tested") <- n_tested
  attr(results_processed, "n_samples") <- ncol(data_matrix)
  attr(results_processed, "design_formula") <- deparse(formula)
  attr(results_processed, "log2fc_threshold") <- log2_fold_change
  attr(results_processed, "fdr_threshold") <- fdr
  attr(results_processed, "normalize_method") <- normalize_method
  attr(results_processed, "robust") <- robust
  attr(results_processed, "filter_min_cpm") <- filter_min_cpm
  attr(results_processed, "library_size_range") <- range(lib_sizes)
  attr(results_processed, "remove_problematic_genes") <- remove_problematic_genes
  if (has_nonfinite && remove_problematic_genes) {
    attr(results_processed, "n_genes_removed") <- length(problematic_gene_names)
    attr(results_processed, "genes_removed") <- problematic_gene_names
  }

  # Interactive summary
  if (interactive()) {
    cat("edgeR Differential Expression Analysis Summary:\n")
    cat("================================================\n")

    if (has_nonfinite && remove_problematic_genes) {
      cat("Genes (original):", n_genes_original, "\n")
      cat("Genes removed (non-finite values):", length(problematic_gene_names), "\n")
      cat("Genes after CPM filter:", nrow(dge), "\n")
    } else {
      cat("Genes input:", n_genes_original, "\n")
      cat("Genes after CPM filter:", nrow(dge), "\n")
    }
    cat("Genes tested:", n_tested, "\n")
    cat("Samples analyzed:", ncol(data_matrix), "\n")
    cat("Design formula:", deparse(formula), "\n")
    cat("Normalization:", normalize_method, "\n")
    cat("Log2FC threshold:", log2_fold_change, "(", round(2^log2_fold_change, 2), "-fold)\n")
    cat("FDR threshold:", fdr, "\n")
    cat("Robust estimation:", ifelse(robust, "Yes", "No"), "\n\n")

    cat("Library size summary:\n")
    cat("  Range:", min(lib_sizes), "to", max(lib_sizes), "counts\n")
    cat("  Median:", median(lib_sizes), "counts\n")

    if (main_factor %in% names(sample_aligned)) {
      cat("  Group sizes:", paste(names(group_counts), "=", group_counts, collapse = ", "), "\n\n")
    }

    cat("Results:\n")
    cat("  Up-regulated genes:", n_upregulated, "\n")
    cat("  Down-regulated genes:", n_downregulated, "\n")
    cat("  Not significant:", sum(results_processed$regulation == "Not significant"), "\n")
    cat("  Total DE genes:", n_total_de,
      paste0("(", round(100 * n_total_de / max(n_tested, 1), 1), "% of tested)"), "\n\n")

    if (n_total_de > 0) {
      cat("Top DE genes by significance:\n")
      top_genes <- results_processed %>%
        dplyr::filter(regulation != "Not significant") %>%
        head(5)

      for (i in seq_len(nrow(top_genes))) {
        direction_symbol <- ifelse(top_genes$log2FoldChange[i] > 0, "↑", "↓")
        cat(sprintf("  %s %s: %s (log2FC = %.2f, padj = %.2e)\n",
          direction_symbol,
          top_genes$gene[i],
          top_genes$regulation[i],
          top_genes$log2FoldChange[i],
          top_genes$padj[i]))
      }

      cat("\nExpression change magnitude distribution:\n")
      magnitude_summary <- results_processed %>%
        dplyr::filter(regulation != "Not significant") %>%
        dplyr::mutate(
          magnitude = dplyr::case_when(
            abs_log2fc >= 3 ~ ">8-fold",
            abs_log2fc >= 2 ~ "4-8 fold",
            abs_log2fc >= 1 ~ "2-4 fold",
            TRUE ~ "<2-fold"
          )
        ) %>%
        dplyr::count(regulation, magnitude) %>%
        tidyr::pivot_wider(names_from = regulation, values_from = n, values_fill = 0)

      print(magnitude_summary)
    }
    cat("\n")
  }

  results_processed
}
