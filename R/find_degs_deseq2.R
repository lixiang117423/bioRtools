#' Identify Differentially Expressed Genes Using DESeq2
#'
#' This function performs differential gene expression analysis using DESeq2's
#' negative binomial generalized linear model. It identifies genes that are
#' significantly up-regulated or down-regulated between experimental conditions,
#' accounting for library size differences and biological variability. DESeq2
#' is the gold standard for RNA-seq differential expression analysis.
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
#' @param formula Design formula specifying the experimental design for DESeq2.
#'   Default is \code{~group}, assuming a "group" column exists in sample metadata.
#'   Common examples:
#'   \itemize{
#'     \item \code{~condition}: Simple two-group comparison
#'     \item \code{~condition + batch}: Control for batch effects
#'     \item \code{~condition + sex + age}: Multiple covariates
#'     \item \code{~condition * time}: Interaction terms for time series
#'   }
#' @param log2FoldChange Minimum absolute log2 fold change threshold for
#'   significance classification (default: 1, equivalent to 2-fold change).
#'   Genes with |log2FC| >= this value AND padj < padj threshold are considered
#'   differentially expressed. Use 0.5 for 1.4-fold, 1.5 for 3-fold changes
#' @param padj Adjusted p-value (FDR) threshold for statistical significance
#'   (default: 0.05). Genes with Benjamini-Hochberg adjusted p-values below
#'   this threshold are considered statistically significant
#' @param shrink.lfc Logical indicating whether to apply log2 fold change
#'   shrinkage (default: TRUE). Shrinkage reduces noise in log2FC estimates
#'   for genes with low counts, providing more accurate effect size estimates
#' @param independent.filtering Logical indicating whether to perform independent
#'   filtering (default: TRUE). Filters out genes with very low counts to
#'   improve multiple testing correction and statistical power
#' @param alpha Alpha level for outlier detection and Cook's distance filtering
#'   (default: 0.1). Lower values are more stringent for outlier detection
#' @param remove_problematic_genes Logical indicating whether to automatically
#'   remove genes with non-finite values (Inf, -Inf, NaN) (default: TRUE).
#'   If TRUE, problematic genes are removed with a warning. If FALSE, the function
#'   will stop with an error message detailing which genes and samples have issues
#'
#' @return A data frame containing differential expression results with columns:
#'   \itemize{
#'     \item \code{gene}: Gene identifiers (from row names of input matrix)
#'     \item \code{baseMean}: Mean normalized counts across all samples
#'     \item \code{log2FoldChange}: Log2 fold change between conditions
#'     \item \code{lfcSE}: Standard error of log2 fold change estimate
#'     \item \code{stat}: Wald test statistic
#'     \item \code{pvalue}: Raw p-value from Wald test
#'     \item \code{padj}: Benjamini-Hochberg adjusted p-value (FDR)
#'     \item \code{regulation}: Gene regulation classification:
#'       \itemize{
#'         \item "Up-regulated": log2FC > threshold & padj < threshold
#'         \item "Down-regulated": log2FC < -threshold & padj < threshold
#'         \item "Not significant": Does not meet both criteria
#'       }
#'     \item \code{fold_change}: Actual fold change (2^log2FoldChange)
#'     \item \code{abs_log2fc}: Absolute log2 fold change for ranking
#'     \item \code{significance_level}: Detailed significance classification
#'   }
#'
#' @details
#' DESeq2 workflow implemented in this function:
#' \enumerate{
#'   \item \strong{Create DESeqDataSet}: Combines count matrix with sample metadata
#'   \item \strong{Size factor estimation}: Normalizes for library size differences
#'   \item \strong{Dispersion estimation}: Estimates gene-wise and fitted dispersions
#'   \item \strong{Negative binomial GLM fitting}: Fits model for each gene
#'   \item \strong{Wald tests}: Tests for differential expression
#'   \item \strong{Multiple testing correction}: Applies Benjamini-Hochberg FDR
#'   \item \strong{Log2FC shrinkage}: Reduces noise in effect size estimates (optional)
#' }
#'
#' Key assumptions and requirements:
#' \itemize{
#'   \item Count data follows negative binomial distribution
#'   \item Genes are independent (after accounting for design factors)
#'   \item Samples are independent biological replicates
#'   \item At least 3 replicates per condition (6+ recommended)
#'   \item Input data are raw integer counts, not normalized values
#' }
#'
#' @note
#' \itemize{
#'   \item Input must be raw counts (integers), not normalized data
#'   \item Genes with very low counts may be filtered automatically
#'   \item DESeq2 handles library size normalization internally
#'   \item Consider batch correction for multi-batch experiments
#'   \item Results are sorted by adjusted p-value and absolute fold change
#' }
#'
#' @references
#' Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold change
#' and dispersion for RNA-seq data with DESeq2. Genome Biology, 15, 550.
#'
#' Anders, S. and Huber, W. (2010) Differential expression analysis for sequence
#' count data. Genome Biology, 11, R106.
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{volcano_plot}} for visualization of DE results
#'   \item \code{\link{enrich_go}} for functional enrichment of DE genes
#'   \item \code{\link{enrich_kegg}} for pathway enrichment analysis
#' }
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' library(bioRtools)
#' library(dplyr)
#'
#' # Example 1: Basic differential expression analysis
#' \dontrun{
#' # Load example RNA-seq data
#' data(df.rnaseq.gene)    # Gene count matrix
#' data(df.rnaseq.sample)  # Sample metadata
#'
#' # Basic DE analysis
#' de_results <- find_degs_deseq2(
#'   data = df.rnaseq.gene,
#'   sample = df.rnaseq.sample
#' )
#'
#' # View results summary
#' print(table(de_results$regulation))
#'
#' # Top up-regulated genes
#' upregulated <- de_results %>%
#'   filter(regulation == "Up-regulated") %>%
#'   arrange(padj) %>%
#'   head(10)
#'
#' print("Top up-regulated genes:")
#' print(upregulated[, c("gene", "log2FoldChange", "padj", "baseMean")])
#'
#' # Top down-regulated genes
#' downregulated <- de_results %>%
#'   filter(regulation == "Down-regulated") %>%
#'   arrange(padj) %>%
#'   head(10)
#'
#' print("Top down-regulated genes:")
#' print(downregulated[, c("gene", "log2FoldChange", "padj", "baseMean")])
#' }
#'
#' # Example 2: Custom thresholds and batch correction
#' \dontrun{
#' # More stringent analysis with batch correction
#' de_strict <- find_degs_deseq2(
#'   data = df.rnaseq.gene,
#'   sample = df.rnaseq.sample,
#'   formula = ~ batch + condition,    # Control for batch effects
#'   log2FoldChange = 1.5,           # 3-fold change threshold
#'   padj = 0.01,                    # 1% FDR
#'   shrink.lfc = TRUE,              # Apply LFC shrinkage
#'   independent.filtering = TRUE     # Enable independent filtering
#' )
#'
#' print(paste("Strict analysis found",
#'   sum(de_strict$regulation != "Not significant"),
#'   "DE genes"))
#'
#' # Compare effect sizes
#' effect_comparison <- de_strict %>%
#'   filter(regulation != "Not significant") %>%
#'   summarise(
#'     mean_abs_lfc = mean(abs_log2fc),
#'     median_abs_lfc = median(abs_log2fc),
#'     max_abs_lfc = max(abs_log2fc)
#'   )
#' print("Effect size summary:")
#' print(effect_comparison)
#' }
#'
#' # Example 3: Time series analysis
#' \dontrun{
#' # Longitudinal RNA-seq experiment
#' # Assume sample metadata has 'timepoint' and 'subject' columns
#' timeseries_de <- find_degs_deseq2(
#'   data = df.rnaseq.gene,
#'   sample = df.rnaseq.sample,
#'   formula = ~ subject + timepoint,  # Control for subject effects
#'   log2FoldChange = 0.5,           # More sensitive for time effects
#'   padj = 0.1                      # Less stringent for discovery
#' )
#'
#' # Identify time-responsive genes
#' time_responsive <- timeseries_de %>%
#'   filter(regulation != "Not significant") %>%
#'   arrange(desc(abs_log2fc))
#'
#' print(paste("Time-responsive genes:", nrow(time_responsive)))
#'
#' # Expression magnitude categories
#' magnitude_summary <- time_responsive %>%
#'   mutate(
#'     magnitude = case_when(
#'       abs_log2fc >= 2 ~ "High (4+ fold)",
#'       abs_log2fc >= 1 ~ "Medium (2-4 fold)",
#'       TRUE ~ "Low (<2 fold)"
#'     )
#'   ) %>%
#'   count(regulation, magnitude)
#'
#' print("Expression change magnitudes:")
#' print(magnitude_summary)
#' }
#'
#' # Example 4: Creating sample count data
#' # This demonstrates the required data format
#' set.seed(123)
#' n_genes <- 1000
#' n_samples <- 12
#'
#' # Simulate count data (negative binomial distribution)
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
#' # Add some differentially expressed genes
#' de_genes <- sample(1:n_genes, 100)  # 10% DE genes
#' up_genes <- de_genes[1:50]
#' down_genes <- de_genes[51:100]
#'
#' # Simulate higher expression in treatment samples (7-12)
#' mock_counts[up_genes, 7:12] <- mock_counts[up_genes, 7:12] * 3
#' mock_counts[down_genes, 7:12] <- round(mock_counts[down_genes, 7:12] * 0.3)
#'
#' # Create sample metadata
#' mock_sample <- data.frame(
#'   row.names = colnames(mock_counts),
#'   group = rep(c("Control", "Treatment"), each = 6),
#'   batch = rep(c("Batch1", "Batch2"), times = 6),
#'   replicate = rep(1:6, times = 2)
#' )
#'
#' \dontrun{
#' # Run DE analysis on simulated data
#' mock_de <- find_degs_deseq2(
#'   data = mock_counts,
#'   sample = mock_sample,
#'   formula = ~group
#' )
#'
#' print("Simulated data DE analysis:")
#' print(table(mock_de$regulation))
#'
#' # Check recovery of simulated DE genes
#' recovered_up <- sum(paste0("Gene_", up_genes) %in%
#'   mock_de$gene[mock_de$regulation == "Up-regulated"])
#' recovered_down <- sum(paste0("Gene_", down_genes) %in%
#'   mock_de$gene[mock_de$regulation == "Down-regulated"])
#'
#' print(paste("Recovered", recovered_up, "of", length(up_genes), "up-regulated genes"))
#' print(paste("Recovered", recovered_down, "of", length(down_genes), "down-regulated genes"))
#' }
#'
#' # Example 5: Quality control and filtering
#' \dontrun{
#' # Pre-analysis quality control
#' raw_counts <- df.rnaseq.gene
#'
#' # Check sequencing depth
#' library_sizes <- colSums(raw_counts)
#' print(paste("Library sizes range:", min(library_sizes), "to", max(library_sizes)))
#'
#' # Filter low-expressed genes (optional pre-filtering)
#' # Keep genes with at least 10 counts in at least 25% of samples
#' min_counts <- 10
#' min_samples <- ceiling(0.25 * ncol(raw_counts))
#' keep_genes <- rowSums(raw_counts >= min_counts) >= min_samples
#'
#' filtered_counts <- raw_counts[keep_genes, ]
#' print(paste("Retained", nrow(filtered_counts), "of", nrow(raw_counts),
#'   "genes after filtering"))
#'
#' # Run analysis on filtered data
#' filtered_de <- find_degs_deseq2(
#'   data = filtered_counts,
#'   sample = df.rnaseq.sample,
#'   formula = ~group
#' )
#'
#' print("Filtered analysis results:")
#' print(table(filtered_de$regulation))
#' }
#'
find_degs_deseq2 <- function(data, sample, formula = ~group, log2FoldChange = 1, padj = 0.05,
                             shrink.lfc = TRUE, independent.filtering = TRUE, alpha = 0.1,
                             remove_problematic_genes = TRUE, pairwise = TRUE,
                             ref_group = NULL) {
  # Input validation
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

  # Validate remove_problematic_genes parameter
  if (!is.logical(remove_problematic_genes) || length(remove_problematic_genes) != 1) {
    stop(err_invalid_input("remove_problematic_genes", "a single logical value (TRUE or FALSE)"))
  }

  # Check if data contains valid count values
  if (any(data < 0, na.rm = TRUE)) {
    min_val <- min(data, na.rm = TRUE)
    stop(err_negative_values("data", min_val))
  }

  # Prepare data matrix: detect and remove non-numeric columns (e.g., gene_id)
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

  # Save original expression values for output (before rounding/filtering)
  expr_df <- as.data.frame(data_matrix)
  expr_df$gene <- rownames(expr_df)
  rownames(expr_df) <- NULL

  # Round non-integer values first (FPKM/TPM → integer for DESeq2)
  if (!is.integer(data_matrix)) {
    if (any(data_matrix != round(data_matrix), na.rm = TRUE)) {
      warning("Non-integer values detected. Rounding to integers. DESeq2 requires raw counts; results from rounded normalized values may be less reliable.")
    }
    data_matrix[!is.finite(data_matrix)] <- 0
    data_matrix <- round(data_matrix)
  }

  has_nonfinite <- any(!is.finite(data_matrix), na.rm = TRUE)

  if (has_nonfinite) {
    # Identify problematic genes and samples
    problem_genes_idx <- which(apply(data_matrix, 1, function(x) any(!is.finite(x))))
    problem_samples_idx <- which(apply(data_matrix, 2, function(x) any(!is.finite(x))))

    gene_names <- if (!is.null(rownames(data_matrix))) rownames(data_matrix) else paste0("Row_", 1:nrow(data_matrix))
    sample_names <- if (!is.null(colnames(data_matrix))) colnames(data_matrix) else paste0("Col_", 1:ncol(data_matrix))

    problematic_gene_names <- gene_names[problem_genes_idx]
    problematic_sample_names <- sample_names[problem_samples_idx]

    if (remove_problematic_genes) {
      # Auto-remove problematic genes
      warn_genes_removed(length(problematic_gene_names), "non-finite values (Inf, -Inf, NaN)")

      # Remove problematic genes
      data_matrix <- data_matrix[-problem_genes_idx, , drop = FALSE]
      if (nrow(data_matrix) == 0) {
        stop("All genes removed after filtering. Check your data - DESeq2 requires raw integer counts, not normalized values (FPKM/TPM/RPKM).")
      }
    } else {
      # Stop with detailed error message
      genes_str <- paste(head(problematic_gene_names, 5), collapse = ", ")
      samples_str <- paste(problematic_sample_names, collapse = ", ")
      stop(err_infinite_values(problem_genes_idx, problem_samples_idx))
    }
  } else {
    # No issues found, use original data
    data_matrix <- data_matrix
  }

  # Convert to integer storage
  storage.mode(data_matrix) <- "integer"

  # Validate sample-column correspondence
  if (ncol(data_matrix) != nrow(sample)) {
    stop("Number of columns in 'data' must equal number of rows in 'sample'")
  }

  # Check sample names matching
  data_samples <- colnames(data_matrix)
  sample_names <- rownames(sample)

  # Auto-detect sample ID column in sample (if no rownames or default rownames)
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
    missing_samples <- setdiff(data_samples, sample_names)
    stop(err_row_col_mismatch("data", "sample metadata"))
  }

  if (!all(sample_names %in% data_samples)) {
    extra_samples <- setdiff(sample_names, data_samples)
    warn_data_processing(paste("Extra samples in metadata (ignored):", paste(extra_samples, collapse = ", ")))
  }

  # Align sample metadata with count data
  sample_aligned <- sample[data_samples, , drop = FALSE]

  # Validate design formula
  if (!inherits(formula, "formula")) {
    stop(err_invalid_input("formula", "a formula object (e.g., ~condition)"))
  }

  validate_formula_vars(formula, sample_aligned)

  formula_vars <- all.vars(formula)
  if (length(formula_vars) == 0) {
    stop(err_invalid_input("formula", "at least one variable"))
  }

  # Check for missing values in design variables
  for (var in formula_vars) {
    if (any(is.na(sample_aligned[[var]]))) {
      stop(err_missing_values(sprintf("design variable '%s'", var), sum(is.na(sample_aligned[[var]]))))
    }
  }

  # Validate threshold parameters
  if (!is.numeric(log2FoldChange) || length(log2FoldChange) != 1 || log2FoldChange < 0) {
    stop("'log2FoldChange' must be a single non-negative number")
  }

  if (!is.numeric(padj) || length(padj) != 1 || padj <= 0 || padj > 1) {
    stop("'padj' must be a single number between 0 and 1")
  }

  if (!is.logical(shrink.lfc) || length(shrink.lfc) != 1) {
    stop("'shrink.lfc' must be a single logical value")
  }

  if (!is.logical(independent.filtering) || length(independent.filtering) != 1) {
    stop("'independent.filtering' must be a single logical value")
  }

  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single number between 0 and 1")
  }

  # Check experimental design adequacy
  main_factor <- formula_vars[length(formula_vars)]  # Usually the last term is main effect
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

  # Check for genes with zero counts across all samples
  zero_count_genes <- rowSums(data_matrix) == 0
  n_zero_genes <- sum(zero_count_genes)

  if (n_zero_genes > 0) {
    warning(paste("Found", n_zero_genes, "genes with zero counts across all samples. These will be filtered."))
    if (n_zero_genes == nrow(data_matrix)) {
      stop("All genes have zero counts. Check data quality.")
    }
  }

  # Check library sizes for extreme differences
  lib_sizes <- colSums(data_matrix)
  lib_size_ratio <- max(lib_sizes) / min(lib_sizes)

  if (lib_size_ratio > 10) {
    warning(paste("Large library size differences detected (max/min ratio:",
      round(lib_size_ratio, 1), "). Consider checking data quality."))
  }

  if (min(lib_sizes) < 1000) {
    warning("Very small library sizes detected (<1000 counts). Results may be unreliable.")
  }

  # --- Pairwise mode --------------------------------------------------------
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

    # Remove zero-count genes globally (once)
    data_matrix <- data_matrix[!zero_count_genes, , drop = FALSE]

    run_one_pair <- function(pair) {
      idx <- sample_aligned[[main_factor]] %in% pair
      d_sub <- data_matrix[, idx, drop = FALSE]
      s_sub <- sample_aligned[idx, , drop = FALSE]
      s_sub[[main_factor]] <- factor(as.character(s_sub[[main_factor]]), levels = pair)

      if (any(table(s_sub[[main_factor]]) < 2)) return(NULL)

      pair_formula <- as.formula(paste0("~", main_factor))
      comp_label <- paste(pair[2], "vs", pair[1])

      dds <- tryCatch(
        DESeq2::DESeqDataSetFromMatrix(countData = d_sub, colData = s_sub, design = pair_formula),
        error = function(e) NULL
      )
      if (is.null(dds)) return(NULL)

      dds <- tryCatch(DESeq2::DESeq(dds, quiet = TRUE), error = function(e) NULL)
      if (is.null(dds)) return(NULL)

      res <- tryCatch({
        if (shrink.lfc) {
          res_raw <- tryCatch(
            DESeq2::lfcShrink(dds,
              coef = DESeq2::resultsNames(dds)[length(DESeq2::resultsNames(dds))],
              type = "apeglm", quiet = TRUE),
            error = function(e) DESeq2::results(dds, alpha = alpha,
              independentFiltering = independent.filtering))
        } else {
          res_raw <- DESeq2::results(dds, alpha = alpha,
            independentFiltering = independent.filtering)
        }
        as.data.frame(res_raw)
      }, error = function(e) NULL)
      if (is.null(res)) return(NULL)

      res <- tibble::rownames_to_column(res, "gene")
      res$padj <- ifelse(is.na(res$padj), 1, res$padj)
      res$regulation <- ifelse(res$log2FoldChange > log2FoldChange & res$padj < padj, "Up",
                        ifelse(res$log2FoldChange < -log2FoldChange & res$padj < padj, "Down", "NS"))
      res$fold_change <- 2^res$log2FoldChange
      res$abs_log2fc <- abs(res$log2FoldChange)
      res$comparison <- comp_label
      res$group <- pair[2]
      res$ref_group <- pair[1]

      # Add group descriptive statistics
      idx_grp <- s_sub[[main_factor]] == pair[2]
      idx_ref <- s_sub[[main_factor]] == pair[1]
      n_grp <- sum(idx_grp)
      n_ref <- sum(idx_ref)

      # Normalized counts for mean/sd
      norm_counts <- DESeq2::counts(dds, normalized = TRUE)
      grp_means <- rowMeans(norm_counts[, idx_grp, drop = FALSE], na.rm = TRUE)
      grp_sds <- apply(norm_counts[, idx_grp, drop = FALSE], 1, sd, na.rm = TRUE)
      ref_means <- rowMeans(norm_counts[, idx_ref, drop = FALSE], na.rm = TRUE)
      ref_sds <- apply(norm_counts[, idx_ref, drop = FALSE], 1, sd, na.rm = TRUE)

      res$group_mean <- round(grp_means[res$gene], 2)
      res$group_sd <- round(grp_sds[res$gene], 2)
      res$group_n <- n_grp
      res$ref_mean <- round(ref_means[res$gene], 2)
      res$ref_sd <- round(ref_sds[res$gene], 2)
      res$ref_n <- n_ref
      res$test_method <- "DESeq2-Wald"
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

    # Join original expression values
    combined <- dplyr::left_join(combined, expr_df, by = "gene")
    combined <- combined[order(combined$comparison, combined$padj, -combined$abs_log2fc), ]

    if (interactive()) {
      cat("DESeq2 Pairwise Differential Expression Summary:\n")
      cat("=================================================\n")
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

  # --- Single comparison (original logic) -----------------------------------

  # Create DESeqDataSet
  tryCatch(
    {
      dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = data_matrix,
        colData = sample_aligned,
        design = formula
      )
    },
    error = function(e) {
      stop(paste("Error creating DESeqDataSet:", e$message))
    })

  # Run DESeq2 analysis
  tryCatch(
    {
      dds_analyzed <- DESeq2::DESeq(dds, quiet = TRUE)
    },
    error = function(e) {
      stop(paste("Error in DESeq2 analysis:", e$message,
        "\nThis may be due to insufficient replicates or design matrix issues."))
    })

  # # Extract results with specified parameters
  # tryCatch({
  #   if (shrink.lfc) {
  #     # Apply LFC shrinkage for more accurate effect size estimates
  #     results_raw <- DESeq2::lfcShrink(
  #       dds_analyzed,
  #       coef = DESeq2::resultsNames(dds_analyzed)[length(DESeq2::resultsNames(dds_analyzed))],
  #       alpha = alpha,
  #       type = "apeglm",  # Recommended shrinkage method
  #       quiet = TRUE
  #     )
  #   } else {
  #     results_raw <- DESeq2::results(
  #       dds_analyzed,
  #       alpha = alpha,
  #       independentFiltering = independent.filtering
  #     )
  #   }

  #   results_df <- as.data.frame(results_raw)
  # }, error = function(e) {
  #   # Fallback to basic results if shrinkage fails
  #   warning("LFC shrinkage failed, using unshrunken results:", e$message)
  #   results_raw <- DESeq2::results(
  #     dds_analyzed,
  #     alpha = alpha,
  #     independentFiltering = independent.filtering
  #   )
  #   results_df <- as.data.frame(results_raw)
  # })
  # 将 tryCatch 的结果赋值给 results_df
  results_df <- tryCatch(
    {
      # 'try' 代码块，用于正常流程
      if (shrink.lfc) {
        # 尝试进行 LFC shrinkage
        message("Attempting LFC shrinkage with apeglm...")
        results_raw <- DESeq2::lfcShrink(
          dds_analyzed,
          coef = DESeq2::resultsNames(dds_analyzed)[length(DESeq2::resultsNames(dds_analyzed))],
          type = "apeglm",
          quiet = TRUE
        )
        message("LFC shrinkage successful.")
      } else {
        # 如果不进行 shrinkage
        results_raw <- DESeq2::results(
          dds_analyzed,
          alpha = alpha,
          independentFiltering = independent.filtering
        )
      }
      # 如果 try 成功，将结果转换为数据框并返回
      as.data.frame(results_raw)

    },
    error = function(e) {
      # 'error' 函数，在 try 失败时执行
      # 打印一个更详细的警告信息
      warning(paste("LFC shrinkage failed, using unshrunken results instead. Reason:", e$message,
        "\nThis can sometimes be caused by package conflicts or issues with the DESeqDataSet object."))

      # 执行备用方案
      results_raw <- DESeq2::results(
        dds_analyzed,
        alpha = alpha,
        independentFiltering = independent.filtering
      )
      # 将备用方案的结果转换为数据框并返回
      as.data.frame(results_raw)
    })

  # Process and enhance results
  results_processed <- results_df %>%
    tibble::rownames_to_column(var = "gene") %>%
    dplyr::mutate(
      # Main classification
      regulation = dplyr::case_when(
        log2FoldChange > !!log2FoldChange & padj < !!padj ~ "Up-regulated",
        log2FoldChange < -!!log2FoldChange & padj < !!padj ~ "Down-regulated",
        TRUE ~ "Not significant"
      ),
      # Additional useful columns
      fold_change = 2^log2FoldChange,
      abs_log2fc = abs(log2FoldChange),
      # Handle missing padj values
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
    # Join original expression values
    dplyr::left_join(expr_df, by = "gene") %>%
    # Sort by significance and effect size
    dplyr::arrange(padj, desc(abs_log2fc))

  # Calculate summary statistics
  n_upregulated <- sum(results_processed$regulation == "Up-regulated", na.rm = TRUE)
  n_downregulated <- sum(results_processed$regulation == "Down-regulated", na.rm = TRUE)
  n_total_de <- n_upregulated + n_downregulated
  n_tested <- sum(!is.na(results_processed$padj))

  # Add comprehensive metadata
  attr(results_processed, "n_genes_input") <- nrow(data_matrix)
  attr(results_processed, "n_genes_original") <- n_genes_original
  attr(results_processed, "n_genes_tested") <- n_tested
  attr(results_processed, "n_samples") <- ncol(data_matrix)
  attr(results_processed, "design_formula") <- deparse(formula)
  attr(results_processed, "log2fc_threshold") <- log2FoldChange
  attr(results_processed, "padj_threshold") <- padj
  attr(results_processed, "shrinkage_applied") <- shrink.lfc
  attr(results_processed, "independent_filtering") <- independent.filtering
  attr(results_processed, "alpha_level") <- alpha
  attr(results_processed, "library_size_range") <- range(lib_sizes)
  attr(results_processed, "remove_problematic_genes") <- remove_problematic_genes
  if (has_nonfinite && remove_problematic_genes) {
    attr(results_processed, "n_genes_removed") <- length(problematic_gene_names)
    attr(results_processed, "genes_removed") <- problematic_gene_names
  }

  # Interactive summary
  if (interactive()) {
    cat("DESeq2 Differential Expression Analysis Summary:\n")
    cat("===============================================\n")

    # Show gene filtering information
    if (has_nonfinite && remove_problematic_genes) {
      cat("Genes (original):", n_genes_original, "\n")
      cat("Genes removed (non-finite values):", length(problematic_gene_names), "\n")
      cat("Genes analyzed:", nrow(data_matrix), "\n")
    } else {
      cat("Genes input:", nrow(data_matrix), "\n")
    }
    cat("Genes tested:", n_tested, "\n")
    cat("Samples analyzed:", ncol(data_matrix), "\n")
    cat("Design formula:", deparse(formula), "\n")
    cat("Log2FC threshold:", log2FoldChange, "(", round(2^log2FoldChange, 2), "-fold)\n")
    cat("FDR threshold:", padj, "\n")
    cat("LFC shrinkage:", ifelse(shrink.lfc, "Applied", "Not applied"), "\n")
    cat("Independent filtering:", ifelse(independent.filtering, "Enabled", "Disabled"), "\n\n")

    # Library size summary
    cat("Library size summary:\n")
    cat("  Range:", min(lib_sizes), "to", max(lib_sizes), "counts\n")
    cat("  Median:", median(lib_sizes), "counts\n")

    if (main_factor %in% names(sample_aligned)) {
      cat("  Group sizes:", paste(names(group_counts), "=", group_counts, collapse = ", "), "\n\n")
    }

    # Results summary
    cat("Results:\n")
    cat("  Up-regulated genes:", n_upregulated, "\n")
    cat("  Down-regulated genes:", n_downregulated, "\n")
    cat("  Not significant:", sum(results_processed$regulation == "Not significant"), "\n")
    cat("  Total DE genes:", n_total_de,
      paste0("(", round(100 * n_total_de / n_tested, 1), "% of tested)"), "\n\n")

    if (n_total_de > 0) {
      cat("Top DE genes by significance:\n")
      top_genes <- results_processed %>%
        dplyr::filter(regulation != "Not significant") %>%
        head(5)

      for (i in 1:nrow(top_genes)) {
        direction_symbol <- ifelse(top_genes$log2FoldChange[i] > 0, "↑", "↓")
        cat(sprintf("  %s %s: %s (log2FC = %.2f, padj = %.2e)\n",
          direction_symbol,
          top_genes$gene[i],
          top_genes$regulation[i],
          top_genes$log2FoldChange[i],
          top_genes$padj[i]))
      }

      # Effect size distribution
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
