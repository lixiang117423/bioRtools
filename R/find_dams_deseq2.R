#' Identify Differentially Abundant Microorganisms/Taxa Using DESeq2
#'
#' This function performs differential abundance analysis on microbiome count data
#' using DESeq2's negative binomial model. It identifies taxa (OTUs, ASVs, or other
#' features) that are significantly enriched or depleted between different conditions
#' or groups. DESeq2 is particularly well-suited for count data with overdispersion
#' commonly observed in microbiome datasets.
#'
#' @param data Integer count matrix where rows represent microbial features (OTUs, ASVs,
#'   taxa, etc.) and columns represent samples. Row names should be feature identifiers
#'   and column names should match sample names in the sample metadata. All values
#'   must be non-negative integers (raw counts, not normalized or transformed)
#' @param sample Data frame containing sample metadata where rows represent samples
#'   and columns contain experimental factors. Row names should match column names
#'   in the count data. Must contain at least one grouping variable for comparison
#' @param formula Design formula specifying the experimental design for DESeq2
#'   analysis. Default is \code{~group}, assuming a column named "group" exists
#'   in the sample metadata. Examples:
#'   \itemize{
#'     \item \code{~treatment}: Simple comparison between treatments
#'     \item \code{~treatment + batch}: Account for batch effects
#'     \item \code{~treatment + time}: Include time as a factor
#'   }
#' @param log2FoldChange Minimum absolute log2 fold change threshold for
#'   significance calling (default: 1, equivalent to 2-fold change). Features
#'   with |log2FC| >= this value and padj < padj threshold are considered
#'   differentially abundant
#' @param padj Adjusted p-value (FDR) threshold for significance (default: 0.05).
#'   Features with adjusted p-values below this threshold are considered significant
#' @param min.count.per.sample Minimum count per sample for feature filtering
#'   (default: 1). Features with counts below this threshold in all samples may
#'   be filtered out to improve statistical power
#' @param min.samples.with.counts Minimum number of samples that must have counts
#'   above min.count.per.sample for a feature to be retained (default: 2)
#'
#' @return A data frame containing differential abundance results with columns:
#'   \itemize{
#'     \item \code{feature_id}: Feature identifiers (formerly OTU column)
#'     \item \code{baseMean}: Mean normalized counts across all samples
#'     \item \code{log2FoldChange}: Log2 fold change between conditions
#'     \item \code{lfcSE}: Standard error of log2 fold change
#'     \item \code{stat}: Wald statistic
#'     \item \code{pvalue}: Raw p-value from Wald test
#'     \item \code{padj}: Benjamini-Hochberg adjusted p-value (FDR)
#'     \item \code{significance}: Significance classification:
#'       \itemize{
#'         \item "Enriched": log2FC > threshold & padj < threshold
#'         \item "Depleted": log2FC < -threshold & padj < threshold  
#'         \item "NS": Not significant
#'       }
#'     \item \code{fold_change}: Actual fold change (2^log2FoldChange)
#'     \item \code{abs_log2fc}: Absolute log2 fold change for ranking
#'   }
#'
#' @details
#' This function implements the DESeq2 workflow for microbiome data:
#' \enumerate{
#'   \item Creates DESeqDataSet from count matrix and sample metadata
#'   \item Estimates size factors for normalization
#'   \item Estimates dispersions (gene-wise and fitted)
#'   \item Fits negative binomial generalized linear models
#'   \item Performs Wald tests for differential abundance
#'   \item Applies Benjamini-Hochberg correction for multiple testing
#' }
#' 
#' DESeq2 assumptions:
#' \itemize{
#'   \item Count data follows negative binomial distribution
#'   \item Features are independent
#'   \item Samples are independent
#'   \item Sequencing depth varies between samples (handled by normalization)
#' }
#'
#' @note
#' \itemize{
#'   \item Input data must be raw counts (integers), not relative abundances
#'   \item Very sparse features (present in few samples) may not test reliably
#'   \item Consider rarefying or filtering low-abundance features before analysis
#'   \item DESeq2 internally normalizes for library size differences
#'   \item For complex designs, ensure sufficient samples per group
#' }
#'
#' @references
#' Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold change
#' and dispersion for RNA-seq data with DESeq2. Genome Biology, 15, 550.
#' 
#' McMurdie, P.J. and Holmes, S. (2014) Waste not, want not: why rarefying
#' microbiome data is inadmissible. PLoS Computational Biology, 10(4), e1003531.
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{find_dams_lefse}} for LEfSe-based differential abundance analysis
#'   \item \code{\link{rarefy_table}} for rarefaction of count data
#'   \item \code{\link{volcano_plot}} for visualization of results
#' }
#'
#' @export
#'
#' @examples
#' library(bioRtools)
#' library(dplyr)
#' 
#' # Example 1: Basic differential abundance analysis
#' \dontrun{
#' # Load example microbiome data
#' data(df.pcoa.otu)     # OTU count table
#' data(df.pcoa.sample)  # Sample metadata
#' 
#' # Transpose count data (DESeq2 expects features as rows)
#' count_matrix <- t(df.pcoa.otu)
#' 
#' # Basic differential abundance analysis
#' dam_results <- find_dams_deseq2(
#'   data = count_matrix,
#'   sample = df.pcoa.sample
#' )
#' 
#' # View results summary
#' print(table(dam_results$significance))
#' 
#' # Top enriched features
#' enriched <- dam_results %>%
#'   filter(significance == "Enriched") %>%
#'   arrange(padj) %>%
#'   head(10)
#' print(enriched[, c("feature_id", "log2FoldChange", "padj", "baseMean")])
#' }
#' 
#' # Example 2: Custom thresholds and design
#' \dontrun{
#' # More stringent analysis
#' dam_strict <- find_dams_deseq2(
#'   data = t(df.pcoa.otu),
#'   sample = df.pcoa.sample,
#'   formula = ~treatment + batch,  # Account for batch effects
#'   log2FoldChange = 1.5,          # 3-fold change threshold
#'   padj = 0.01,                   # 1% FDR threshold
#'   min.count.per.sample = 5,      # Higher count threshold
#'   min.samples.with.counts = 3    # Present in at least 3 samples
#' )
#' 
#' print(paste("Strict analysis found", sum(dam_strict$significance != "NS"), 
#'             "significant features"))
#' }
#' 
#' # Example 3: Time series analysis
#' \dontrun{
#' # Assume sample metadata has 'timepoint' and 'subject' columns
#' # For longitudinal microbiome data
#' longitudinal_results <- find_dams_deseq2(
#'   data = t(df.pcoa.otu),
#'   sample = df.pcoa.sample,
#'   formula = ~timepoint + subject,  # Account for subject variation
#'   log2FoldChange = 0.5,            # More sensitive threshold
#'   padj = 0.1                       # Less stringent for exploratory analysis
#' )
#' 
#' # Focus on taxa that change over time
#' time_responsive <- longitudinal_results %>%
#'   filter(significance != "NS") %>%
#'   arrange(abs_log2fc) %>%
#'   tail(20)  # Top 20 by effect size
#' 
#' print("Top time-responsive taxa:")
#' print(time_responsive[, c("feature_id", "log2FoldChange", "padj")])
#' }
#' 
#' # Example 4: Creating simple count data for testing
#' # This demonstrates the required data format
#' set.seed(123)
#' n_features <- 100
#' n_samples <- 20
#' 
#' # Create mock count data
#' mock_counts <- matrix(
#'   rnbinom(n_features * n_samples, size = 1, mu = 10),
#'   nrow = n_features,
#'   ncol = n_samples,
#'   dimnames = list(
#'     paste0("OTU_", 1:n_features),
#'     paste0("Sample_", 1:n_samples)
#'   )
#' )
#' 
#' # Create mock sample metadata
#' mock_sample <- data.frame(
#'   row.names = colnames(mock_counts),
#'   group = rep(c("Control", "Treatment"), each = 10),
#'   batch = rep(c("Batch1", "Batch2"), times = 10),
#'   subject_id = paste0("Subject_", 1:20)
#' )
#' 
#' \dontrun{
#' # Run analysis on mock data
#' mock_results <- find_dams_deseq2(
#'   data = mock_counts,
#'   sample = mock_sample,
#'   formula = ~group
#' )
#' 
#' print("Mock data analysis results:")
#' print(table(mock_results$significance))
#' }
#'
find_dams_deseq2 <- function(data, sample, formula = ~group, log2FoldChange = 1, padj = 0.05,
                             min.count.per.sample = 1, min.samples.with.counts = 2) {
  
  # Input validation
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("'data' must be a matrix or data frame")
  }
  
  if (!is.data.frame(sample)) {
    stop("'sample' must be a data frame")
  }
  
  if (nrow(data) == 0 || ncol(data) == 0) {
    stop("'data' cannot be empty")
  }
  
  if (nrow(sample) == 0) {
    stop("'sample' cannot be empty")
  }
  
  # Check if data contains only non-negative integers
  if (any(data < 0, na.rm = TRUE)) {
    stop("'data' must contain only non-negative values (counts)")
  }
  
  if (any(!is.finite(as.matrix(data)), na.rm = TRUE)) {
    stop("'data' contains non-finite values (Inf, -Inf, NaN)")
  }
  
  # Convert to integer matrix if not already
  if (!is.integer(as.matrix(data)[1,1])) {
    warning("Converting count data to integers. Ensure input data represents raw counts.")
    data <- round(as.matrix(data))
    storage.mode(data) <- "integer"
  }
  
  # Check sample-column correspondence
  if (ncol(data) != nrow(sample)) {
    stop("Number of columns in 'data' must equal number of rows in 'sample'")
  }
  
  # Check if sample names match between data and metadata
  data_samples <- colnames(data)
  sample_names <- rownames(sample)
  
  if (is.null(data_samples) || is.null(sample_names)) {
    stop("Both 'data' (columns) and 'sample' (rows) must have names")
  }
  
  if (!all(data_samples %in% sample_names)) {
    missing_samples <- setdiff(data_samples, sample_names)
    stop(paste("Sample(s) missing from metadata:", paste(missing_samples, collapse = ", ")))
  }
  
  if (!all(sample_names %in% data_samples)) {
    extra_samples <- setdiff(sample_names, data_samples)
    warning(paste("Extra samples in metadata (will be ignored):", paste(extra_samples, collapse = ", ")))
  }
  
  # Align sample metadata with count data
  sample_aligned <- sample[data_samples, , drop = FALSE]
  
  # Validate formula design variables
  formula_vars <- all.vars(formula)
  missing_vars <- setdiff(formula_vars, names(sample_aligned))
  if (length(missing_vars) > 0) {
    stop(paste("Formula variables not found in sample metadata:", paste(missing_vars, collapse = ", ")))
  }
  
  # Check for missing values in design variables
  for (var in formula_vars) {
    if (any(is.na(sample_aligned[[var]]))) {
      stop(paste("Missing values found in design variable:", var))
    }
  }
  
  # Validate thresholds
  if (!is.numeric(log2FoldChange) || length(log2FoldChange) != 1 || log2FoldChange < 0) {
    stop("'log2FoldChange' must be a single non-negative number")
  }
  
  if (!is.numeric(padj) || length(padj) != 1 || padj <= 0 || padj > 1) {
    stop("'padj' must be a single number between 0 and 1")
  }
  
  if (!is.numeric(min.count.per.sample) || min.count.per.sample < 0) {
    stop("'min.count.per.sample' must be non-negative")
  }
  
  if (!is.numeric(min.samples.with.counts) || min.samples.with.counts < 1) {
    stop("'min.samples.with.counts' must be at least 1")
  }
  
  # Optional filtering of low-count features
  if (min.count.per.sample > 0 && min.samples.with.counts > 0) {
    feature_filter <- rowSums(data >= min.count.per.sample) >= min.samples.with.counts
    n_filtered <- sum(!feature_filter)
    
    if (n_filtered > 0) {
      warning(paste("Filtering", n_filtered, "low-count features"))
      data <- data[feature_filter, , drop = FALSE]
    }
    
    if (nrow(data) == 0) {
      stop("No features remain after filtering. Consider reducing filter thresholds.")
    }
  }
  
  # Check for sufficient features
  if (nrow(data) < 2) {
    stop("At least 2 features required for differential abundance analysis")
  }
  
  # Create DESeqDataSet
  tryCatch({
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = data,
      colData = sample_aligned,
      design = formula
    )
  }, error = function(e) {
    stop(paste("Error creating DESeqDataSet:", e$message))
  })
  
  # Check for sufficient samples per group
  design_vars <- all.vars(formula)
  main_factor <- design_vars[1]  # Assume first variable is main factor
  
  if (main_factor %in% names(sample_aligned)) {
    group_counts <- table(sample_aligned[[main_factor]])
    if (any(group_counts < 2)) {
      warning(paste("Some groups have fewer than 2 samples. Results may be unreliable."))
    }
  }
  
  # Run DESeq2 analysis
  tryCatch({
    dds_analyzed <- DESeq2::DESeq(dds, quiet = TRUE)
  }, error = function(e) {
    stop(paste("Error in DESeq2 analysis:", e$message))
  })
  
  # Extract results
  tryCatch({
    results_raw <- DESeq2::results(dds_analyzed)
    results_df <- as.data.frame(results_raw)
  }, error = function(e) {
    stop(paste("Error extracting DESeq2 results:", e$message))
  })
  
  # Process results
  results_processed <- results_df %>%
    tibble::rownames_to_column(var = "feature_id") %>%
    dplyr::mutate(
      # Classification based on thresholds
      significance = dplyr::case_when(
        log2FoldChange > !!log2FoldChange & padj < !!padj ~ "Enriched",
        log2FoldChange < -!!log2FoldChange & padj < !!padj ~ "Depleted", 
        TRUE ~ "NS"
      ),
      # Additional useful columns
      fold_change = 2^log2FoldChange,
      abs_log2fc = abs(log2FoldChange),
      # Handle missing padj values
      padj = ifelse(is.na(padj), 1, padj)
    ) %>%
    # Arrange by significance and effect size
    dplyr::arrange(padj, desc(abs_log2fc))
  
  # Add analysis metadata as attributes
  attr(results_processed, "n_features_input") <- nrow(data) + ifelse(exists("n_filtered"), n_filtered, 0)
  attr(results_processed, "n_features_tested") <- nrow(data)
  attr(results_processed, "n_samples") <- ncol(data)
  attr(results_processed, "design_formula") <- deparse(formula)
  attr(results_processed, "log2fc_threshold") <- log2FoldChange
  attr(results_processed, "padj_threshold") <- padj
  attr(results_processed, "filtering_applied") <- min.count.per.sample > 0
  
  # Generate summary statistics
  n_enriched <- sum(results_processed$significance == "Enriched", na.rm = TRUE)
  n_depleted <- sum(results_processed$significance == "Depleted", na.rm = TRUE)
  n_total_sig <- n_enriched + n_depleted
  
  # Print summary if interactive
  if (interactive()) {
    cat("DESeq2 Differential Abundance Analysis Summary:\n")
    cat("==============================================\n")
    cat("Features tested:", nrow(results_processed), "\n")
    cat("Samples analyzed:", ncol(data), "\n")
    cat("Design formula:", deparse(formula), "\n")
    cat("Log2FC threshold:", log2FoldChange, "(", round(2^log2FoldChange, 2), "-fold)\n")
    cat("FDR threshold:", padj, "\n\n")
    
    cat("Results:\n")
    cat("  Enriched features:", n_enriched, "\n")
    cat("  Depleted features:", n_depleted, "\n") 
    cat("  Non-significant:", sum(results_processed$significance == "NS", na.rm = TRUE), "\n")
    cat("  Total significant:", n_total_sig, 
        paste0("(", round(100 * n_total_sig / nrow(results_processed), 1), "%)"), "\n\n")
    
    if (n_total_sig > 0) {
      cat("Top significant features by effect size:\n")
      top_features <- results_processed %>%
        dplyr::filter(significance != "NS") %>%
        dplyr::arrange(desc(abs_log2fc)) %>%
        head(5)
      
      for (i in 1:nrow(top_features)) {
        cat(sprintf("%s: %s (log2FC = %.2f, padj = %.2e)\n",
                    top_features$feature_id[i],
                    top_features$significance[i], 
                    top_features$log2FoldChange[i],
                    top_features$padj[i]))
      }
    }
    cat("\n")
  }
  
  return(results_processed)
}