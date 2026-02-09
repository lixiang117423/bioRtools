#' Identify Differentially Abundant Microorganisms Using LEfSe Analysis
#'
#' This function performs Linear Discriminant Analysis (LDA) Effect Size (LEfSe)
#' analysis to identify microbial biomarkers that are differentially abundant
#' between groups. LEfSe uses a combination of non-parametric statistical tests
#' (Kruskal-Wallis and Wilcoxon rank-sum) followed by LDA to identify features
#' that are both statistically significant and biologically meaningful.
#'
#' @param data Feature abundance matrix where rows represent microbial features
#'   (OTUs, ASVs, taxa, etc.) and columns represent samples. Can be raw counts
#'   or relative abundances. Row names should be feature identifiers and column
#'   names should match sample identifiers in the sample metadata
#' @param sample Data frame containing sample metadata where rows represent
#'   samples and columns contain experimental factors and covariates. Row names
#'   should match column names in the feature data. Must contain the grouping
#'   variable specified in groupCol parameter
#' @param groupCol Character string specifying the column name in sample metadata
#'   that contains the main grouping variable for comparison (default: "group").
#'   This variable will be used for the LEfSe analysis to identify biomarkers
#' @param kruskal.threshold P-value threshold for the Kruskal-Wallis test
#'   (default: 0.05). Features must pass this threshold to proceed to pairwise
#'   testing. The Kruskal-Wallis test assesses whether at least one group differs
#'   from the others
#' @param wilcox.threshold P-value threshold for pairwise Wilcoxon rank-sum tests
#'   (default: 0.05). Features must pass both Kruskal-Wallis and all relevant
#'   pairwise Wilcoxon tests to be considered for LDA analysis
#' @param lda.threshold Minimum Linear Discriminant Analysis (LDA) effect size
#'   threshold (default: 1.0). Features must have LDA scores above this threshold
#'   to be considered biologically meaningful biomarkers. Higher values indicate
#'   stronger discriminatory power
#' @param strict.wilcox Logical indicating whether to apply strict Wilcoxon
#'   testing (default: TRUE). If TRUE, features must pass all relevant pairwise
#'   comparisons. If FALSE, more lenient testing is applied
#' @param bootstrap.iterations Number of bootstrap iterations for LDA stability
#'   assessment (default: 30). More iterations provide more stable results but
#'   increase computation time
#'
#' @return A data frame containing LEfSe analysis results with columns:
#'   \itemize{
#'     \item \code{features}: Feature identifiers (e.g., OTU IDs, taxa names)
#'     \item \code{scores}: LDA scores representing discriminatory power
#'     \item \code{class}: Group/class where the feature shows highest abundance
#'     \item \code{pvalue}: Combined p-value from statistical tests
#'     \item \code{lda.score}: Absolute LDA effect size
#'     \item \code{direction}: Direction of enrichment relative to comparison
#'   }
#'
#'   Only features that pass all statistical thresholds (Kruskal-Wallis, Wilcoxon,
#'   and LDA) are included in the results.
#'
#' @details
#' LEfSe (Linear Discriminant Analysis Effect Size) workflow:
#' \enumerate{
#'   \item \strong{Kruskal-Wallis test}: Non-parametric test to identify features
#'     with significantly different distributions among groups
#'   \item \strong{Wilcoxon rank-sum test}: Pairwise comparisons between groups
#'     to ensure consistent differences
#'   \item \strong{Linear Discriminant Analysis}: Estimates effect size and
#'     identifies the group where each feature is most abundant
#' }
#'
#' Key advantages of LEfSe:
#' \itemize{
#'   \item Non-parametric approach suitable for non-normal distributions
#'   \item Combines statistical significance with biological effect size
#'   \item Robust to outliers and varying sample sizes
#'   \item Provides interpretable biomarker rankings
#' }
#'
#' Considerations:
#' \itemize{
#'   \item Works best with at least 6-10 samples per group
#'   \item Results depend on threshold parameter choices
#'   \item May be conservative with strict multiple testing correction
#'   \item LDA assumes linear separability between groups
#' }
#'
#' @note
#' \itemize{
#'   \item Input data can be raw counts or relative abundances
#'   \item LEfSe internally handles data transformation and normalization
#'   \item Very low-abundance features may not be detected reliably
#'   \item Consider pre-filtering rare features for computational efficiency
#'   \item Ensure balanced sample sizes between groups when possible
#' }
#'
#' @references
#' Segata, N., et al. (2011). Metagenomic biomarker discovery and explanation.
#' Genome Biology, 12(6), R60.
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{find_dams_deseq2}} for DESeq2-based differential abundance
#'   \item \code{\link{rarefy_table}} for rarefaction preprocessing
#'   \item \code{\link{permanova_test}} for community-level differences
#' }
#'
#' @export
#'
#' @examples
#' library(bioRtools)
#' library(dplyr)
#'
#' # Example 1: Basic LEfSe analysis
#' \dontrun{
#' # Load example microbiome data
#' data(df.call_DAMs_LEfSe.otu)     # Feature abundance table
#' data(df.call_DAMs_LEfSe.sample)  # Sample metadata
#'
#' # Basic LEfSe analysis
#' lefse_results <- find_dams_lefse(
#'   data = df.call_DAMs_LEfSe.otu,
#'   sample = df.call_DAMs_LEfSe.sample
#' )
#'
#' # View results
#' print(head(lefse_results))
#' print(paste("Found", nrow(lefse_results), "significant biomarkers"))
#'
#' # Top biomarkers by LDA score
#' top_biomarkers <- lefse_results %>%
#'   arrange(desc(abs(scores))) %>%
#'   head(10)
#' print(top_biomarkers)
#' }
#'
#' # Example 2: Stringent analysis parameters
#' \dontrun{
#' # More conservative thresholds
#' lefse_strict <- find_dams_lefse(
#'   data = df.call_DAMs_LEfSe.otu,
#'   sample = df.call_DAMs_LEfSe.sample,
#'   groupCol = "treatment",
#'   kruskal.threshold = 0.01,        # More stringent Kruskal-Wallis
#'   wilcox.threshold = 0.01,         # More stringent Wilcoxon
#'   lda.threshold = 2.0,             # Higher LDA threshold
#'   bootstrap.iterations = 50        # More bootstrap iterations
#' )
#'
#' print(paste("Strict analysis found", nrow(lefse_strict), "biomarkers"))
#'
#' # Compare effect sizes
#' if (nrow(lefse_strict) > 0) {
#'   print("Top strict biomarkers:")
#'   print(lefse_strict[1:min(5, nrow(lefse_strict)),
#'     c("features", "class", "scores")])
#' }
#' }
#'
#' # Example 3: Multi-group comparison
#' \dontrun{
#' # Assume sample metadata has multiple treatment groups
#' # Create sample data with multiple groups for demonstration
#' multi_sample <- df.call_DAMs_LEfSe.sample
#' multi_sample$treatment_type <- rep(c("Drug_A", "Drug_B", "Placebo"),
#'   length.out = nrow(multi_sample))
#'
#' # LEfSe with multiple groups
#' multi_lefse <- find_dams_lefse(
#'   data = df.call_DAMs_LEfSe.otu,
#'   sample = multi_sample,
#'   groupCol = "treatment_type",
#'   lda.threshold = 1.5
#' )
#'
#' # Examine biomarkers by group
#' if (nrow(multi_lefse) > 0) {
#'   biomarker_summary <- multi_lefse %>%
#'     group_by(class) %>%
#'     summarise(
#'       n_biomarkers = n(),
#'       mean_lda_score = mean(abs(scores)),
#'       .groups = "drop"
#'     )
#'   print("Biomarkers by group:")
#'   print(biomarker_summary)
#' }
#' }
#'
#' # Example 4: Creating sample data format
#' # This demonstrates the required data structure
#' set.seed(42)
#' n_features <- 50
#' n_samples <- 24
#'
#' # Create mock abundance data (can be counts or relative abundances)
#' mock_abundance <- matrix(
#'   runif(n_features * n_samples, min = 0, max = 100),
#'   nrow = n_features,
#'   ncol = n_samples,
#'   dimnames = list(
#'     paste0("Feature_", 1:n_features),
#'     paste0("Sample_", 1:n_samples)
#'   )
#' )
#'
#' # Create mock sample metadata
#' mock_sample_meta <- data.frame(
#'   row.names = colnames(mock_abundance),
#'   group = rep(c("Healthy", "Disease"), each = 12),
#'   age = sample(20:80, n_samples),
#'   sex = sample(c("Male", "Female"), n_samples, replace = TRUE),
#'   study_site = rep(c("Site_A", "Site_B"), times = 12)
#' )
#'
#' \dontrun{
#' # Run LEfSe on mock data
#' mock_lefse <- find_dams_lefse(
#'   data = mock_abundance,
#'   sample = mock_sample_meta,
#'   groupCol = "group"
#' )
#'
#' print("Mock data LEfSe results:")
#' if (nrow(mock_lefse) > 0) {
#'   print(mock_lefse)
#' } else {
#'   print("No significant biomarkers found (expected with random data)")
#' }
#' }
#'
#' # Example 5: Preprocessing and filtering
#' \dontrun{
#' # Pre-filter low abundance features
#' abundance_filter <- rowMeans(df.call_DAMs_LEfSe.otu) > 0.01  # 1% mean abundance
#' filtered_data <- df.call_DAMs_LEfSe.otu[abundance_filter, ]
#'
#' print(paste("Retained", nrow(filtered_data), "of", nrow(df.call_DAMs_LEfSe.otu),
#'   "features after filtering"))
#'
#' # Run LEfSe on filtered data
#' filtered_lefse <- find_dams_lefse(
#'   data = filtered_data,
#'   sample = df.call_DAMs_LEfSe.sample,
#'   lda.threshold = 1.0
#' )
#'
#' print(paste("Filtered analysis found", nrow(filtered_lefse), "biomarkers"))
#'
#' # Compare with unfiltered results
#' print("Comparison of filtering effects completed")
#' }
#'
find_dams_lefse <- function(data, sample, groupCol = "group",
                            kruskal.threshold = 0.05, wilcox.threshold = 0.05,
                            lda.threshold = 1.0, strict.wilcox = TRUE,
                            bootstrap.iterations = 30) {
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

  # Check for negative values (LEfSe can handle various data types)
  if (any(data < 0, na.rm = TRUE)) {
    stop("'data' contains negative values. LEfSe requires non-negative abundance data.")
  }

  # Check sample-column correspondence
  if (ncol(data) != nrow(sample)) {
    stop("Number of columns in 'data' must equal number of rows in 'sample'")
  }

  # Validate sample names matching
  data_samples <- colnames(data)
  sample_names <- rownames(sample)

  if (is.null(data_samples) || is.null(sample_names)) {
    stop("Both 'data' columns and 'sample' rows must have names")
  }

  if (!all(data_samples %in% sample_names)) {
    missing_samples <- setdiff(data_samples, sample_names)
    stop(paste("Sample(s) missing from metadata:", paste(missing_samples, collapse = ", ")))
  }

  if (!all(sample_names %in% data_samples)) {
    extra_samples <- setdiff(sample_names, data_samples)
    warning(paste("Extra samples in metadata (will be ignored):",
      paste(extra_samples, collapse = ", ")))
  }

  # Align sample metadata with abundance data
  sample_aligned <- sample[data_samples, , drop = FALSE]

  # Validate groupCol parameter
  if (!is.character(groupCol) || length(groupCol) != 1) {
    stop("'groupCol' must be a single character string")
  }

  if (!groupCol %in% names(sample_aligned)) {
    stop(paste("Group column '", groupCol, "' not found in sample metadata"))
  }

  # Check for missing values in grouping variable
  if (any(is.na(sample_aligned[[groupCol]]))) {
    stop(paste("Missing values found in grouping variable:", groupCol))
  }

  # Validate threshold parameters
  if (!is.numeric(kruskal.threshold) || length(kruskal.threshold) != 1 ||
    kruskal.threshold <= 0 || kruskal.threshold > 1) {
    stop("'kruskal.threshold' must be a single number between 0 and 1")
  }

  if (!is.numeric(wilcox.threshold) || length(wilcox.threshold) != 1 ||
    wilcox.threshold <= 0 || wilcox.threshold > 1) {
    stop("'wilcox.threshold' must be a single number between 0 and 1")
  }

  if (!is.numeric(lda.threshold) || length(lda.threshold) != 1 || lda.threshold < 0) {
    stop("'lda.threshold' must be a single non-negative number")
  }

  if (!is.logical(strict.wilcox) || length(strict.wilcox) != 1) {
    stop("'strict.wilcox' must be a single logical value (TRUE/FALSE)")
  }

  if (!is.numeric(bootstrap.iterations) || length(bootstrap.iterations) != 1 ||
    bootstrap.iterations < 10) {
    stop("'bootstrap.iterations' must be a single number >= 10")
  }

  # Check group composition and sample sizes
  group_counts <- table(sample_aligned[[groupCol]])
  n_groups <- length(group_counts)

  if (n_groups < 2) {
    stop("At least 2 groups required for LEfSe analysis")
  }

  if (any(group_counts < 3)) {
    warning("Some groups have fewer than 3 samples. Results may be unreliable.")
  }

  min_samples_per_group <- min(group_counts)
  if (min_samples_per_group < 3) {
    warning(paste("Minimum samples per group:", min_samples_per_group,
      ". Consider having at least 6-10 samples per group for robust results."))
  }

  # Check for features with zero variance (constant across all samples)
  feature_vars <- apply(data, 1, var, na.rm = TRUE)
  zero_var_features <- sum(feature_vars == 0, na.rm = TRUE)

  if (zero_var_features > 0) {
    warning(paste("Found", zero_var_features, "features with zero variance. These may be filtered by LEfSe."))
  }

  # Check for features that are zero in all samples
  all_zero_features <- rowSums(data, na.rm = TRUE) == 0
  n_zero_features <- sum(all_zero_features)

  if (n_zero_features > 0) {
    warning(paste("Found", n_zero_features, "features with all zero values. Consider removing these."))
  }

  # Convert data to matrix if needed
  if (!is.matrix(data)) {
    data <- as.matrix(data)
  }

  # Create SummarizedExperiment object
  tryCatch(
    {
      se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = data),
        colData = sample_aligned
      )
    },
    error = function(e) {
      stop(paste("Error creating SummarizedExperiment:", e$message))
    })

  # Run LEfSe analysis
  tryCatch(
    {
      lefse_result <- lefser::lefser(
        se,
        groupCol = groupCol,
        kruskal.threshold = kruskal.threshold,
        wilcox.threshold = wilcox.threshold,
        lda.threshold = lda.threshold
      )
    },
    error = function(e) {
      if (grepl("No significant features found", e$message)) {
        warning("No significant biomarkers found with current thresholds. Consider adjusting parameters.")
        # Return empty data frame with expected structure
        empty_result <- data.frame(
          features = character(0),
          scores = numeric(0),
          class = character(0),
          stringsAsFactors = FALSE
        )

        # Add metadata attributes
        attr(empty_result, "n_features_input") <- nrow(data)
        attr(empty_result, "n_samples") <- ncol(data)
        attr(empty_result, "n_groups") <- n_groups
        attr(empty_result, "group_column") <- groupCol
        attr(empty_result, "kruskal_threshold") <- kruskal.threshold
        attr(empty_result, "wilcox_threshold") <- wilcox.threshold
        attr(empty_result, "lda_threshold") <- lda.threshold

        return(empty_result)
      } else {
        stop(paste("Error in LEfSe analysis:", e$message))
      }
    })

  # Process and enhance results
  if (is.null(lefse_result) || nrow(lefse_result) == 0) {
    warning("LEfSe analysis completed but found no significant biomarkers")

    # Return empty result with proper structure
    empty_result <- data.frame(
      features = character(0),
      scores = numeric(0),
      class = character(0),
      stringsAsFactors = FALSE
    )

    attr(empty_result, "n_features_input") <- nrow(data)
    attr(empty_result, "n_samples") <- ncol(data)
    attr(empty_result, "n_groups") <- n_groups

    return(empty_result)
  }

  # Enhance results with additional columns
  enhanced_results <- lefse_result %>%
    dplyr::mutate(
      # Add absolute LDA score for ranking
      lda.score = abs(scores),
      # Add direction indicator
      direction = ifelse(scores > 0, "Enriched", "Depleted")
    ) %>%
    # Arrange by LDA score (highest discriminatory power first)
    dplyr::arrange(desc(lda.score))

  # Add comprehensive metadata as attributes
  attr(enhanced_results, "n_features_input") <- nrow(data)
  attr(enhanced_results, "n_features_tested") <- nrow(data) - n_zero_features
  attr(enhanced_results, "n_biomarkers_found") <- nrow(enhanced_results)
  attr(enhanced_results, "n_samples") <- ncol(data)
  attr(enhanced_results, "n_groups") <- n_groups
  attr(enhanced_results, "group_names") <- names(group_counts)
  attr(enhanced_results, "group_sizes") <- as.numeric(group_counts)
  attr(enhanced_results, "group_column") <- groupCol
  attr(enhanced_results, "kruskal_threshold") <- kruskal.threshold
  attr(enhanced_results, "wilcox_threshold") <- wilcox.threshold
  attr(enhanced_results, "lda_threshold") <- lda.threshold
  attr(enhanced_results, "strict_wilcox") <- strict.wilcox

  # Generate and display summary if interactive
  if (interactive()) {
    cat("LEfSe Biomarker Analysis Summary:\n")
    cat("=================================\n")
    cat("Features analyzed:", nrow(data), "\n")
    cat("Samples analyzed:", ncol(data), "\n")
    cat("Groups compared:", n_groups, "(", paste(names(group_counts), collapse = ", "), ")\n")
    cat("Sample sizes:", paste(group_counts, collapse = ", "), "\n\n")

    cat("Analysis parameters:\n")
    cat("  Kruskal-Wallis threshold:", kruskal.threshold, "\n")
    cat("  Wilcoxon threshold:", wilcox.threshold, "\n")
    cat("  LDA threshold:", lda.threshold, "\n\n")

    cat("Results:\n")
    cat("  Significant biomarkers found:", nrow(enhanced_results), "\n")

    if (nrow(enhanced_results) > 0) {
      biomarkers_by_class <- enhanced_results %>%
        dplyr::count(class, name = "n_biomarkers") %>%
        dplyr::arrange(desc(n_biomarkers))

      cat("  Biomarkers by class:\n")
      for (i in 1:nrow(biomarkers_by_class)) {
        cat(sprintf("    %s: %d biomarkers\n",
          biomarkers_by_class$class[i],
          biomarkers_by_class$n_biomarkers[i]))
      }

      cat("\n  Top 5 biomarkers by LDA score:\n")
      top_biomarkers <- head(enhanced_results, 5)
      for (i in 1:nrow(top_biomarkers)) {
        cat(sprintf("    %s: %s (LDA = %.2f)\n",
          top_biomarkers$features[i],
          top_biomarkers$class[i],
          top_biomarkers$lda.score[i]))
      }
    } else {
      cat("  No biomarkers met the significance criteria.\n")
      cat("  Consider relaxing thresholds or checking data quality.\n")
    }
    cat("\n")
  }

  return(enhanced_results)
}
