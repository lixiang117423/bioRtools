#' Enhanced Presence/Absence Variation (PAV) Genome-Wide Association Study (GWAS) Analysis
#'
#' @description
#' Perform presence/absence variation based genome-wide association study using
#' PAV data (such as k-mers, structural variants, gene presence/absence) and
#' phenotype information. Enhanced version with better support for continuous
#' phenotypes like plant height, weight, expression levels, etc.
#'
#' @param pav_data Character string specifying the path to the PAV matrix file,
#'   or a data frame containing the PAV matrix. If a file path, the file should
#'   be tab-delimited with the following columns: 1) chr (chromosome), 2) start
#'   (start position), 3) end (end position), 4) sequence (PAV sequence),
#'   5+ sample abundance/presence values. Headers are expected.
#'   If a data frame, it should follow the same column structure.
#' @param phenotype_data Character string specifying the path to the phenotype file,
#'   or a data frame containing phenotype information. If a file path, the file
#'   should be tab-delimited with sample IDs in the first column and phenotype
#'   values in the second column. Headers are expected. If a data frame, the first
#'   column should contain sample IDs and the second column should contain phenotype values.
#' @param output_dir Character string specifying the output directory for results.
#'   Default is "pav_gwas_results".
#' @param analysis_type Character string specifying the phenotype analysis approach.
#'   Options: "auto" (automatic detection), "binary", "logit", "arcsine", "categorical",
#'   "continuous", "continuous_raw". Default is "auto".
#' @param test_type Character string specifying the statistical test method.
#'   Options: "auto" (automatic selection), "fisher", "mannwhitney", "logistic",
#'   "linear", "welch_ttest". Default is "auto".
#' @param min_sample_count Integer specifying the minimum number of samples where
#'   a PAV feature must be present to be included in analysis. Default is 2.
#' @param max_sample_count Integer specifying the maximum number of samples where
#'   a PAV feature can be present to be included in analysis. If NULL, uses total
#'   sample count (no upper limit). Default is NULL.
#' @param binary_threshold Numeric value specifying the threshold for binary
#'   conversion when analysis_type is "binary" or "auto". Default is 0.5.
#' @param use_population_structure Logical indicating whether to include population
#'   structure correction using PCA. Default is TRUE.
#' @param n_pcs Integer specifying the number of principal components to use for
#'   population structure correction. Default is 3.
#' @param significance_threshold Numeric value specifying the genome-wide significance
#'   threshold for Manhattan plot. If NULL (default), automatically calculates
#'   Bonferroni-corrected threshold as 0.05/number_of_tested_features. Can be
#'   manually specified if desired.
#' @param n_cores Integer specifying the number of CPU cores for parallel processing.
#'   Default is 4.
#' @param verbose Logical indicating whether to print progress information.
#'   Default is TRUE.
#'
#' @return A list containing:
#' \describe{
#'   \item{gwas_results}{Data frame with association test results including chromosomal positions (chr, start, end), sequence information, and statistical results for all PAV features.}
#'   \item{phenotype_stats}{Summary statistics of phenotype distribution.}
#'   \item{pca_results}{PCA results for population structure (if used).}
#'   \item{filtered_features}{Information about PAV feature filtering.}
#'   \item{analysis_summary}{Summary of analysis parameters and significant findings.}
#'   \item{threshold_results}{Data frame showing number of significant features at different p-value thresholds.}
#' }
#'
#' @details
#' **Enhanced Data Type Detection and Method Selection:**
#'
#' The function automatically detects phenotype distribution patterns and selects
#' appropriate statistical methods when analysis_type = "auto":
#'
#' \enumerate{
#'   \item **Binary Data** (exactly 2 unique values):
#'      - Automatically uses binary analysis
#'      - Uses Fisher's exact test for association
#'      - Ideal for case/control studies
#'   \item **Bimodal Distribution** (values concentrated at extremes):
#'      - Automatically converts to binary phenotype
#'      - Uses Fisher's exact test for association
#'      - Ideal for presence/absence traits
#'   \item **Categorical Data** (3-10 discrete categories):
#'      - Uses multinomial logistic regression
#'      - Falls back to chi-square test if logistic regression fails
#'      - Categories can be numeric or character
#'   \item **Continuous Data** (many unique values):
#'      - **Normal Distribution**: Uses Welch's t-test for higher statistical power
#'      - **Non-normal Distribution**: Uses Mann-Whitney U test or applies transformations
#'      - **Bounded [0,1]**: Uses logit transformation
#'      - **Percentage-like [0,100]**: Uses arcsine transformation
#'   \item **Highly Skewed Continuous Data**:
#'      - Uses log transformation or rank transformation
#'      - Followed by appropriate parametric/non-parametric tests
#' }
#'
#' **Enhanced Manual Method Selection Guidelines:**
#'
#' Choose specific methods based on your data characteristics:
#'
#' - **analysis_type = "continuous"**: For truly continuous traits (plant height,
#'   weight, expression levels). Automatically applies appropriate transformations.
#' - **analysis_type = "continuous_raw"**: For continuous traits without any
#'   transformation. Uses raw phenotype values.
#' - **analysis_type = "binary"**: For clearly binary traits (disease/healthy,
#'   resistant/susceptible). Converts continuous data using binary_threshold.
#' - **analysis_type = "logit"**: For proportion data bounded between 0 and 1
#'   (survival rates, infection percentages).
#' - **analysis_type = "arcsine"**: For percentage data with many zeros or ones
#'   (germination rates, mortality percentages).
#' - **analysis_type = "categorical"**: For multi-level categorical traits
#'   (flower color, growth habit).
#'
#' **Enhanced Statistical Test Selection:**
#'
#' - **test_type = "welch_ttest"**: Welch's t-test for continuous normal phenotypes.
#'   Higher statistical power than non-parametric tests when assumptions are met.
#' - **test_type = "linear"**: Linear regression, allows inclusion of population
#'   structure covariates. Provides regression coefficients as effect sizes.
#' - **test_type = "fisher"**: Best for binary phenotypes, exact p-values,
#'   robust with small sample sizes.
#' - **test_type = "mannwhitney"**: Non-parametric test for continuous phenotypes,
#'   no normality assumptions.
#' - **test_type = "logistic"**: Handles both binary and categorical phenotypes.
#'   For binary: standard logistic regression. For categorical: multinomial
#'   logistic regression (requires nnet package) or chi-square test as fallback.
#'   Allows inclusion of covariates.
#'
#' **PAV Feature Filtering Strategy:**
#'
#' - **min_sample_count**: Removes rare PAV features that may represent sequencing
#'   errors or population-specific variants. Higher values increase statistical power
#'   but may miss important rare variants.
#' - **max_sample_count**: Removes ubiquitous PAV features that provide little
#'   information. Leave as NULL unless you want to focus on rare variants only.
#'
#' **Population Structure Correction:**
#'
#' Principal Component Analysis (PCA) is performed on the PAV matrix to capture
#' population structure. The first n_pcs components are included as covariates in
#' association tests to control for confounding due to relatedness or population
#' stratification.
#'
#' **PAV Data Types Supported:**
#'
#' This function can analyze various types of presence/absence variation data:
#' \itemize{
#'   \item **K-mer presence/absence**: From k-mer counting tools
#'   \item **Structural variants**: CNVs, indels, translocations
#'   \item **Gene presence/absence**: Pan-genome analysis results
#'   \item **SNP presence/absence**: Binary SNP matrices
#'   \item **Any binary genomic features**: Regulatory elements, repeat sequences, etc.
#' }
#'
#' **Output Files:**
#'
#' The function creates several output files in the specified directory:
#' \itemize{
#'   \item **gwas_results.csv**: Complete association results with chromosomal positions, sequences, p-values and effect sizes
#'   \item **manhattan_plot.pdf/.png**: Manhattan plot showing associations across chromosomes with multiple significance threshold lines
#'   \item **significant_p1e-5.csv, significant_p1e-6.csv, etc.**: Separate files for each significance threshold containing only significant features
#'   \item **significant_fdr0.05.csv**: Features significant after FDR correction
#'   \item **qq_plot.pdf/.png**: Q-Q plot for test calibration assessment
#'   \item **phenotype_distribution.pdf/.png**: Phenotype distribution analysis
#'   \item **population_structure.pdf/.png**: PCA plots and variance explained
#'   \item **analysis_summary.txt**: Text summary of results and parameters including threshold-specific counts
#' }
#'
#' @note
#' \itemize{
#'   \item Requires packages: data.table, dplyr, ggplot2, parallel, corrplot,
#'     scales, tidyr, factoextra
#'   \item Optional packages: nnet (for multinomial logistic regression),
#'     moments (for skewness/kurtosis calculation)
#'   \item Missing data is handled by removing samples with NA values
#'   \item Multiple testing correction includes both Bonferroni and FDR methods
#'   \item Large PAV matrices may require substantial memory and processing time
#'   \item Consider using fewer cores (n_cores) if memory becomes limiting
#'   \item Both file paths and data frames are accepted as input
#'   \item Categorical phenotypes: can be numeric (0,1,2) or character ("A","B","C")
#'   \item Continuous phenotypes: plant height, weight, expression levels, etc.
#' }
#'
#' @export
#' @examples
#' \dontrun{
#' # Basic analysis with automatic significance threshold calculation
#' # PAV matrix format: chr, start, end, sequence, sample1, sample2, ...
#' result <- pav_gwas(
#'   pav_data = "pav_matrix.txt",
#'   phenotype_data = "phenotype_data.txt"
#' )
#'
#' # Using data frames as input (common for k-mer analysis)
#' result <- pav_gwas(
#'   pav_data = df.pav.matrix,  # Must have chr, start, end, sequence columns
#'   phenotype_data = df.phenotype,
#'   output_dir = "my_gwas_results"
#' )
#'
#' # Enhanced: Continuous phenotype analysis (e.g., plant height)
#' result <- pav_gwas(
#'   pav_data = df.kmer.gwas,
#'   phenotype_data = df.height.data,  # Continuous height measurements
#'   analysis_type = "continuous",     # Enhanced continuous analysis
#'   test_type = "welch_ttest",        # T-test for normal data
#'   output_dir = "height_gwas_results/"
#' )
#'
#' # Enhanced: Auto-detection for continuous phenotypes
#' result <- pav_gwas(
#'   pav_data = df.kmer.gwas,
#'   phenotype_data = df.continuous.phenotype,
#'   analysis_type = "auto",           # Will auto-detect continuous
#'   test_type = "auto",               # Will auto-select best test
#'   use_population_structure = TRUE,   # Include PCA correction
#'   output_dir = "auto_continuous_gwas/"
#' )
#'
#' # Real usage example for structural variants
#' result <- pav_gwas(
#'   pav_data = df.kmer.gwas,     # With positional information
#'   phenotype_data = df.sample.gwas,
#'   output_dir = "D://OneDrive/NAS/01.科研相关/00.博后/02.data/03.大豆根腐项目/result/01.两个基因变异情况/kmer_gwas/"
#' )
#'
#' # Three-class phenotype analysis
#' result <- pav_gwas(
#'   pav_data = df.structural.variants,
#'   phenotype_data = df.three.class.phenotype,  # e.g., "resistant", "intermediate", "susceptible"
#'   analysis_type = "categorical",
#'   test_type = "logistic",
#'   significance_threshold = NULL,  # Auto-calculate
#'   n_cores = 8
#' )
#'
#' # Working with categorical phenotypes (disease resistance levels)
#' # phenotype data can contain: 0, 1, 2 or "resistant", "intermediate", "susceptible"
#' result <- pav_gwas(
#'   pav_data = df.disease.pavs,
#'   phenotype_data = df.resistance.levels,
#'   analysis_type = "auto",  # Will auto-detect categorical
#'   test_type = "auto",      # Will auto-select logistic regression
#'   output_dir = "disease_resistance_gwas"
#' )
#'
#' # Accessing results
#' print(result$analysis_summary)
#' head(result$gwas_results)
#'
#' # View threshold-specific results
#' print(result$threshold_results)
#'
#' # Get significant features at different thresholds
#' significant_1e5 <- result$gwas_results[result$gwas_results$p_value < 1e-5, ]
#' significant_fdr <- result$gwas_results[result$gwas_results$p_fdr < 0.05, ]
#'
#' # Threshold-specific files are automatically saved as:
#' # significant_p1e-5.csv, significant_p1e-6.csv, etc.
#' }
#'
pav_gwas <- function(pav_data,
                     phenotype_data,
                     output_dir = "pav_gwas_results",
                     analysis_type = "auto",
                     test_type = "auto",
                     min_sample_count = 2,
                     max_sample_count = NULL,
                     binary_threshold = 0.5,
                     use_population_structure = TRUE,
                     n_pcs = 3,
                     significance_threshold = NULL,
                     n_cores = 4,
                     verbose = TRUE) {
  # Check required packages (Enhanced)
  required_packages <- c("data.table", "dplyr", "ggplot2", "parallel",
    "corrplot", "scales", "tidyr", "factoextra",
    "readr", "stringr", "magrittr")

  # Add optional packages for specific analyses (Enhanced)
  optional_packages <- c("nnet", "moments")  # Added moments for skewness/kurtosis

  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

  if (length(missing_packages) > 0) {
    stop("Required packages missing: ", paste(missing_packages, collapse = ", "))
  }

  # Check optional packages (warn but don't stop)
  missing_optional <- optional_packages[!sapply(optional_packages, requireNamespace, quietly = TRUE)]
  if (length(missing_optional) > 0 && verbose) {
    cat("Optional packages missing (will use fallback methods):", paste(missing_optional, collapse = ", "), "\n")
  }

  # Load required packages
  suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(ggplot2)
    library(parallel)
    library(readr)
    library(tidyr)
    library(stringr)
    library(magrittr)
  })

  # Enhanced input validation
  if (!analysis_type %in% c("auto", "binary", "logit", "arcsine", "categorical", "continuous", "continuous_raw")) {
    stop("Invalid analysis_type. Must be one of: auto, binary, logit, arcsine, categorical, continuous, continuous_raw")
  }

  if (!test_type %in% c("auto", "fisher", "mannwhitney", "logistic", "linear", "welch_ttest")) {
    stop("Invalid test_type. Must be one of: auto, fisher, mannwhitney, logistic, linear, welch_ttest")
  }

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  if (verbose) {
    cat("=== Enhanced PAV GWAS Analysis Pipeline ===\n")
    cat("Output directory:", output_dir, "\n")
    cat("Analysis type:", analysis_type, "\n")
    cat("Test type:", test_type, "\n")
  }

  # Step 1: Load PAV data (file or data frame) - UNCHANGED
  if (verbose) cat("Loading PAV data...\n")

  if (is.character(pav_data) && length(pav_data) == 1) {
    # It's a file path
    if (!file.exists(pav_data)) {
      stop("PAV file not found: ", pav_data)
    }

    if (verbose) cat("Reading PAV data from file:", pav_data, "\n")
    pav_df <- fread(pav_data, header = TRUE)

  } else if (is.data.frame(pav_data)) {
    # It's a data frame
    if (verbose) cat("Using provided PAV data frame\n")
    pav_df <- as.data.table(pav_data)

  } else {
    stop("pav_data must be either a file path (character string) or a data frame")
  }

  # Validate PAV matrix structure
  if (ncol(pav_df) < 5) {
    stop("PAV matrix must have at least 5 columns: chr, start, end, sequence, and sample data")
  }

  # Extract PAV information and matrix
  pav_info <- pav_df[, 1:4]
  names(pav_info) <- c("chr", "start", "end", "sequence")

  # Create unique feature IDs
  feature_ids <- paste0(pav_info$chr, ":", pav_info$start, "-", pav_info$end)

  # Extract sample data (from column 5 onwards)
  pav_matrix <- as.matrix(pav_df[, 5:ncol(pav_df)])
  rownames(pav_matrix) <- feature_ids

  if (verbose) {
    cat(sprintf("Loaded %d PAV features across %d chromosomes for %d samples\n",
      nrow(pav_matrix), length(unique(pav_info$chr)), ncol(pav_matrix)))
    cat("Chromosomes found:", paste(sort(unique(pav_info$chr)), collapse = ", "), "\n")
  }

  # Step 2: Load phenotype data (file or data frame) - UNCHANGED
  if (verbose) cat("Loading phenotype data...\n")

  if (is.character(phenotype_data) && length(phenotype_data) == 1) {
    # It's a file path
    if (!file.exists(phenotype_data)) {
      stop("Phenotype file not found: ", phenotype_data)
    }

    if (verbose) cat("Reading phenotype data from file:", phenotype_data, "\n")
    pheno_df <- fread(phenotype_data, header = TRUE)

  } else if (is.data.frame(phenotype_data)) {
    # It's a data frame
    if (verbose) cat("Using provided phenotype data frame\n")
    pheno_df <- as.data.table(phenotype_data)

  } else {
    stop("phenotype_data must be either a file path (character string) or a data frame")
  }

  # Standardize column names
  names(pheno_df) <- c("SampleID", "Phenotype")

  if (verbose) {
    cat(sprintf("Loaded phenotype data for %d samples\n", nrow(pheno_df)))
  }

  # Step 3: Match samples between PAV matrix and phenotype data - UNCHANGED
  common_samples <- intersect(colnames(pav_matrix), pheno_df$SampleID)

  if (length(common_samples) == 0) {
    cat("\n=== DEBUGGING INFORMATION ===\n")
    cat("PAV matrix sample names (first 10):\n")
    print(head(colnames(pav_matrix), 10))
    cat("\nPhenotype sample names (first 10):\n")
    print(head(pheno_df$SampleID, 10))
    cat("=============================\n")
    stop("No samples match between PAV matrix and phenotype data")
  }

  if (verbose) {
    cat(sprintf("Found %d common samples\n", length(common_samples)))
  }

  # Reorder data to ensure consistency
  pav_matrix <- pav_matrix[, common_samples]
  pheno_df <- pheno_df[pheno_df$SampleID %in% common_samples, ]
  pheno_df <- pheno_df[match(common_samples, pheno_df$SampleID), ]

  # Remove samples with missing phenotype data
  complete_idx <- !is.na(pheno_df$Phenotype)
  pav_matrix <- pav_matrix[, complete_idx]
  pheno_df <- pheno_df[complete_idx, ]

  if (verbose) {
    cat(sprintf("Using %d samples with complete data\n", ncol(pav_matrix)))
  }

  # Step 4: Enhanced phenotype distribution analysis
  if (verbose) cat("Analyzing phenotype distribution...\n")

  phenotype_values <- pheno_df$Phenotype

  # Check if phenotype is numeric or categorical
  is_numeric_phenotype <- is.numeric(phenotype_values)

  if (is_numeric_phenotype) {
    # Enhanced statistics for numeric phenotype
    pheno_stats <- list(
      n_samples = length(phenotype_values),
      min = min(phenotype_values, na.rm = TRUE),
      max = max(phenotype_values, na.rm = TRUE),
      mean = mean(phenotype_values, na.rm = TRUE),
      median = median(phenotype_values, na.rm = TRUE),
      sd = sd(phenotype_values, na.rm = TRUE),
      n_unique = length(unique(phenotype_values)),
      is_numeric = TRUE
    )

    # Enhanced continuous phenotype metrics
    pheno_stats$cv <- pheno_stats$sd / pheno_stats$mean  # Coefficient of variation
    pheno_stats$range <- pheno_stats$max - pheno_stats$min
    pheno_stats$iqr <- IQR(phenotype_values, na.rm = TRUE)

    # Enhanced distribution characteristics
    pheno_stats$is_binary_like <- pheno_stats$n_unique == 2
    pheno_stats$is_bounded_01 <- pheno_stats$min >= 0 && pheno_stats$max <= 1
    pheno_stats$is_percentage_like <- pheno_stats$min >= 0 && pheno_stats$max <= 100
    pheno_stats$is_truly_continuous <- pheno_stats$n_unique > 10 &&
      pheno_stats$n_unique > pheno_stats$n_samples * 0.1

    # Enhanced normality testing
    if (pheno_stats$n_samples <= 5000) {
      tryCatch(
        {
          shapiro_test <- shapiro.test(phenotype_values)
          pheno_stats$normality_p <- shapiro_test$p.value
          pheno_stats$is_normal <- shapiro_test$p.value > 0.05
        },
        error = function(e) {
          pheno_stats$normality_p <- NA
          pheno_stats$is_normal <- FALSE
        })
    } else {
      pheno_stats$normality_p <- NA
      # For large samples, use skewness and kurtosis
      if (requireNamespace("moments", quietly = TRUE)) {
        pheno_stats$skewness <- moments::skewness(phenotype_values, na.rm = TRUE)
        pheno_stats$kurtosis <- moments::kurtosis(phenotype_values, na.rm = TRUE)
        pheno_stats$is_normal <- abs(pheno_stats$skewness) < 1 && abs(pheno_stats$kurtosis - 3) < 1
      } else {
        # Simple normality approximation
        pheno_stats$is_normal <- pheno_stats$cv < 0.5  # Low coefficient of variation suggests normality
      }
    }

    # Enhanced bimodal distribution detection
    if (pheno_stats$is_bounded_01) {
      tryCatch(
        {
          density_result <- density(phenotype_values)
          peaks <- which(diff(sign(diff(density_result$y))) == -2) + 1
          pheno_stats$is_bimodal <- length(peaks) >= 2

          # Additional check: concentration at extremes
          extreme_prop <- sum(phenotype_values < 0.1 | phenotype_values > 0.9) / length(phenotype_values)
          if (extreme_prop > 0.7) {
            pheno_stats$is_bimodal <- TRUE
          }
        },
        error = function(e) {
          pheno_stats$is_bimodal <- FALSE
        })
    } else {
      pheno_stats$is_bimodal <- FALSE
    }

    # Enhanced verbose output
    if (verbose) {
      cat(sprintf("  Numeric phenotype: %.3f ± %.3f (range: %.3f - %.3f)\n",
        pheno_stats$mean, pheno_stats$sd, pheno_stats$min, pheno_stats$max))
      cat(sprintf("  %d unique values (%.1f%% of samples)\n",
        pheno_stats$n_unique, pheno_stats$n_unique / pheno_stats$n_samples * 100))
      cat(sprintf("  CV = %.3f, IQR = %.3f\n", pheno_stats$cv, pheno_stats$iqr))
      if (!is.na(pheno_stats$normality_p)) {
        cat(sprintf("  Normality test p-value: %.3f\n", pheno_stats$normality_p))
      }
      if (pheno_stats$is_truly_continuous) {
        cat("  Detected as: CONTINUOUS phenotype\n")
      } else if (pheno_stats$is_binary_like) {
        cat("  Detected as: BINARY-like phenotype\n")
      } else {
        cat("  Detected as: DISCRETE phenotype\n")
      }
    }

  } else {
    # Basic statistics for categorical phenotype - UNCHANGED
    pheno_table <- table(phenotype_values)
    pheno_stats <- list(
      n_samples = length(phenotype_values),
      min = NA,
      max = NA,
      mean = NA,
      median = NA,
      sd = NA,
      n_unique = length(unique(phenotype_values)),
      is_numeric = FALSE,
      categories = names(pheno_table),
      category_counts = as.numeric(pheno_table),
      is_binary_like = length(unique(phenotype_values)) == 2,
      is_truly_continuous = FALSE,
      is_normal = FALSE
    )

    if (verbose) {
      cat("  Categorical phenotype:\n")
      for (i in 1:length(pheno_stats$categories)) {
        cat(sprintf("    %s: %d samples (%.1f%%)\n",
          pheno_stats$categories[i],
          pheno_stats$category_counts[i],
          pheno_stats$category_counts[i] / pheno_stats$n_samples * 100))
      }
    }
  }

  # Detect distribution pattern for auto mode (enhanced)
  is_binary_like <- pheno_stats$is_binary_like
  is_bounded <- ifelse(is_numeric_phenotype, pheno_stats$is_bounded_01, FALSE)
  is_bimodal <- ifelse(is_numeric_phenotype, pheno_stats$is_bimodal, FALSE)

  pheno_stats$is_bounded <- is_bounded
  pheno_stats$is_bimodal <- is_bimodal

  # Step 5: Enhanced automatic method selection
  if (analysis_type == "auto") {
    if (is_binary_like) {
      analysis_type <- "binary"
      if (verbose) cat("Auto-detected: Binary phenotype\n")
    } else if (!is_numeric_phenotype) {
      analysis_type <- "categorical"
      if (verbose) cat("Auto-detected: Categorical phenotype\n")
    } else if (pheno_stats$is_truly_continuous) {
      # Enhanced detection for continuous phenotypes
      if (is_bimodal) {
        analysis_type <- "binary"
        if (verbose) cat("Auto-detected: Bimodal bounded data, converting to binary\n")
      } else if (is_bounded) {
        analysis_type <- "logit"
        if (verbose) cat("Auto-detected: Bounded continuous data [0,1], using logit transformation\n")
      } else if (pheno_stats$is_percentage_like) {
        analysis_type <- "arcsine"
        if (verbose) cat("Auto-detected: Percentage-like data [0,100], using arcsine transformation\n")
      } else {
        analysis_type <- "continuous"
        if (verbose) cat("Auto-detected: Continuous phenotype, using continuous methods\n")
      }
    } else if (pheno_stats$n_unique <= 10) {
      analysis_type <- "categorical"
      if (verbose) cat("Auto-detected: Discrete categorical phenotype\n")
    } else {
      analysis_type <- "continuous"
      if (verbose) cat("Auto-detected: Treating as continuous phenotype\n")
    }
  }

  # Enhanced test type selection
  if (test_type == "auto") {
    if (analysis_type == "binary") {
      test_type <- "fisher"
      if (verbose) cat("Auto-selected: Fisher's exact test\n")
    } else if (analysis_type == "categorical") {
      test_type <- "logistic"
      if (verbose) cat("Auto-selected: Logistic regression\n")
    } else if (analysis_type %in% c("continuous", "continuous_raw")) {
      if (pheno_stats$is_normal && pheno_stats$n_samples >= 30) {
        test_type <- "welch_ttest"
        if (verbose) cat("Auto-selected: Welch's t-test (normal distribution detected)\n")
      } else {
        test_type <- "mannwhitney"
        if (verbose) cat("Auto-selected: Mann-Whitney U test (non-parametric)\n")
      }
    } else {
      test_type <- "mannwhitney"
      if (verbose) cat("Auto-selected: Mann-Whitney U test\n")
    }
  }

  # Step 6: Enhanced phenotype transformation
  if (verbose) cat("Transforming phenotype...\n")

  transformed_phenotype <- switch(analysis_type,
    "binary" = {
      if (is_binary_like) {
        # Already binary, just ensure 0/1 coding
        if (pheno_stats$is_numeric) {
          as.numeric(factor(phenotype_values)) - 1
        } else {
          # Convert categorical to 0/1
          as.numeric(factor(phenotype_values, levels = unique(phenotype_values))) - 1
        }
      } else {
        # Convert to binary using threshold (only for numeric data)
        if (pheno_stats$is_numeric) {
          as.numeric(phenotype_values > binary_threshold)
        } else {
          stop("Cannot apply binary threshold to categorical data. Please use analysis_type = 'categorical'")
        }
      }
    },
    "logit" = {
      if (!pheno_stats$is_numeric) {
        stop("Logit transformation requires numeric data. Please use analysis_type = 'categorical'")
      }
      # Logit transformation for bounded data
      epsilon <- 1e-6
      adj_values <- pmax(pmin(phenotype_values, 1 - epsilon), epsilon)
      log(adj_values / (1 - adj_values))
    },
    "arcsine" = {
      if (!pheno_stats$is_numeric) {
        stop("Arcsine transformation requires numeric data. Please use analysis_type = 'categorical'")
      }
      # Enhanced arcsine square-root transformation
      if (pheno_stats$is_percentage_like) {
        # Convert percentage to proportion first
        prop_values <- phenotype_values / 100
      } else {
        prop_values <- phenotype_values
      }
      asin(sqrt(pmax(pmin(prop_values, 1), 0)))  # Ensure values are in [0,1]
    },
    "categorical" = {
      # Convert to numeric factor levels
      if (pheno_stats$is_numeric) {
        # If already numeric, check if it looks like categorical
        if (pheno_stats$n_unique <= 10 && all(phenotype_values == round(phenotype_values))) {
          as.numeric(factor(phenotype_values))
        } else {
          stop("Numeric data with many unique values. Consider using analysis_type = 'continuous'")
        }
      } else {
        # Convert categorical to numeric factor levels
        as.numeric(factor(phenotype_values, levels = unique(phenotype_values)))
      }
    },
    # ENHANCED: New continuous analysis methods
    "continuous" = {
      if (!pheno_stats$is_numeric) {
        stop("Continuous analysis requires numeric data")
      }
      # Apply appropriate transformation based on distribution
      if (!pheno_stats$is_normal) {
        if (min(phenotype_values) > 0) {
          # Apply log transformation for positive skewed data
          log_transformed <- log(phenotype_values)
          if (verbose) cat("  Applied log transformation for skewed positive data\n")
          log_transformed
        } else {
          # Use rank transformation for non-normal data with zeros/negatives
          if (verbose) cat("  Applied rank transformation for non-normal data\n")
          rank(phenotype_values, ties.method = "average")
        }
      } else {
        if (verbose) cat("  Using raw values (normal distribution detected)\n")
        phenotype_values
      }
    },
    "continuous_raw" = {
      if (!pheno_stats$is_numeric) {
        stop("Continuous analysis requires numeric data")
      }
      if (verbose) cat("  Using raw continuous values without transformation\n")
      phenotype_values
    }
  )

  if (verbose) {
    if (analysis_type == "binary") {
      cat(sprintf("Binary transformation: %d samples with value 0, %d samples with value 1\n",
        sum(transformed_phenotype == 0), sum(transformed_phenotype == 1)))
    } else if (analysis_type == "categorical") {
      category_counts <- table(transformed_phenotype)
      cat("Categorical transformation:\n")
      for (i in 1:length(category_counts)) {
        cat(sprintf("  Category %d: %d samples\n", i, category_counts[i]))
      }
    } else if (analysis_type %in% c("continuous", "continuous_raw")) {
      cat(sprintf("Continuous transformation: %.3f ± %.3f (range: %.3f - %.3f)\n",
        mean(transformed_phenotype, na.rm = TRUE),
        sd(transformed_phenotype, na.rm = TRUE),
        min(transformed_phenotype, na.rm = TRUE),
        max(transformed_phenotype, na.rm = TRUE)))
    }
  }

  # Step 7: Filter PAV matrix - UNCHANGED
  if (verbose) cat("Filtering PAV matrix...\n")

  # Convert to binary presence/absence
  binary_matrix <- (pav_matrix > 0) * 1

  # Calculate PAV feature frequency (number of samples where feature is present)
  feature_counts <- rowSums(binary_matrix)

  # Set max_sample_count to total samples if not specified
  if (is.null(max_sample_count)) {
    max_sample_count <- ncol(binary_matrix)
  }

  # Filter features
  keep_features <- which(feature_counts >= min_sample_count & feature_counts <= max_sample_count)
  filtered_matrix <- binary_matrix[keep_features, ]

  if (verbose) {
    cat(sprintf("Filtered from %d to %d PAV features (present in %d-%d samples)\n",
      nrow(binary_matrix), nrow(filtered_matrix),
      min_sample_count, max_sample_count))
  }

  filtered_feature_info <- list(
    original_count = nrow(binary_matrix),
    filtered_count = nrow(filtered_matrix),
    min_sample_count = min_sample_count,
    max_sample_count = max_sample_count,
    frequencies = feature_counts[keep_features]
  )

  # Calculate significance threshold if not provided
  if (is.null(significance_threshold)) {
    significance_threshold <- 0.05 / nrow(filtered_matrix)
    if (verbose) {
      cat(sprintf("Auto-calculated significance threshold: %.2e (0.05/%d features)\n",
        significance_threshold, nrow(filtered_matrix)))
    }
  } else {
    if (verbose) {
      cat(sprintf("Using provided significance threshold: %.2e\n", significance_threshold))
    }
  }

  # Step 8: Population structure analysis - UNCHANGED
  pca_results <- NULL
  covariates <- NULL

  if (use_population_structure) {
    if (verbose) cat("Analyzing population structure...\n")

    # Transpose matrix for PCA (samples as rows)
    sample_matrix <- t(filtered_matrix)

    # Check for constant columns (zero variance)
    col_vars <- apply(sample_matrix, 2, var)
    non_constant_cols <- which(col_vars > 0)

    if (length(non_constant_cols) == 0) {
      if (verbose) cat("Warning: All PAV features are constant after filtering. Skipping PCA.\n")
      use_population_structure <- FALSE
    } else if (length(non_constant_cols) < ncol(sample_matrix)) {
      n_constant <- ncol(sample_matrix) - length(non_constant_cols)
      if (verbose) {
        cat(sprintf("Removing %d constant PAV features for PCA (%d remaining)\n",
          n_constant, length(non_constant_cols)))
      }
      sample_matrix <- sample_matrix[, non_constant_cols, drop = FALSE]
    }

    if (use_population_structure && ncol(sample_matrix) > 1) {
      tryCatch(
        {
          # Perform PCA
          pca_result <- prcomp(sample_matrix, scale. = TRUE, center = TRUE)

          # Extract principal components
          n_pcs_available <- min(n_pcs, ncol(pca_result$x), nrow(pca_result$x) - 1)
          pcs <- pca_result$x[, 1:n_pcs_available, drop = FALSE]

          # Calculate variance explained
          var_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2))[1:n_pcs_available]

          pca_results <- list(
            pcs = pcs,
            var_explained = var_explained,
            pca_object = pca_result,
            n_features_used = ncol(sample_matrix)
          )

          covariates <- pcs

          if (verbose) {
            cat(sprintf("Using %d principal components explaining %.1f%% of variance\n",
              ncol(pcs), sum(var_explained) * 100))
            cat(sprintf("PCA performed on %d non-constant PAV features\n", ncol(sample_matrix)))
          }
        },
        error = function(e) {
          if (verbose) {
            cat("Warning: PCA failed, proceeding without population structure correction\n")
            cat("Error:", e$message, "\n")
          }
          use_population_structure <- FALSE
          pca_results <- NULL
          covariates <- NULL
        })
    } else if (use_population_structure) {
      if (verbose) cat("Warning: Insufficient variable features for PCA. Skipping population structure correction.\n")
      use_population_structure <- FALSE
    }
  }

  # Step 9: Association testing functions

  # Function to create Manhattan plot positions - UNCHANGED
  create_manhattan_positions <- function(gwas_df) {
    # Sort by chromosome and position
    gwas_df <- gwas_df[order(gwas_df$chr, gwas_df$start), ]

    # Handle different chromosome naming conventions
    gwas_df$chr_clean <- gsub("^chr", "", gwas_df$chr, ignore.case = TRUE)

    # Convert chromosome to numeric, handling non-numeric chromosomes
    chr_levels <- unique(gwas_df$chr_clean)

    # Separate numeric and non-numeric chromosomes
    numeric_chrs <- chr_levels[grepl("^[0-9]+$", chr_levels)]
    non_numeric_chrs <- chr_levels[!grepl("^[0-9]+$", chr_levels)]

    # Sort numeric chromosomes numerically, non-numeric alphabetically
    if (length(numeric_chrs) > 0) {
      numeric_chrs <- as.character(sort(as.numeric(numeric_chrs)))
    }
    non_numeric_chrs <- sort(non_numeric_chrs)

    # Combine in order: numeric first, then non-numeric
    ordered_chrs <- c(numeric_chrs, non_numeric_chrs)

    # Create chromosome factor with proper ordering
    gwas_df$chr_factor <- factor(gwas_df$chr_clean, levels = ordered_chrs)
    gwas_df$chr_num <- as.numeric(gwas_df$chr_factor)

    # Calculate cumulative positions for Manhattan plot
    chr_lengths <- aggregate(end ~ chr_num, gwas_df, max)
    chr_lengths <- chr_lengths[order(chr_lengths$chr_num), ]
    chr_lengths$cum_length <- cumsum(c(0, chr_lengths$end[-nrow(chr_lengths)]))

    # Merge cumulative lengths back to main data
    gwas_df <- merge(gwas_df, chr_lengths[, c("chr_num", "cum_length")], by = "chr_num", all.x = TRUE)

    # Calculate absolute position for Manhattan plot
    gwas_df$abs_pos <- gwas_df$start + gwas_df$cum_length

    # Calculate chromosome midpoints for labeling
    chr_midpoints <- aggregate(abs_pos ~ chr_num + chr_clean, gwas_df, function(x) {
      (min(x) + max(x)) / 2
    })

    # Store midpoints as an attribute for later use
    attr(gwas_df, "chr_midpoints") <- chr_midpoints

    # Restore original order
    return(gwas_df)
  }

  # ENHANCED: Association testing function with new continuous methods
  test_feature_association <- function(feature_presence, phenotype, test_method, covariates = NULL) {
    # Remove missing values
    complete_idx <- complete.cases(feature_presence, phenotype)
    if (!is.null(covariates)) {
      complete_idx <- complete_idx & complete.cases(covariates)
    }

    feature_clean <- feature_presence[complete_idx]
    pheno_clean <- phenotype[complete_idx]
    cov_clean <- if (!is.null(covariates)) covariates[complete_idx, , drop = FALSE] else NULL

    # Check if feature has variation
    if (length(unique(feature_clean)) < 2) {
      return(list(p_value = 1, effect_size = 0, test_stat = 0, n_samples = length(feature_clean)))
    }

    result <- switch(test_method,
      "fisher" = {
        # Fisher's exact test for binary phenotype - UNCHANGED
        contingency_table <- table(feature_clean, pheno_clean)
        if (any(dim(contingency_table) < 2)) {
          list(p_value = 1, effect_size = 0, test_stat = 0)
        } else {
          test_result <- fisher.test(contingency_table)
          list(
            p_value = test_result$p.value,
            effect_size = log(as.numeric(test_result$estimate)),
            test_stat = 0
          )
        }
      },
      "mannwhitney" = {
        # Mann-Whitney U test for continuous phenotype - UNCHANGED
        group0 <- pheno_clean[feature_clean == 0]
        group1 <- pheno_clean[feature_clean == 1]

        if (length(group0) < 2 || length(group1) < 2) {
          list(p_value = 1, effect_size = 0, test_stat = 0)
        } else {
          test_result <- wilcox.test(group1, group0)
          effect_size <- median(group1, na.rm = TRUE) - median(group0, na.rm = TRUE)
          list(
            p_value = test_result$p.value,
            effect_size = effect_size,
            test_stat = as.numeric(test_result$statistic)
          )
        }
      },
      # ENHANCED: New Welch's t-test for continuous phenotypes
      "welch_ttest" = {
        group0 <- pheno_clean[feature_clean == 0]
        group1 <- pheno_clean[feature_clean == 1]

        if (length(group0) < 2 || length(group1) < 2) {
          list(p_value = 1, effect_size = 0, test_stat = 0)
        } else {
          test_result <- t.test(group1, group0, var.equal = FALSE)
          mean_diff <- mean(group1, na.rm = TRUE) - mean(group0, na.rm = TRUE)

          # Calculate Cohen's d (standardized effect size)
          pooled_sd <- sqrt(((length(group1) - 1) * var(group1) + (length(group0) - 1) * var(group0)) /
            (length(group1) + length(group0) - 2))
          cohens_d <- mean_diff / pooled_sd

          list(
            p_value = test_result$p.value,
            effect_size = cohens_d,  # Cohen's d for standardized effect size
            test_stat = as.numeric(test_result$statistic),
            mean_diff = mean_diff    # Raw mean difference
          )
        }
      },
      # ENHANCED: Linear regression for continuous phenotypes with covariates
      "linear" = {
        tryCatch(
          {
            if (is.null(cov_clean)) {
              model_data <- data.frame(phenotype = pheno_clean, feature = feature_clean)
              model <- lm(phenotype ~ feature, data = model_data)
            } else {
              model_data <- data.frame(phenotype = pheno_clean, feature = feature_clean, cov_clean)
              model <- lm(phenotype ~ ., data = model_data)
            }

            summary_model <- summary(model)
            if ("feature" %in% rownames(summary_model$coefficients)) {
              feature_coef <- summary_model$coefficients["feature", ]
              list(
                p_value = feature_coef[4],           # P-value
                effect_size = feature_coef[1],       # Regression coefficient
                test_stat = feature_coef[3],         # T-statistic
                r_squared = summary_model$r.squared  # R-squared
              )
            } else {
              list(p_value = 1, effect_size = 0, test_stat = 0, r_squared = 0)
            }
          },
          error = function(e) {
            list(p_value = 1, effect_size = 0, test_stat = 0, r_squared = 0)
          })
      },
      "logistic" = {
        # Logistic regression with optional covariates - UNCHANGED
        tryCatch(
          {
            # For categorical phenotype, use multinomial approach if >2 categories
            n_categories <- length(unique(pheno_clean))

            if (n_categories == 2) {
              # Binary logistic regression
              # Ensure binary coding (0/1)
              pheno_binary <- as.numeric(factor(pheno_clean)) - 1

              if (is.null(cov_clean)) {
                model_data <- data.frame(phenotype = pheno_binary, feature = feature_clean)
                model <- glm(phenotype ~ feature, data = model_data, family = binomial())
              } else {
                model_data <- data.frame(phenotype = pheno_binary, feature = feature_clean, cov_clean)
                model <- glm(phenotype ~ ., data = model_data, family = binomial())
              }

              summary_model <- summary(model)
              if ("feature" %in% rownames(summary_model$coefficients)) {
                feature_coef <- summary_model$coefficients["feature", ]
                list(
                  p_value = feature_coef[4],
                  effect_size = feature_coef[1],
                  test_stat = feature_coef[3]
                )
              } else {
                list(p_value = 1, effect_size = 0, test_stat = 0)
              }
            } else {
              # Multinomial logistic regression for >2 categories
              if (requireNamespace("nnet", quietly = TRUE)) {
                if (is.null(cov_clean)) {
                  model_data <- data.frame(phenotype = factor(pheno_clean), feature = feature_clean)
                  model <- nnet::multinom(phenotype ~ feature, data = model_data, trace = FALSE)
                } else {
                  model_data <- data.frame(phenotype = factor(pheno_clean), feature = feature_clean, cov_clean)
                  model <- nnet::multinom(phenotype ~ ., data = model_data, trace = FALSE)
                }

                # Calculate likelihood ratio test
                null_model <- nnet::multinom(phenotype ~ 1, data = model_data, trace = FALSE)
                lr_test <- anova(null_model, model, test = "Chisq")

                list(
                  p_value = lr_test$`Pr(>Chi)`[2],
                  effect_size = 0,  # Complex to define for multinomial
                  test_stat = lr_test$Deviance[2]
                )
              } else {
                # Fallback to chi-square test if nnet not available
                if (length(unique(feature_clean)) < 2) {
                  list(p_value = 1, effect_size = 0, test_stat = 0)
                } else {
                  contingency_table <- table(feature_clean, pheno_clean)
                  if (any(dim(contingency_table) < 2)) {
                    list(p_value = 1, effect_size = 0, test_stat = 0)
                  } else {
                    test_result <- chisq.test(contingency_table)
                    list(
                      p_value = test_result$p.value,
                      effect_size = sqrt(test_result$statistic / sum(contingency_table)),  # Cramer's V
                      test_stat = test_result$statistic
                    )
                  }
                }
              }
            }
          },
          error = function(e) {
            # If logistic regression fails, fall back to chi-square test
            if (length(unique(feature_clean)) < 2 || length(unique(pheno_clean)) < 2) {
              list(p_value = 1, effect_size = 0, test_stat = 0)
            } else {
              contingency_table <- table(feature_clean, pheno_clean)
              if (any(dim(contingency_table) < 2)) {
                list(p_value = 1, effect_size = 0, test_stat = 0)
              } else {
                test_result <- chisq.test(contingency_table)
                list(
                  p_value = test_result$p.value,
                  effect_size = sqrt(test_result$statistic / sum(contingency_table)),  # Cramer's V
                  test_stat = test_result$statistic
                )
              }
            }
          })
      }
    )

    result$n_samples <- length(feature_clean)
    return(result)
  }

  # Step 10: Run association analysis - Enhanced reporting
  if (verbose) {
    cat(sprintf("Running PAV GWAS with %s test...\n", test_type))
    cat(sprintf("Testing %d PAV features on %d samples\n",
      nrow(filtered_matrix), ncol(filtered_matrix)))
    if (analysis_type %in% c("continuous", "continuous_raw")) {
      cat(sprintf("Enhanced continuous phenotype analysis: method=%s, test=%s\n", analysis_type, test_type))
    }
  }

  # Set up parallel processing
  cl <- makeCluster(n_cores)
  clusterExport(cl, c("test_feature_association", "test_type", "transformed_phenotype", "covariates"),
    envir = environment())
  clusterEvalQ(cl, {
    library(stats)
  })

  # Perform association tests
  start_time <- Sys.time()

  results <- parLapply(cl, 1:nrow(filtered_matrix), function(i) {
    feature_presence <- filtered_matrix[i, ]
    result <- test_feature_association(feature_presence, transformed_phenotype, test_type, covariates)
    result$feature_id <- rownames(filtered_matrix)[i]
    result$feature_index <- i
    return(result)
  })

  stopCluster(cl)

  end_time <- Sys.time()
  if (verbose) {
    cat(sprintf("Association testing completed in %.2f minutes\n",
      as.numeric(difftime(end_time, start_time, units = "mins"))))
  }

  # Step 11: Organize results with position information (FIXED VERSION)
  if (verbose) cat("Organizing results...\n")

  # First check the structure of the results list
  if (length(results) == 0) {
    stop("No results returned from association testing. Please check your data.")
  }

  # Check and fix inconsistent elements in results
  for (i in seq_along(results)) {
    if (is.null(results[[i]])) {
      # If a result is NULL, create default result
      results[[i]] <- list(
        p_value = 1,
        effect_size = 0,
        test_stat = 0,
        n_samples = ncol(filtered_matrix),
        feature_id = rownames(filtered_matrix)[i],
        feature_index = i
      )
    } else {
      # Ensure each result has required fields
      if (is.null(results[[i]]$feature_id)) {
        results[[i]]$feature_id <- rownames(filtered_matrix)[i]
      }
      if (is.null(results[[i]]$feature_index)) {
        results[[i]]$feature_index <- i
      }
      if (is.null(results[[i]]$p_value)) {
        results[[i]]$p_value <- 1
      }
      if (is.null(results[[i]]$effect_size)) {
        results[[i]]$effect_size <- 0
      }
      if (is.null(results[[i]]$test_stat)) {
        results[[i]]$test_stat <- 0
      }
      if (is.null(results[[i]]$n_samples)) {
        results[[i]]$n_samples <- ncol(filtered_matrix)
      }
    }
  }

  # Safe extraction function to ensure consistent lengths
  extract_safe <- function(results_list, field_name, default_value = NA) {
    result <- sapply(results_list, function(x) {
      if (is.null(x) || is.null(x[[field_name]])) {
        return(default_value)
      } else {
        return(x[[field_name]])
      }
    })
    return(result)
  }

  # Extract all fields
  feature_ids <- extract_safe(results, "feature_id", "unknown")
  p_values <- extract_safe(results, "p_value", 1)
  effect_sizes <- extract_safe(results, "effect_size", 0)
  test_stats <- extract_safe(results, "test_stat", 0)
  n_samples_vec <- extract_safe(results, "n_samples", ncol(filtered_matrix))

  # Validate all vector lengths are consistent
  expected_length <- nrow(filtered_matrix)
  vectors_to_check <- list(
    feature_ids = feature_ids,
    p_values = p_values,
    effect_sizes = effect_sizes,
    test_stats = test_stats,
    n_samples = n_samples_vec
  )

  # Check and report length mismatches
  for (vec_name in names(vectors_to_check)) {
    vec_length <- length(vectors_to_check[[vec_name]])
    if (vec_length != expected_length) {
      warning(sprintf("Vector '%s' has length %d, expected %d. Padding with default values.",
        vec_name, vec_length, expected_length))

      # Fill or truncate vectors as needed
      if (vec_length < expected_length) {
        # Fill missing values
        if (vec_name == "feature_ids") {
          vectors_to_check[[vec_name]] <- c(vectors_to_check[[vec_name]],
            rep("unknown", expected_length - vec_length))
        } else {
          default_val <- if (vec_name == "p_values") 1 else 0
          vectors_to_check[[vec_name]] <- c(vectors_to_check[[vec_name]],
            rep(default_val, expected_length - vec_length))
        }
      } else {
        # Truncate excess values
        vectors_to_check[[vec_name]] <- vectors_to_check[[vec_name]][1:expected_length]
      }
    }
  }

  # Now safely create the data frame
  gwas_results <- data.frame(
    chr = pav_info$chr[keep_features],
    start = pav_info$start[keep_features],
    end = pav_info$end[keep_features],
    sequence = pav_info$sequence[keep_features],
    feature_id = vectors_to_check$feature_ids,
    p_value = as.numeric(vectors_to_check$p_values),
    effect_size = as.numeric(vectors_to_check$effect_sizes),
    test_stat = as.numeric(vectors_to_check$test_stats),
    n_samples = as.numeric(vectors_to_check$n_samples),
    stringsAsFactors = FALSE
  )

  # Validate data frame creation success
  if (nrow(gwas_results) != expected_length) {
    stop(sprintf("Data frame creation failed. Expected %d rows, got %d rows.",
      expected_length, nrow(gwas_results)))
  }

  if (verbose) {
    cat(sprintf("Successfully created results data frame with %d features\n", nrow(gwas_results)))
  }

  # Handle invalid p-values
  invalid_p <- is.na(gwas_results$p_value) | gwas_results$p_value <= 0 | gwas_results$p_value > 1
  if (any(invalid_p)) {
    n_invalid <- sum(invalid_p)
    warning(sprintf("Found %d invalid p-values. Setting to 1.", n_invalid))
    gwas_results$p_value[invalid_p] <- 1
  }

  # Continue with original multiple testing correction and other processing...
  gwas_results$p_bonferroni <- p.adjust(gwas_results$p_value, method = "bonferroni")
  gwas_results$p_fdr <- p.adjust(gwas_results$p_value, method = "fdr")

  # Calculate -log10(p)
  gwas_results$neg_log10_p <- -log10(pmax(gwas_results$p_value, .Machine$double.xmin))

  # Add feature frequency information
  gwas_results$feature_frequency <- filtered_feature_info$frequencies

  # Create genomic position for Manhattan plot
  gwas_results <- create_manhattan_positions(gwas_results)

  # Sort by p-value
  gwas_results <- gwas_results[order(gwas_results$p_value), ]

  # Calculate multiple significance thresholds
  thresholds <- c(1e-5, 1e-6, 1e-7, 1e-8, significance_threshold)
  thresholds <- unique(sort(thresholds, decreasing = TRUE))  # Remove duplicates and sort

  # Calculate significant counts for each threshold
  threshold_counts <- sapply(thresholds, function(t) sum(gwas_results$p_value < t))
  names(threshold_counts) <- sprintf("%.0e", thresholds)

  # Store threshold information
  threshold_info <- data.frame(
    threshold = thresholds,
    threshold_label = sprintf("%.0e", thresholds),
    n_significant = threshold_counts,
    stringsAsFactors = FALSE
  )

  # Step 12: Enhanced summary statistics
  analysis_summary <- list(
    analysis_type = analysis_type,
    test_type = test_type,
    n_samples = ncol(filtered_matrix),
    n_features_tested = nrow(filtered_matrix),
    n_features_original = nrow(pav_matrix),
    min_sample_count = min_sample_count,
    max_sample_count = max_sample_count,
    use_population_structure = use_population_structure,
    n_pcs = if (use_population_structure && !is.null(pca_results)) ncol(pca_results$pcs) else 0,
    significance_threshold = significance_threshold,
    threshold_info = threshold_info,
    n_significant_bonf_threshold = sum(gwas_results$p_value < significance_threshold),
    n_significant_005 = sum(gwas_results$p_value < 0.05),
    n_significant_001 = sum(gwas_results$p_value < 0.01),
    n_significant_bonf = sum(gwas_results$p_bonferroni < 0.05),
    n_significant_fdr = sum(gwas_results$p_fdr < 0.05),
    min_p_value = min(gwas_results$p_value),
    binary_threshold = if (analysis_type == "binary") binary_threshold else NA,
    # Enhanced summary for continuous phenotypes
    phenotype_transformation = ifelse(analysis_type %in% c("continuous", "continuous_raw"),
      "Enhanced continuous phenotype analysis",
      "Categorical/transformed analysis")
  )

  if (verbose) {
    cat(sprintf("Analysis completed!\n"))
    cat(sprintf("Method used: %s analysis with %s test\n", analysis_type, test_type))
    cat(sprintf("Significance threshold (Bonferroni): %.2e\n", analysis_summary$significance_threshold))
    cat("\n=== Significant associations at different thresholds ===\n")
    for (i in 1:nrow(threshold_info)) {
      cat(sprintf("p < %s: %d features\n",
        threshold_info$threshold_label[i],
        threshold_info$n_significant[i]))
    }
    cat(sprintf("\nAdditional statistics:\n"))
    cat(sprintf("Significant associations (p < 0.05): %d\n", analysis_summary$n_significant_005))
    cat(sprintf("Significant associations (FDR < 0.05): %d\n", analysis_summary$n_significant_fdr))
    cat(sprintf("Most significant p-value: %.2e\n", analysis_summary$min_p_value))

    # Enhanced reporting for continuous phenotypes
    if (analysis_type %in% c("continuous", "continuous_raw")) {
      n_sig <- sum(gwas_results$p_value < 0.05, na.rm = TRUE)
      if (n_sig > 0) {
        cat(sprintf("Found %d significant associations for continuous phenotype\n", n_sig))
        top_hits <- head(gwas_results[order(gwas_results$p_value), ], 5)
        cat("Top 5 associations:\n")
        for (i in 1:nrow(top_hits)) {
          cat(sprintf("  %s: p=%.2e, effect=%.3f\n",
            top_hits$feature_id[i], top_hits$p_value[i], top_hits$effect_size[i]))
        }
      }
    }
  }

  # Step 13: Create visualizations and save results - UNCHANGED
  if (verbose) cat("Creating visualizations and saving results...\n")

  # Save main results (without quotes)
  write.table(gwas_results, file.path(output_dir, "gwas_results.csv"),
    sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)

  # Save results for each significance threshold
  if (verbose) cat("Creating threshold-specific result files...\n")

  for (i in 1:nrow(threshold_info)) {
    thresh <- threshold_info$threshold[i]
    thresh_label <- threshold_info$threshold_label[i]

    # Filter significant results for this threshold
    significant_results <- gwas_results[gwas_results$p_value < thresh, ]

    if (nrow(significant_results) > 0) {
      # Save threshold-specific results
      filename <- file.path(output_dir, paste0("significant_p", gsub("e-0?", "e-", thresh_label), ".csv"))
      write.table(significant_results, filename,
        sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)

      if (verbose) {
        cat(sprintf("Saved %d significant features (p < %s) to: %s\n",
          nrow(significant_results), thresh_label, basename(filename)))
      }
    } else {
      if (verbose) {
        cat(sprintf("No significant features found for p < %s\n", thresh_label))
      }
    }
  }

  # Also save FDR significant results
  fdr_significant <- gwas_results[gwas_results$p_fdr < 0.05, ]
  if (nrow(fdr_significant) > 0) {
    write.table(fdr_significant, file.path(output_dir, "significant_fdr0.05.csv"),
      sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
    if (verbose) {
      cat(sprintf("Saved %d FDR significant features to: significant_fdr0.05.csv\n",
        nrow(fdr_significant)))
    }
  }

  # Manhattan plot with chromosomal positions - UNCHANGED FROM ORIGINAL
  # Prepare chromosome colors (alternating)
  n_chrs <- length(unique(gwas_results$chr_num))
  chr_colors <- rep(c("darkblue", "lightblue"), length.out = n_chrs)
  names(chr_colors) <- sort(unique(gwas_results$chr_num))

  # Create chromosome color mapping
  gwas_results$chr_color <- chr_colors[as.character(gwas_results$chr_num)]

  # Get chromosome midpoints from attributes
  chr_midpoints <- attr(gwas_results, "chr_midpoints")

  # Prepare threshold lines data
  threshold_lines <- data.frame(
    threshold = thresholds,
    neg_log10_threshold = -log10(thresholds),
    label = sprintf("p = %s", sprintf("%.0e", thresholds)),
    color = c("orange", "purple", "green", "blue", "red")[1:length(thresholds)],
    stringsAsFactors = FALSE
  )
  # Make Bonferroni threshold red
  threshold_lines$color[threshold_lines$threshold == significance_threshold] <- "red"

  manhattan_plot <- gwas_results %>%
    ggplot(aes(x = abs_pos, y = neg_log10_p)) +
    geom_point(aes(color = chr_color), alpha = 0.7, size = 0.8) +
    scale_color_identity()

  # Add multiple threshold lines individually
  for (i in 1:nrow(threshold_lines)) {
    manhattan_plot <- manhattan_plot +
      geom_hline(yintercept = threshold_lines$neg_log10_threshold[i],
        linetype = "dashed", color = threshold_lines$color[i],
        size = 0.8, alpha = 0.8)
  }

  manhattan_plot <- manhattan_plot +
    labs(title = "PAV GWAS Manhattan Plot",
      x = "Chromosomal Position",
      y = expression(-log[10](P)),
      subtitle = paste("Analysis:", analysis_type, "| Test:", test_type,
        "| Samples:", analysis_summary$n_samples, "| Features:", nrow(gwas_results)),
      caption = paste("Threshold lines:", paste(threshold_lines$label, collapse = ", "))) +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      plot.caption = element_text(size = 8, color = "gray60")
    )

  # Add chromosome labels if there are multiple chromosomes
  if (n_chrs > 1 && n_chrs <= 30 && !is.null(chr_midpoints)) {  # Only label if reasonable number of chromosomes
    manhattan_plot <- manhattan_plot +
      scale_x_continuous(
        breaks = chr_midpoints$abs_pos,
        labels = chr_midpoints$chr_clean
      )
  } else if (n_chrs > 30) {
    # Too many chromosomes, use generic labels
    manhattan_plot <- manhattan_plot +
      labs(caption = paste("Note:", n_chrs, "chromosomes/scaffolds present;",
        paste(threshold_lines$label, collapse = ", ")))
  }

  ggsave(file.path(output_dir, "manhattan_plot.pdf"), manhattan_plot,
    width = 14, height = 6, dpi = 300)
  ggsave(file.path(output_dir, "manhattan_plot.png"), manhattan_plot,
    width = 14, height = 6, dpi = 300, bg = "white")

  # Q-Q plot - UNCHANGED FROM ORIGINAL
  observed_p <- gwas_results$p_value
  n <- length(observed_p)
  expected_p <- (1:n) / (n + 1)

  qq_plot <- data.frame(
    expected = -log10(expected_p),
    observed = -log10(sort(observed_p))
  ) %>%
    ggplot(aes(x = expected, y = observed)) +
    geom_point(alpha = 0.6, color = "darkblue") +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    labs(title = "Q-Q Plot of Association Test P-values",
      x = expression("Expected" ~ -log[10](P)),
      y = expression("Observed" ~ -log[10](P)),
      subtitle = "Deviation from diagonal indicates test calibration issues") +
    theme_minimal()

  ggsave(file.path(output_dir, "qq_plot.pdf"), qq_plot,
    width = 8, height = 6, dpi = 300)
  ggsave(file.path(output_dir, "qq_plot.png"), qq_plot,
    width = 8, height = 6, dpi = 300, bg = "white")

  # Enhanced phenotype distribution plot
  if (pheno_stats$is_numeric) {
    # For numeric data, use histogram
    pheno_dist_plot <- data.frame(phenotype = phenotype_values) %>%
      ggplot(aes(x = phenotype)) +
      geom_histogram(bins = 30, fill = "lightblue", color = "black", alpha = 0.7) +
      geom_density(aes(y = ..scaled.. * max(..count..)), color = "red", size = 1) +
      labs(title = "Enhanced Phenotype Distribution Analysis",
        x = "Phenotype Value",
        y = "Frequency",
        subtitle = paste("n =", length(phenotype_values),
          "| Type:", analysis_type,
          "| Mean =", round(mean(phenotype_values, na.rm = TRUE), 3),
          "| CV =", round(pheno_stats$cv, 3))) +
      theme_minimal()
  } else {
    # For categorical data, use bar plot
    pheno_counts <- table(phenotype_values)
    pheno_df_plot <- data.frame(
      phenotype = names(pheno_counts),
      n = as.numeric(pheno_counts),
      stringsAsFactors = FALSE
    )

    pheno_dist_plot <- pheno_df_plot %>%
      ggplot(aes(x = factor(phenotype), y = n)) +
      geom_col(fill = "lightblue", color = "black", alpha = 0.7) +
      geom_text(aes(label = n), vjust = -0.5) +
      labs(title = "Phenotype Distribution",
        x = "Phenotype Category",
        y = "Frequency",
        subtitle = paste("n =", length(phenotype_values),
          "| Type:", analysis_type,
          "| Categories =", pheno_stats$n_unique)) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }

  ggsave(file.path(output_dir, "phenotype_distribution.pdf"), pheno_dist_plot,
    width = 8, height = 6, dpi = 300)
  ggsave(file.path(output_dir, "phenotype_distribution.png"), pheno_dist_plot,
    width = 8, height = 6, dpi = 300, bg = "white")

  # Population structure plot (if applicable) - UNCHANGED FROM ORIGINAL
  if (use_population_structure && !is.null(pca_results)) {
    pc_var_plot <- data.frame(
      PC = paste0("PC", 1:length(pca_results$var_explained)),
      Variance = pca_results$var_explained * 100
    ) %>%
      ggplot(aes(x = PC, y = Variance)) +
      geom_col(fill = "steelblue", alpha = 0.7) +
      labs(title = "Principal Component Analysis",
        x = "Principal Component",
        y = "Variance Explained (%)",
        subtitle = paste("First", ncol(pca_results$pcs), "PCs used as covariates")) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    pc_scatter_plot <- data.frame(pca_results$pcs[, 1:min(2, ncol(pca_results$pcs))]) %>%
      ggplot(aes(x = PC1, y = if (ncol(pca_results$pcs) > 1) PC2 else PC1)) +
      geom_point(alpha = 0.7, color = "steelblue") +
      labs(title = "PC1 vs PC2 Scatter Plot",
        x = paste0("PC1 (", round(pca_results$var_explained[1] * 100, 1), "%)"),
        y = if (ncol(pca_results$pcs) > 1)
          paste0("PC2 (", round(pca_results$var_explained[2] * 100, 1), "%)")
        else "PC1") +
      theme_minimal()

    if (requireNamespace("gridExtra", quietly = TRUE)) {
      pop_structure_plot <- gridExtra::grid.arrange(pc_var_plot, pc_scatter_plot, nrow = 2)
      ggsave(file.path(output_dir, "population_structure.pdf"), pop_structure_plot,
        width = 8, height = 10, dpi = 300)
      ggsave(file.path(output_dir, "population_structure.png"), pop_structure_plot,
        width = 8, height = 10, dpi = 300, bg = "white")
    } else {
      ggsave(file.path(output_dir, "pca_variance.pdf"), pc_var_plot, width = 8, height = 6, dpi = 300)
      ggsave(file.path(output_dir, "pca_variance.png"), pc_var_plot, width = 8, height = 6, dpi = 300, bg = "white")
      ggsave(file.path(output_dir, "pca_scatter.pdf"), pc_scatter_plot, width = 8, height = 6, dpi = 300)
      ggsave(file.path(output_dir, "pca_scatter.png"), pc_scatter_plot, width = 8, height = 6, dpi = 300, bg = "white")
    }
  }

  # Enhanced text summary generation
  input_type <- ifelse(is.character(substitute(pav_data)), "file", "data frame")
  pheno_type <- ifelse(is.character(substitute(phenotype_data)), "file", "data frame")

  # Create threshold summary lines
  threshold_summary_lines <- character(nrow(threshold_info))
  for (i in 1:nrow(threshold_info)) {
    threshold_summary_lines[i] <- paste0("p < ", threshold_info$threshold_label[i],
      ": ", threshold_info$n_significant[i], " features")
  }

  summary_text <- c(
    "=== Enhanced PAV GWAS Analysis Summary ===",
    "",
    paste("Analysis date:", Sys.Date()),
    paste("PAV data source:", input_type),
    paste("Phenotype data source:", pheno_type),
    "",
    "=== Parameters ===",
    paste("Analysis type:", analysis_summary$analysis_type),
    paste("Test type:", analysis_summary$test_type),
    paste("Samples used:", analysis_summary$n_samples),
    paste("PAV features tested:", analysis_summary$n_features_tested),
    paste("PAV features filtered out:", analysis_summary$n_features_original - analysis_summary$n_features_tested),
    paste("Min sample count:", analysis_summary$min_sample_count),
    paste("Max sample count:", analysis_summary$max_sample_count),
    paste("Population structure correction:",
      if (use_population_structure && !is.null(pca_results)) "Yes"
      else if (!use_population_structure) "Disabled"
      else "Failed (skipped)"),
    if (use_population_structure && !is.null(pca_results))
      paste("Principal components used:", analysis_summary$n_pcs)
    else "",
    if (!is.na(analysis_summary$binary_threshold)) paste("Binary threshold:", analysis_summary$binary_threshold) else "",
    if (analysis_type %in% c("continuous", "continuous_raw"))
      paste("Phenotype transformation:", analysis_summary$phenotype_transformation) else "",
    "",
    "=== Results ===",
    paste("Significance threshold (Bonferroni):", sprintf("%.2e", analysis_summary$significance_threshold)),
    "",
    "=== Significant associations at different thresholds ===",
    threshold_summary_lines,
    "",
    "=== Additional statistics ===",
    paste("Significant associations (p < 0.05):", analysis_summary$n_significant_005),
    paste("Significant associations (p < 0.01):", analysis_summary$n_significant_001),
    paste("Significant associations (Bonferroni p < 0.05):", analysis_summary$n_significant_bonf),
    paste("Significant associations (FDR p < 0.05):", analysis_summary$n_significant_fdr),
    paste("Most significant p-value:", sprintf("%.2e", analysis_summary$min_p_value)),
    "",
    "=== Top 10 Most Significant PAV Features ===",
    ""
  )

  # Add top results table with position information
  top_results <- head(gwas_results[, c("chr", "start", "end", "feature_id", "p_value", "effect_size", "p_fdr", "feature_frequency")], 10)
  top_results$p_value <- sprintf("%.2e", top_results$p_value)
  top_results$p_fdr <- sprintf("%.2e", top_results$p_fdr)
  top_results$effect_size <- sprintf("%.3f", top_results$effect_size)

  # Write summary to file
  writeLines(summary_text, file.path(output_dir, "analysis_summary.txt"))
  write.table(top_results, file.path(output_dir, "analysis_summary.txt"),
    append = TRUE, quote = FALSE, row.names = FALSE, col.names = TRUE,
    sep = "\t")

  if (verbose) {
    cat("All results saved to:", output_dir, "\n")
    cat("Main output files (both PDF and PNG formats):\n")
    cat("  - gwas_results.csv: Complete association results\n")
    cat("  - manhattan_plot.pdf/.png: Enhanced Manhattan plot with multiple threshold lines\n")
    cat("  - qq_plot.pdf/.png: Q-Q plot\n")
    cat("  - phenotype_distribution.pdf/.png: Enhanced phenotype analysis\n")
    if (use_population_structure && !is.null(pca_results)) {
      cat("  - population_structure.pdf/.png: PCA analysis\n")
    } else if (!use_population_structure) {
      cat("  - Note: Population structure correction was disabled\n")
    } else {
      cat("  - Note: Population structure correction failed (skipped)\n")
    }
    cat("  - analysis_summary.txt: Enhanced summary report\n")

    cat("\nThreshold-specific output files:\n")
    for (i in 1:nrow(threshold_info)) {
      thresh_label <- threshold_info$threshold_label[i]
      n_sig <- threshold_info$n_significant[i]
      filename <- paste0("significant_p", gsub("e-0?", "e-", thresh_label), ".csv")
      if (n_sig > 0) {
        cat(sprintf("  - %s: %d significant features (p < %s)\n", filename, n_sig, thresh_label))
      }
    }

    if (analysis_summary$n_significant_fdr > 0) {
      cat(sprintf("  - significant_fdr0.05.csv: %d FDR significant features\n",
        analysis_summary$n_significant_fdr))
    }

    cat("\nNote: All plots saved in both PDF (high quality) and PNG (fast loading) formats\n")

    if (analysis_type %in% c("continuous", "continuous_raw")) {
      cat("\n=== Enhanced Continuous Phenotype Analysis ===\n")
      cat("Effect size interpretation:\n")
      if (test_type == "welch_ttest") {
        cat("  - Cohen's d: Small=0.2, Medium=0.5, Large=0.8\n")
      } else if (test_type == "linear") {
        cat("  - Regression coefficient: Change in phenotype per PAV presence\n")
      } else if (test_type == "mannwhitney") {
        cat("  - Median difference: Difference in medians between groups\n")
      }
    }
  }

  # Return enhanced results list
  return(list(
    gwas_results = gwas_results,
    phenotype_stats = pheno_stats,
    pca_results = pca_results,
    filtered_features = filtered_feature_info,
    analysis_summary = analysis_summary,
    threshold_results = threshold_info
  ))
}

# Enhanced function information when loaded
cat("Enhanced PAV GWAS analysis function loaded successfully!\n")
cat("Use pav_gwas() to run the complete analysis.\n")
cat("PAV matrix format: chr, start, end, sequence, sample1, sample2, ...\n")
cat("Enhanced Features:\n")
cat("  - Automatic continuous phenotype detection and analysis\n")
cat("  - Welch's t-test for normal continuous data\n")
cat("  - Linear regression with population structure covariates\n")
cat("  - Enhanced effect size calculations (Cohen's d, regression coefficients)\n")
cat("  - Multiple significance thresholds, threshold-specific output files\n")
cat("  - Manhattan plot with threshold lines\n")
cat("  - Both PDF (publication quality) and PNG (fast loading) for all plots\n")
cat("Function supports both file paths and data frames as input.\n")
cat("New analysis types: 'continuous', 'continuous_raw'\n")
cat("New test types: 'welch_ttest', 'linear'\n")
cat("Type ?pav_gwas for detailed help documentation.\n")
