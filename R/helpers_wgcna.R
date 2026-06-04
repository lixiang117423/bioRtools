#' WGCNA Helper Functions
#'
#' Internal helper functions to support modular WGCNA workflow.
#' These functions break down the complex WGCNA analysis into manageable steps.
#'
#' @keywords internal
#' @noRd

# ============================================================================
# WGCNA Input Validation
# ============================================================================

#' Validate WGCNA input parameters
#'
#' @param expression_data Expression data matrix/data frame
#' @param sample_info Sample information data frame
#' @param trait_columns Character vector of trait column names
#' @param min_gene_sd Minimum gene standard deviation
#' @param tom_type Type of TOM matrix
#' @param correlation_method Correlation method
#'
#' @return Invisibly returns TRUE if valid
#' @keywords internal
#' @noRd
validate_wgcna_inputs <- function(expression_data, sample_info, trait_columns,
                                  min_gene_sd = 0.1, tom_type = "signed",
                                  correlation_method = "pearson") {
  if (missing(expression_data) || is.null(expression_data)) {
    stop(err_missing_required("expression_data"))
  }

  if (missing(sample_info) || is.null(sample_info)) {
    stop(err_missing_required("sample_info"))
  }

  if (missing(trait_columns) || length(trait_columns) == 0) {
    stop(err_missing_required("trait_columns"))
  }

  if (!is.character(trait_columns)) {
    stop(err_invalid_input("trait_columns", "a character vector"))
  }

  if (!is.numeric(min_gene_sd) || min_gene_sd < 0) {
    stop(err_invalid_threshold("min_gene_sd", min_gene_sd, "[0, ∞)"))
  }

  if (!tom_type %in% c("signed", "unsigned", "signed hybrid")) {
    stop(err_invalid_input("tom_type", "one of: 'signed', 'unsigned', 'signed hybrid'"))
  }

  if (!correlation_method %in% c("pearson", "spearman", "kendall")) {
    stop(err_invalid_input("correlation_method", "one of: 'pearson', 'spearman', 'kendall'"))
  }

  invisible(TRUE)
}

# ============================================================================
# WGCNA Data Preparation
# ============================================================================

#' Filter low-variance genes for WGCNA
#'
#' @param expression_data Expression matrix with genes as rows, samples as columns
#' @param min_sd Minimum standard deviation threshold
#' @param verbose Print progress messages
#'
#' @return Filtered expression matrix
#' @keywords internal
#' @noRd
filter_wgcna_genes <- function(expression_data, min_sd = 0.1, verbose = TRUE) {
  gene_sds <- apply(expression_data, 1, stats::sd, na.rm = TRUE)
  keep_genes <- gene_sds >= min_sd

  n_removed <- sum(!keep_genes)
  n_kept <- sum(keep_genes)

  if (verbose) {
    cat(sprintf("Filtering genes by standard deviation (min_sd = %.2f):\n", min_sd))
    cat(sprintf("  Removed: %d genes\n", n_removed))
    cat(sprintf("  Kept: %d genes\n\n", n_kept))
  }

  if (n_kept == 0) {
    stop("No genes remain after filtering. Consider lowering min_sd threshold.")
  }

  expression_data[keep_genes, ]
}

#' Prepare trait data for WGCNA
#'
#' @param sample_info Sample information data frame
#' @param trait_columns Column names of traits
#' @param verbose Print progress messages
#'
#' @return Data frame of traits with samples as rows
#' @keywords internal
#' @noRd
prepare_wgcna_traits <- function(sample_info, trait_columns, verbose = TRUE) {
  # Validate trait columns exist
  missing_cols <- setdiff(trait_columns, colnames(sample_info))
  if (length(missing_cols) > 0) {
    stop(err_missing_design_variable(missing_cols[1], colnames(sample_info)))
  }

  traits <- sample_info[, trait_columns, drop = FALSE]

  # Convert to numeric if possible
  traits <- traits %>%
    dplyr::mutate(dplyr::across(dplyr::everything(),
      ~ tryCatch(as.numeric(.), warning = function(w) .)))

  if (verbose) {
    cat("Prepared trait data:\n")
    print(str(traits))
    cat("\n")
  }

  traits
}

# ============================================================================
# WGCNA Soft-thresholding
# ============================================================================

#' Perform soft-thresholding analysis
#'
#' Tests a range of soft-thresholding powers to find optimal scale-free topology.
#'
#' @param expression_data Expression matrix (genes x samples)
#' @param power_vector Vector of powers to test
#' @param correlation_method Correlation method ("pearson", "spearman", etc.)
#' @param tom_type Type of TOM matrix
#' @param verbose Print progress
#'
#' @return Power selection results
#' @keywords internal
#' @noRd
analyze_softthreshold_power <- function(expression_data, power_vector = 1:30,
                                        correlation_method = "pearson",
                                        tom_type = "signed", verbose = TRUE) {
  if (verbose) {
    cat("Analyzing soft-thresholding powers for scale-free topology...\n")
    cat("This may take a few minutes for large datasets.\n\n")
  }

  # Use WGCNA's built-in function
  powers_analysis <- WGCNA::pickSoftThreshold(
    expression_data,
    powerVector = power_vector,
    corFnc = correlation_method,
    networkType = tom_type,
    verbose = if (verbose) 2 else 0
  )

  if (verbose) {
    cat("Soft-thresholding analysis complete.\n")
    cat(sprintf("Recommended power: %d (R² = %.3f)\n\n",
      powers_analysis$powerEstimate,
      powers_analysis$fitIndices[
        which.min(abs(powers_analysis$fitIndices$SFT.R.sq - 0.9)),
        "SFT.R.sq"
      ]))
  }

  powers_analysis
}

# ============================================================================
# Network Construction
# ============================================================================

#' Build correlation network from expression data
#'
#' @param expression_data Expression matrix
#' @param soft_power Soft-thresholding power
#' @param correlation_method Correlation method
#' @param tom_type Type of TOM matrix
#' @param max_block_size Maximum block size for calculation
#' @param verbose Print progress
#'
#' @return Adjacency matrix or TOM matrix
#' @keywords internal
#' @noRd
build_correlation_network <- function(expression_data, soft_power = 6,
                                      correlation_method = "pearson",
                                      tom_type = "signed",
                                      max_block_size = NULL, verbose = TRUE) {
  if (verbose) {
    cat(sprintf("Building %s TOM network with power = %d...\n",
      tom_type, soft_power))
  }

  # Create adjacency matrix
  adjacency <- WGCNA::adjacency(
    expression_data,
    power = soft_power,
    corFnc = correlation_method,
    type = tom_type
  )

  if (verbose) {
    cat(sprintf("Adjacency matrix created (%d x %d).\n",
      nrow(adjacency), ncol(adjacency)))
  }

  adjacency
}

# ============================================================================
# Module Detection
# ============================================================================

#' Detect co-expression modules
#'
#' @param tom_matrix TOM matrix (topological overlap matrix)
#' @param min_module_size Minimum module size
#' @param merge_cutoff Merge height cutoff for similar modules
#' @param verbose Print progress
#'
#' @return List with module assignments and dendrogram
#' @keywords internal
#' @noRd
detect_modules <- function(tom_matrix, min_module_size = 50,
                          merge_cutoff = 0.25, verbose = TRUE) {
  if (verbose) {
    cat(sprintf("Detecting modules with min_module_size = %d...\n",
      min_module_size))
  }

  # Calculate dissimilarity
  dissimilarity <- 1 - tom_matrix

  # Hierarchical clustering
  gene_tree <- stats::hclust(stats::as.dist(dissimilarity), method = "average")

  # Cut dendrogram to detect modules
  dynamic_colors <- WGCNA::cutreeDynamic(
    dendro = gene_tree,
    distM = dissimilarity,
    deepSplit = 2,
    pamStage = FALSE,
    minClusterSize = min_module_size
  )

  if (verbose) {
    n_modules <- length(unique(dynamic_colors[dynamic_colors != 0]))
    cat(sprintf("Detected %d modules.\n\n", n_modules))
  }

  list(
    module_colors = dynamic_colors,
    dendrogram = gene_tree,
    dissimilarity = dissimilarity
  )
}

# ============================================================================
# Module Analysis
# ============================================================================

#' Calculate module-trait associations
#'
#' @param expression_data Expression matrix
#' @param module_colors Module color assignments
#' @param traits Trait data frame
#' @param correlation_method Correlation method
#' @param verbose Print progress
#'
#' @return Data frame of module-trait correlations and p-values
#' @keywords internal
#' @noRd
calculate_module_trait_correlations <- function(expression_data, module_colors,
                                                 traits, correlation_method = "pearson",
                                                 verbose = TRUE) {
  if (verbose) {
    cat("Calculating module-trait associations...\n")
  }

  # Calculate module eigengenes
  module_eigengenes <- WGCNA::moduleEigengenes(
    expr = expression_data,
    colors = module_colors,
    corFnc = correlation_method
  )$eigengenes

  # Calculate correlations
  module_trait_cor <- stats::cor(module_eigengenes, traits, use = "complete.obs",
                                 method = correlation_method)
  module_trait_pval <- matrix(NA, nrow = nrow(module_trait_cor),
                              ncol = ncol(module_trait_cor))

  for (i in 1:nrow(module_trait_cor)) {
    for (j in 1:ncol(module_trait_cor)) {
      pval <- stats::cor.test(module_eigengenes[, i], traits[, j],
                             method = correlation_method)$p.value
      module_trait_pval[i, j] <- pval
    }
  }

  rownames(module_trait_cor) <- colnames(module_eigengenes)
  rownames(module_trait_pval) <- colnames(module_eigengenes)
  colnames(module_trait_cor) <- colnames(traits)
  colnames(module_trait_pval) <- colnames(traits)

  if (verbose) {
    cat("Module-trait correlations calculated.\n\n")
  }

  list(
    correlations = module_trait_cor,
    pvalues = module_trait_pval,
    eigengenes = module_eigengenes
  )
}

# ============================================================================
# Output Export
# ============================================================================

#' Save WGCNA results to files
#'
#' @param results WGCNA results list
#' @param output_dir Output directory path
#' @param save_network Save network files
#' @param verbose Print progress
#'
#' @return Invisibly returns TRUE
#' @keywords internal
#' @noRd
export_wgcna_results <- function(results, output_dir, save_network = TRUE,
                                verbose = TRUE) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  if (verbose) {
    cat("Exporting WGCNA results to:", output_dir, "\n")
  }

  # Save module assignments
  if (!is.null(results$module_colors)) {
    module_df <- data.frame(
      gene = names(results$module_colors),
      module = results$module_colors
    )
    utils::write.csv(module_df,
      file.path(output_dir, "module_assignments.csv"),
      row.names = FALSE)
  }

  # Save module-trait correlations
  if (!is.null(results$module_trait_correlations)) {
    utils::write.csv(results$module_trait_correlations,
      file.path(output_dir, "module_trait_correlations.csv"))
  }

  if (verbose) {
    cat("Results exported successfully.\n\n")
  }

  invisible(TRUE)
}
