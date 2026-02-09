#' Perform rarefaction to standardize sampling depth across samples
#'
#' @description
#' Perform rarefaction analysis to standardize the sampling depth across all
#' samples in a community data matrix. This is essential for comparing diversity
#' metrics across samples with different sequencing depths in microbiome studies.
#'
#' @param data Community data matrix where rows represent samples and columns
#'   represent species/OTUs/taxa. The data should contain count data (integers).
#' @param value Subsample size for rarefying community. Can be:
#'   \itemize{
#'     \item NULL (default): Uses the minimum row sum across all samples
#'     \item Single numeric value: Same rarefaction depth for all samples
#'     \item Numeric vector: Different rarefaction depth for each sample
#'   }
#' @param seed Integer value for setting random seed to ensure reproducible
#'   results. Default is 123.
#' @param min_reads Minimum number of reads required for a sample to be included
#'   in rarefaction. Samples with fewer reads will be removed. Default is 1000.
#' @param remove_empty Logical indicating whether to remove taxa (columns) that
#'   become absent after rarefaction. Default is TRUE.
#' @param verbose Logical indicating whether to print progress information.
#'   Default is TRUE.
#'
#' @return A list containing four components:
#' \describe{
#'   \item{data.rarefied}{A data frame containing the rarefied community matrix
#'     with the same structure as input but standardized sampling depth.}
#'   \item{rarefaction.info}{A summary data frame with original and rarefied
#'     read counts, number of taxa, and other statistics for each sample.}
#'   \item{removed.samples}{Character vector of sample names that were removed
#'     due to insufficient reads (if any).}
#'   \item{rarefaction.depth}{The actual rarefaction depth used for analysis.}
#' }
#'
#' @details
#' Rarefaction randomly subsamples reads from each sample to a standardized depth,
#' enabling fair comparison of diversity metrics across samples. The function:
#' \itemize{
#'   \item Removes samples with insufficient reads (below min_reads)
#'   \item Performs random subsampling using vegan::rrarefy()
#'   \item Optionally removes taxa that become absent after rarefaction
#'   \item Provides comprehensive summary statistics
#' }
#'
#' @note
#' \itemize{
#'   \item Input data should be raw count data (not relative abundances)
#'   \item Rarefaction involves random sampling, so results may vary between runs
#'     unless seed is set
#'   \item Consider the trade-off between retaining samples and rarefaction depth
#'   \item Some studies prefer alternative normalization methods to rarefaction
#' }
#'
#' @export
#'
#' @examples
#' library(dplyr)
#' library(bioRtools)
#'
#' # Load example microbiome data
#' data(df.permanova.otu)
#'
#' # Basic rarefaction using minimum sample depth
#' rarefy_result <- rarefy_table(data = df.permanova.otu)
#'
#' # View rarefied data
#' head(rarefy_result$data.rarefied)
#'
#' # Check rarefaction summary
#' rarefy_result$rarefaction.info
#'
#' # Custom rarefaction depth
#' rarefy_custom <- rarefy_table(
#'   data = df.permanova.otu,
#'   value = 5000,
#'   seed = 456
#' )
#'
#' # Check rarefaction depth used
#' rarefy_custom$rarefaction.depth
#'
#' # Multiple rarefaction depths (if you have multiple sample types)
#' # rarefy_multi <- rarefy_table(
#' #   data = df.permanova.otu,
#' #   value = c(rep(3000, 5), rep(4000, 4))  # Different depths per sample
#' # )
#'
#' # Transpose data if needed (samples as columns, taxa as rows)
#' df.permanova.otu %>%
#'   t() %>%
#'   as.data.frame() %>%
#'   rarefy_table(value = 3000) -> df.transposed.rarefied
#'
rarefy_table <- function(data,
                         value = NULL,
                         seed = 123,
                         min_reads = 1000,
                         remove_empty = TRUE,
                         verbose = TRUE) {
  # Input validation
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("'data' must be a data frame or matrix")
  }

  # Convert to data frame if matrix
  if (is.matrix(data)) {
    data <- as.data.frame(data)
  }

  # Check if data contains only numeric values
  if (!all(sapply(data, is.numeric))) {
    stop("All columns in 'data' must be numeric")
  }

  # Check for negative values
  if (any(data < 0, na.rm = TRUE)) {
    stop("Data contains negative values. Count data should be non-negative integers")
  }

  # Remove any samples with all zeros
  data_clean <- data[rowSums(data, na.rm = TRUE) > 0, ]

  if (nrow(data_clean) == 0) {
    stop("No samples contain any reads")
  }

  if (nrow(data_clean) < nrow(data)) {
    removed_empty <- rownames(data)[!rownames(data) %in% rownames(data_clean)]
    if (verbose) {
      message("Removed ", length(removed_empty), " samples with zero reads: ",
        paste(removed_empty, collapse = ", "))
    }
  }

  # Calculate sample read counts
  sample_sums <- rowSums(data_clean, na.rm = TRUE)

  # Remove samples below minimum threshold
  sufficient_samples <- sample_sums >= min_reads
  removed_samples <- rownames(data_clean)[!sufficient_samples]

  if (sum(sufficient_samples) == 0) {
    stop("No samples have sufficient reads (>= ", min_reads, ")")
  }

  if (length(removed_samples) > 0 && verbose) {
    message("Removed ", length(removed_samples), " samples with < ", min_reads,
      " reads: ", paste(removed_samples, collapse = ", "))
  }

  data_filtered <- data_clean[sufficient_samples, ]
  sample_sums_filtered <- sample_sums[sufficient_samples]

  # Determine rarefaction depth
  if (is.null(value)) {
    rarefaction_depth <- min(sample_sums_filtered)
    if (verbose) {
      message("Using minimum sample depth for rarefaction: ", rarefaction_depth)
    }
  } else {
    if (length(value) == 1) {
      rarefaction_depth <- value
      if (verbose) {
        message("Using custom rarefaction depth: ", rarefaction_depth)
      }

      # Check if any samples have fewer reads than rarefaction depth
      insufficient <- sample_sums_filtered < rarefaction_depth
      if (any(insufficient)) {
        warning("Some samples have fewer reads than the rarefaction depth. ",
          "These samples will be randomly subsampled with replacement.")
      }
    } else if (length(value) == nrow(data_filtered)) {
      rarefaction_depth <- value
      if (verbose) {
        message("Using individual rarefaction depths for each sample")
      }
    } else {
      stop("'value' must be NULL, a single number, or a vector with length equal to number of samples")
    }
  }

  # Set seed for reproducibility
  set.seed(seed)

  # Perform rarefaction
  if (verbose) {
    message("Performing rarefaction...")
  }

  tryCatch(
    {
      if (length(rarefaction_depth) == 1) {
        rarefied_data <- vegan::rrarefy(data_filtered, sample = rarefaction_depth)
      } else {
        rarefied_data <- vegan::rrarefy(data_filtered, sample = rarefaction_depth)
      }
    },
    error = function(e) {
      stop("Rarefaction failed: ", e$message)
    })

  # Convert to data frame
  rarefied_df <- as.data.frame(rarefied_data)

  # Remove empty taxa if requested
  if (remove_empty) {
    taxa_sums <- colSums(rarefied_df)
    non_empty_taxa <- taxa_sums > 0

    if (sum(non_empty_taxa) < ncol(rarefied_df)) {
      removed_taxa_count <- sum(!non_empty_taxa)
      rarefied_df <- rarefied_df[, non_empty_taxa]

      if (verbose) {
        message("Removed ", removed_taxa_count, " taxa that became absent after rarefaction")
      }
    }
  }

  # Create summary information
  rarefaction_info <- data.frame(
    sample = rownames(rarefied_df),
    original_reads = sample_sums_filtered,
    rarefied_reads = rowSums(rarefied_df),
    original_taxa = apply(data_filtered, 1, function(x) sum(x > 0)),
    rarefied_taxa = apply(rarefied_df, 1, function(x) sum(x > 0)),
    rarefaction_depth = if (length(rarefaction_depth) == 1) {
      rep(rarefaction_depth, nrow(rarefied_df))
    } else {
      rarefaction_depth
    },
    stringsAsFactors = FALSE
  )

  # Calculate additional statistics
  rarefaction_info$taxa_lost <- rarefaction_info$original_taxa - rarefaction_info$rarefied_taxa
  rarefaction_info$taxa_retention_rate <- round(
    rarefaction_info$rarefied_taxa / rarefaction_info$original_taxa * 100, 2
  )

  if (verbose) {
    message("Rarefaction completed successfully!")
    message("Final dataset: ", nrow(rarefied_df), " samples × ", ncol(rarefied_df), " taxa")
    if (length(unique(rarefaction_info$rarefaction_depth)) == 1) {
      message("Rarefaction depth: ", unique(rarefaction_info$rarefaction_depth))
    }
  }

  # Return comprehensive results following project conventions
  return(list(
    data.rarefied = rarefied_df,
    rarefaction.info = rarefaction_info,
    removed.samples = removed_samples,
    rarefaction.depth = if (length(rarefaction_depth) == 1) {
      rarefaction_depth
    } else {
      rarefaction_depth
    }
  ))
}
