#' Get Haplotypes from pheatmap Result
#'
#' Extract row haplotype clusters from a pheatmap result object based on hierarchical clustering tree.
#'
#' @param pheatmap_result A pheatmap result object
#' @param cluster_n Integer, number of haplotype clusters to cut the tree into
#'
#' @return A data frame with three columns:
#'   \itemize{
#'     \item \code{position}: The position (order) of the sample in the heatmap (1-indexed)
#'     \item \code{sample}: The sample ID
#'     \item \code{hap}: The haplotype ID (e.g., "Hap1", "Hap2")
#'   }
#'   The data frame is sorted by position in the heatmap order.
#' @export
#'
#' @examples
#' \dontrun{
#' # Create a simple heatmap
#' mat <- matrix(rnorm(100), nrow = 10)
#' rownames(mat) <- paste0("Gene", 1:10)
#' colnames(mat) <- paste0("Sample", 1:10)
#' ph_result <- pheatmap::pheatmap(mat)
#'
#' # Extract haplotypes
#' hap_df <- get_hap_from_heatmap(ph_result, cluster_n = 3)
#' # View results
#' print(hap_df)
#' # Filter specific haplotype
#' hap1_items <- hap_df[hap_df$hap == "Hap1", ]
#' }
#'
get_hap_from_heatmap <- function(pheatmap_result, cluster_n) {
  # Validate inputs
  if (!"pheatmap" %in% class(pheatmap_result)) {
    stop("Input must be a pheatmap result object")
  }

  if (!is.numeric(cluster_n) || length(cluster_n) != 1 || cluster_n < 1) {
    stop("cluster_n must be a positive integer")
  }

  cluster_n <- as.integer(cluster_n)

  # Extract row clustering tree from pheatmap result
  row_tree <- pheatmap_result$tree_row

  if (is.null(row_tree)) {
    stop("No row clustering tree found in pheatmap result. ",
         "Make sure to run pheatmap with cluster_rows = TRUE")
  }

  # Cut the tree into clusters
  cluster_assignments <- stats::cutree(row_tree, k = cluster_n)

  # Get row IDs (names) - use indices if no names
  row_ids <- names(cluster_assignments)
  if (is.null(row_ids) || all(row_ids == "")) {
    row_ids <- as.character(seq_along(cluster_assignments))
  }

  # Get the order of rows in the heatmap (1-indexed positions)
  heatmap_order <- row_tree$order

  # Create a mapping from ID to heatmap position
  # The order is based on indices, so we need to map correctly
  id_to_position <- setNames(
    as.character(seq_along(heatmap_order)),
    row_ids[heatmap_order]
  )

  # Actually, let's rethink this:
  # row_tree$order gives the indices in the original data that appear in order on heatmap
  # So position 1 on heatmap corresponds to row_ids[row_tree$order[1]]
  # We need to create: for each ID, what is its position on heatmap

  # Create data frame with sample and haplotype assignment
  cluster_df <- data.frame(
    sample = row_ids,
    hap = cluster_assignments,
    stringsAsFactors = FALSE
  )

  # Add heatmap position for each sample
  # heatmap_order[i] is the index of the row that appears at position i
  cluster_df$position <- NA_integer_
  for (pos in seq_along(heatmap_order)) {
    idx <- heatmap_order[pos]
    cluster_df$position[cluster_df$sample == row_ids[idx]] <- pos
  }

  # Add haplotype name column
  cluster_df$hap <- paste0("Hap", cluster_df$hap)

  # Sort by position only (to display in heatmap order)
  cluster_df <- cluster_df[order(cluster_df$position), ]

  # Keep only the columns we want, in the desired order
  result <- cluster_df[, c("position", "sample", "hap"), drop = FALSE]
  rownames(result) <- NULL

  # Print summary
  message(sprintf("\nSuccessfully extracted %d haplotype(s):", cluster_n))
  for (cl in sort(unique(cluster_df$hap))) {
    cl_data <- result[result$hap == cl, ]
    message(sprintf("  %s: %d items (positions %d - %d)",
                    cl,
                    nrow(cl_data),
                    min(cl_data$position),
                    max(cl_data$position)))
  }

  return(result)
}
