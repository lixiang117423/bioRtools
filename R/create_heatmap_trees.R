#' Create Hierarchical Clustering Trees for Heatmap Annotation
#'
#' Generate row and column dendrograms from a data frame for use with
#' complex heatmap visualizations. This function performs hierarchical
#' clustering and returns ggtree objects that can be used to annotate
#' heatmaps with phylogenetic-style trees.
#'
#' @param data A data frame containing the data to be clustered
#' @param row_column Column name in `data` containing row identifiers
#' @param col_column Column name in `data` containing column identifiers
#' @param value_column Column name in `data` containing the values for clustering
#' @param distance_method Distance measure method (default: "euclidean").
#'   Options include "euclidean", "maximum", "manhattan", "canberra",
#'   "binary", "minkowski", "correlation", etc.
#' @param hclust_method Hierarchical clustering method (default: "complete").
#'   Options include "complete", "single", "average", "ward.D", "ward.D2",
#'   "mcquitty", "median", "centroid"
#'
#' @return A list containing:
#'   \code{tree_row}: ggtree object for row clustering
#'
#'   \code{tree_col}: ggtree object for column clustering
#'
#'   \code{row_hclust}: hclust object for row clustering
#'
#'   \code{col_hclust}: hclust object for column clustering
#'
#'   \code{row_order}: integer vector with row clustering order
#'
#'   \code{col_order}: integer vector with column clustering order
#'
#'   Each tree can be added to a plot using `+` operator or used with
#'   ggtree annotation functions.
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' \dontrun{
#' library(bioRtools)
#' library(dplyr)
#' library(ggplot2)
#'
#' # Create example data
#' set.seed(123)
#' example_data <- tibble::tibble(
#'   gene = rep(paste0("Gene", 1:10), each = 5),
#'   sample = rep(paste0("Sample", 1:5), 10),
#'   expression = rnorm(50, mean = 10, sd = 2)
#' )
#'
#' # Create clustering trees
#' trees <- create_heatmap_trees(
#'   data = example_data,
#'   row_column = "gene",
#'   col_column = "sample",
#'   value_column = "expression"
#' )
#'
#' # View row tree
#' print(trees$tree_row)
#'
#' # View column tree
#' print(trees$tree_col)
#'
#' # Use with different clustering methods
#' trees_correlation <- create_heatmap_trees(
#'   data = example_data,
#'   row_column = "gene",
#'   col_column = "sample",
#'   value_column = "expression",
#'   distance_method = "correlation",
#'   hclust_method = "average"
#' )
#'
#' print(trees_correlation$tree_row)
#' }
create_heatmap_trees <- function(data,
                                  row_column,
                                  col_column,
                                  value_column,
                                  distance_method = "euclidean",
                                  hclust_method = "complete") {

  # Input validation
  if (missing(data) || is.null(data)) {
    stop("'data' is required")
  }

  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }

  if (nrow(data) == 0) {
    stop("'data' is empty")
  }

  # Validate column names
  if (missing(row_column) || missing(col_column) || missing(value_column)) {
    stop("'row_column', 'col_column', and 'value_column' must all be specified")
  }

  if (!is.character(row_column) || length(row_column) != 1) {
    stop("'row_column' must be a single character string")
  }

  if (!is.character(col_column) || length(col_column) != 1) {
    stop("'col_column' must be a single character string")
  }

  if (!is.character(value_column) || length(value_column) != 1) {
    stop("'value_column' must be a single character string")
  }

  required_cols <- c(row_column, col_column, value_column)
  missing_cols <- setdiff(required_cols, names(data))

  if (length(missing_cols) > 0) {
    stop(paste(
      "Columns not found in 'data':",
      paste(missing_cols, collapse = ", ")
    ))
  }

  # Validate method parameters
  valid_distances <- c(
    "euclidean", "maximum", "manhattan", "canberra", "binary",
    "minkowski", "correlation", "spearman", "kendall"
  )

  if (!is.character(distance_method) || length(distance_method) != 1) {
    stop("'distance_method' must be a single character string")
  }

  if (!distance_method %in% valid_distances) {
    stop(paste(
      "'distance_method' must be one of:",
      paste(valid_distances, collapse = ", ")
    ))
  }

  valid_hclust <- c(
    "ward.D", "ward.D2", "single", "complete", "average",
    "mcquitty", "median", "centroid"
  )

  if (!is.character(hclust_method) || length(hclust_method) != 1) {
    stop("'hclust_method' must be a single character string")
  }

  if (!hclust_method %in% valid_hclust) {
    stop(paste(
      "'hclust_method' must be one of:",
      paste(valid_hclust, collapse = ", ")
    ))
  }

  # Reshape data to matrix format
  df_matrix <- data %>%
    dplyr::select({{ row_column }}, {{ col_column }}, {{ value_column }}) %>%
    tidyr::pivot_wider(
      id_cols = {{ row_column }},
      names_from = {{ col_column }},
      values_from = {{ value_column }}
    )

  # Convert row names to character
  row_names <- as.character(df_matrix[[1]])
  df_matrix <- df_matrix %>%
    tibble::column_to_rownames(var = names(df_matrix)[1])

  # Convert to numeric matrix
  tryCatch({
    df_matrix <- as.matrix(df_matrix)
    mode(df_matrix) <- "numeric"
  }, error = function(e) {
    stop("Failed to convert data to numeric matrix. Check that all values are numeric.")
  })

  # Check for missing values
  if (any(is.na(df_matrix))) {
    warning("Data contains NA values. These will be replaced with 0 for distance calculation.")
    df_matrix[is.na(df_matrix)] <- 0
  }

  # Special handling for correlation-based distance
  if (distance_method == "correlation" || distance_method == "spearman" ||
      distance_method == "kendall") {

    # Function to calculate correlation-based distance
    cor_distance <- function(x) {
      if (distance_method == "spearman") {
        cor_matrix <- cor(x, method = "spearman", use = "complete.obs")
      } else if (distance_method == "kendall") {
        cor_matrix <- cor(x, method = "kendall", use = "complete.obs")
      } else {
        cor_matrix <- cor(x, method = "pearson", use = "complete.obs")
      }
      as.dist(1 - abs(cor_matrix))
    }

    # Row clustering
    row_dist <- cor_distance(t(df_matrix))
    tree_row <- stats::hclust(row_dist, method = hclust_method)

    # Column clustering
    col_dist <- cor_distance(df_matrix)
    tree_col <- stats::hclust(col_dist, method = hclust_method)

  } else {
    # Standard distance-based clustering

    # Row clustering
    row_dist <- stats::dist(df_matrix, method = distance_method)
    tree_row <- stats::hclust(row_dist, method = hclust_method)

    # Column clustering
    col_dist <- stats::dist(t(df_matrix), method = distance_method)
    tree_col <- stats::hclust(col_dist, method = hclust_method)
  }

  # Convert to ggtree objects
  tree_row_gg <- ggtree::ggtree(tree_row, branch.length = "none") +
    ggtree::geom_tiplab(size = 3, align = TRUE)

  tree_col_gg <- ggtree::ggtree(tree_col, branch.length = "none") +
    ggtree::geom_tiplab(size = 3, align = TRUE)

  # Return results
  result <- list(
    tree_row = tree_row_gg,
    tree_col = tree_col_gg,
    row_hclust = tree_row,
    col_hclust = tree_col,
    row_order = tree_row$order,
    col_order = tree_col$order
  )

  result
}
