#' Intelligently reorder data for optimal heatmap visualization
#'
#' @description
#' Reorder rows and columns of data in long format to optimize heatmap
#' visualization by grouping features based on their peak expression patterns.
#' This function identifies where each feature reaches its maximum value and
#' groups similar patterns together to reveal biological clusters and trends.
#'
#' @param data A data frame in long format (tidy format) containing the data
#'   to be reordered. Must contain columns specified by \code{col}, \code{row},
#'   and \code{value} parameters.
#' @param col Character string specifying the column name containing sample
#'   identifiers (typically samples, conditions, or time points). This will
#'   become the columns in the heatmap.
#' @param row Character string specifying the column name containing feature
#'   identifiers (typically genes, metabolites, or other measured entities).
#'   This will become the rows in the heatmap.
#' @param value Character string specifying the column name containing the
#'   numerical values to be displayed in the heatmap (e.g., expression levels,
#'   abundances, intensities).
#' @param method Character string specifying the reordering method. Options:
#'   \itemize{
#'     \item "peak" (default): Order by peak expression patterns
#'     \item "hierarchical": Use hierarchical clustering
#'     \item "mean": Order by mean values
#'   }
#' @param ties_method Character string specifying how to handle ties when
#'   finding peak values. Options: "first", "last", "random". Default is "first".
#' @param na_handling Character string specifying how to handle missing values.
#'   Options: "remove", "zero", "mean". Default is "remove".
#' @param verbose Logical indicating whether to print progress information.
#'   Default is TRUE.
#'
#' @return A list containing four components:
#' \describe{
#'   \item{data.reordered}{A data frame with the same structure as input but
#'     with reordered factor levels for optimal heatmap visualization.}
#'   \item{column.order}{A character vector showing the optimal column order
#'     based on the number of features that peak in each column.}
#'   \item{row.order}{A character vector showing the optimal row order based
#'     on peak expression patterns.}
#'   \item{peak.summary}{A summary data frame showing where each feature peaks
#'     and how many features peak in each column.}
#' }
#'
#' @details
#' The reordering algorithm works by:
#' \enumerate{
#'   \item Identifying where each feature (row) reaches its maximum value
#'   \item Counting how many features peak in each sample (column)
#'   \item Ordering columns by the number of features that peak there
#'   \item Ordering rows by their peak expression patterns
#' }
#'
#' This approach reveals:
#' \itemize{
#'   \item Co-expression patterns and modules
#'   \item Temporal expression cascades in time-series data
#'   \item Treatment-specific responses
#'   \item Functional groupings of features
#' }
#'
#' @note
#' \itemize{
#'   \item Input data must be in long (tidy) format, not wide format
#'   \item The function assumes higher values represent higher expression/activity
#'   \item For count data, consider log-transformation before reordering
#'   \item Missing values are handled according to the na_handling parameter
#' }
#'
#' @export
#'
#' @examples
#' library(dplyr)
#' library(bioRtools)
#'
#' # Load example metabolomics data
#' data(df.reorder2heatmap)
#'
#' # Basic reordering for heatmap visualization
#' reorder_result <- reorder_heatmap(
#'   data = df.reorder2heatmap,
#'   col = "sample",
#'   row = "meta",
#'   value = "value"
#' )
#'
#' # View reordered data
#' head(reorder_result$data.reordered)
#'
#' # Check peak summary
#' reorder_result$peak.summary
#'
#' # View column and row orders
#' reorder_result$column.order
#' head(reorder_result$row.order)
#'
#' # Use the reordered data for heatmap
#' library(ggplot2)
#' reorder_result$data.reordered %>%
#'   ggplot(aes(x = col, y = row, fill = value)) +
#'   geom_tile() +
#'   scale_fill_viridis_c() +
#'   theme_minimal() +
#'   theme(axis.text.x = element_text(angle = 45, hjust = 1))
#'
#' # Alternative reordering methods
#' reorder_hierarchical <- reorder_heatmap(
#'   data = df.reorder2heatmap,
#'   col = "sample",
#'   row = "meta",
#'   value = "value",
#'   method = "hierarchical"
#' )
#'
#' # Create example time-series data
#' time_series_data <- data.frame(
#'   gene = rep(paste0("Gene_", 1:20), each = 6),
#'   timepoint = rep(paste0("T", 1:6), 20),
#'   expression = c(
#'     # Early genes
#'     rep(c(10, 8, 4, 2, 1, 1), 5),
#'     # Mid genes
#'     rep(c(2, 6, 10, 8, 3, 1), 5),
#'     # Late genes
#'     rep(c(1, 2, 3, 6, 10, 8), 5),
#'     # Constitutive genes
#'     rep(c(5, 5, 5, 5, 5, 5), 5)
#'   ) + rnorm(120, 0, 0.5)
#' )
#'
#' # Reorder time-series data
#' time_reordered <- reorder_heatmap(
#'   data = time_series_data,
#'   col = "timepoint",
#'   row = "gene",
#'   value = "expression"
#' )
#'
reorder_heatmap <- function(data,
                            col,
                            row,
                            value,
                            method = "peak",
                            ties_method = "first",
                            na_handling = "remove",
                            verbose = TRUE) {
  # Input validation
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }

  # Check if required columns exist
  required_cols <- c(col, row, value)
  missing_cols <- required_cols[!required_cols %in% colnames(data)]

  if (length(missing_cols) > 0) {
    stop("Required columns missing from data: ", paste(missing_cols, collapse = ", "))
  }

  # Validate parameters
  valid_methods <- c("peak", "hierarchical", "mean")
  if (!method %in% valid_methods) {
    stop("'method' must be one of: ", paste(valid_methods, collapse = ", "))
  }

  valid_ties <- c("first", "last", "random")
  if (!ties_method %in% valid_ties) {
    stop("'ties_method' must be one of: ", paste(valid_ties, collapse = ", "))
  }

  valid_na <- c("remove", "zero", "mean")
  if (!na_handling %in% valid_na) {
    stop("'na_handling' must be one of: ", paste(valid_na, collapse = ", "))
  }

  if (verbose) {
    message("Starting heatmap reordering with method: ", method)
    message("Input data: ", nrow(data), " rows, ",
      length(unique(data[[col]])), " columns, ",
      length(unique(data[[row]])), " features")
  }

  # Standardize column names for internal processing
  data_processed <- data %>%
    dplyr::select(
      col_var = !!rlang::sym(col),
      row_var = !!rlang::sym(row),
      value_var = !!rlang::sym(value)
    )

  # Handle missing values
  original_rows <- nrow(data_processed)

  if (na_handling == "remove") {
    data_processed <- data_processed %>%
      dplyr::filter(!is.na(value_var))
  } else if (na_handling == "zero") {
    data_processed <- data_processed %>%
      dplyr::mutate(value_var = ifelse(is.na(value_var), 0, value_var))
  } else if (na_handling == "mean") {
    mean_value <- mean(data_processed$value_var, na.rm = TRUE)
    data_processed <- data_processed %>%
      dplyr::mutate(value_var = ifelse(is.na(value_var), mean_value, value_var))
  }

  if (verbose && nrow(data_processed) < original_rows) {
    message("Removed ", original_rows - nrow(data_processed), " rows with missing values")
  }

  # Check if data remains
  if (nrow(data_processed) == 0) {
    stop("No data remaining after handling missing values")
  }

  # Apply reordering method
  if (method == "peak") {
    result <- reorder_by_peak(data_processed, ties_method, verbose)
  } else if (method == "hierarchical") {
    result <- reorder_by_clustering(data_processed, verbose)
  } else if (method == "mean") {
    result <- reorder_by_mean(data_processed, verbose)
  }

  # Restore original column names
  final_data <- result$data %>%
    dplyr::rename(
      !!col := col_var,
      !!row := row_var,
      !!value := value_var
    )

  if (verbose) {
    message("Reordering completed successfully!")
  }

  # Return comprehensive results following project conventions
  return(list(
    data.reordered = final_data,
    column.order = result$column_order,
    row.order = result$row_order,
    peak.summary = result$summary
  ))
}

# Helper function for peak-based reordering
reorder_by_peak <- function(data, ties_method, verbose) {
  # Find peak values for each feature
  peak_data <- data %>%
    dplyr::group_by(row_var) %>%
    dplyr::slice_max(
      order_by = value_var,
      n = 1,
      with_ties = ifelse(ties_method == "first", FALSE, TRUE)
    ) %>%
    dplyr::ungroup()

  # Handle ties if method is not "first"
  if (ties_method == "last") {
    peak_data <- peak_data %>%
      dplyr::group_by(row_var) %>%
      dplyr::slice_tail(n = 1) %>%
      dplyr::ungroup()
  } else if (ties_method == "random") {
    peak_data <- peak_data %>%
      dplyr::group_by(row_var) %>%
      dplyr::slice_sample(n = 1) %>%
      dplyr::ungroup()
  }

  # Count peaks per column
  peak_counts <- peak_data %>%
    dplyr::count(col_var, name = "peak_count") %>%
    dplyr::arrange(dplyr::desc(peak_count))

  # Create column order
  column_order <- peak_counts$col_var

  # Create summary
  peak_summary <- peak_data %>%
    dplyr::select(row_var, peaked_at = col_var) %>%
    dplyr::left_join(peak_counts, by = c("peaked_at" = "col_var"))

  # Reorder data
  reordered_data <- data %>%
    dplyr::left_join(peak_summary, by = "row_var") %>%
    dplyr::mutate(
      col_var = factor(col_var, levels = column_order),
      peaked_at = factor(peaked_at, levels = column_order),
      row_order = as.numeric(peaked_at),
      row_var = stats::reorder(row_var, row_order)
    ) %>%
    dplyr::select(col_var, row_var, value_var) %>%
    dplyr::arrange(col_var, row_var)

  return(list(
    data = reordered_data,
    column_order = as.character(column_order),
    row_order = levels(reordered_data$row_var),
    summary = peak_summary %>% dplyr::select(-row_order)
  ))
}

# Helper function for hierarchical clustering
reorder_by_clustering <- function(data, verbose) {
  # Convert to wide format for clustering
  wide_data <- data %>%
    tidyr::pivot_wider(
      names_from = col_var,
      values_from = value_var,
      values_fill = 0
    ) %>%
    tibble::column_to_rownames("row_var")

  # Cluster rows
  if (nrow(wide_data) > 1) {
    row_dist <- dist(wide_data)
    row_clust <- hclust(row_dist)
    row_order <- rownames(wide_data)[row_clust$order]
  } else {
    row_order <- rownames(wide_data)
  }

  # Cluster columns
  if (ncol(wide_data) > 1) {
    col_dist <- dist(t(wide_data))
    col_clust <- hclust(col_dist)
    column_order <- colnames(wide_data)[col_clust$order]
  } else {
    column_order <- colnames(wide_data)
  }

  # Reorder data
  reordered_data <- data %>%
    dplyr::mutate(
      col_var = factor(col_var, levels = column_order),
      row_var = factor(row_var, levels = row_order)
    ) %>%
    dplyr::arrange(col_var, row_var)

  # Create summary
  summary_data <- data.frame(
    method = "hierarchical",
    n_rows = length(row_order),
    n_cols = length(column_order)
  )

  return(list(
    data = reordered_data,
    column_order = column_order,
    row_order = row_order,
    summary = summary_data
  ))
}

# Helper function for mean-based reordering
reorder_by_mean <- function(data, verbose) {
  # Calculate mean values
  row_means <- data %>%
    dplyr::group_by(row_var) %>%
    dplyr::summarise(mean_value = mean(value_var, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean_value))

  col_means <- data %>%
    dplyr::group_by(col_var) %>%
    dplyr::summarise(mean_value = mean(value_var, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean_value))

  # Create orders
  row_order <- row_means$row_var
  column_order <- col_means$col_var

  # Reorder data
  reordered_data <- data %>%
    dplyr::mutate(
      col_var = factor(col_var, levels = column_order),
      row_var = factor(row_var, levels = row_order)
    ) %>%
    dplyr::arrange(col_var, row_var)

  # Create summary
  summary_data <- data.frame(
    method = "mean",
    row_means = "ordered by descending mean",
    col_means = "ordered by descending mean"
  )

  return(list(
    data = reordered_data,
    column_order = as.character(column_order),
    row_order = as.character(row_order),
    summary = summary_data
  ))
}
