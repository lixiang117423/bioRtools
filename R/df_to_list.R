#' Convert Data Frame to Named List for Set Analysis and Venn Diagrams
#'
#' This function converts a data frame with group-value pairs into a named list
#' where each element contains all values belonging to a specific group. This
#' format is commonly used for Venn diagram plotting, set analysis, and overlap
#' studies in biological data analysis.
#'
#' @param data A data frame containing at least two columns: one for grouping
#'   and one for values
#' @param group Column name (as string) that contains the group identifiers.
#'   This will become the names of the list elements
#' @param value Column name (as string) that contains the values to be grouped.
#'   These values will be collected into vectors for each group
#'
#' @return A named list where:
#'   \itemize{
#'     \item Each element name corresponds to a unique group from the group column
#'     \item Each element contains a vector of all values belonging to that group
#'     \item Missing values (NA) are automatically removed from each group
#'     \item Empty groups (after NA removal) are retained as empty vectors
#'   }
#'
#' @details
#' This function is particularly useful for preparing data for:
#' \itemize{
#'   \item Venn diagram creation (using packages like VennDiagram, ggVennDiagram)
#'   \item Set overlap analysis
#'   \item Gene/feature intersection studies
#'   \item Comparative analysis between different experimental conditions
#' }
#'
#' The function handles missing values by removing them from the final lists,
#' ensuring clean input for downstream analysis tools. The order of groups
#' in the resulting list follows the order of first appearance in the input data.
#'
#' @note
#' \itemize{
#'   \item The input data frame should contain unique group-value combinations
#'   \item Duplicate values within the same group will be preserved
#'   \item Missing values in either group or value columns are handled gracefully
#' }
#'
#' @export
#'
#' @examples
#' library(bioRtools)
#' library(dplyr)
#'
#' # Example 1: Basic usage with gene lists
#' gene_data <- data.frame(
#'   condition = c(rep("Treatment_A", 5), rep("Treatment_B", 5), rep("Control", 4)),
#'   gene_id = c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5",
#'     "GENE3", "GENE4", "GENE6", "GENE7", "GENE8",
#'     "GENE1", "GENE9", "GENE10", "GENE11")
#' )
#'
#' gene_list <- df_to_list(
#'   data = gene_data,
#'   group = "condition",
#'   value = "gene_id"
#' )
#' print(gene_list)
#' print(lengths(gene_list))  # Number of genes per condition
#'
#' # Example 2: Creating overlapping sets (original example improved)
#' df_A <- data.frame(group = "A", value = 1:10)
#' df_B <- data.frame(group = "B", value = 3:12)
#' df_C <- data.frame(group = "C", value = 7:20)
#'
#' combined_data <- dplyr::bind_rows(df_A, df_B, df_C)
#'
#' overlap_sets <- df_to_list(
#'   data = combined_data,
#'   group = "group",
#'   value = "value"
#' )
#' print(overlap_sets)
#'
#' # Find intersections
#' intersection_AB <- intersect(overlap_sets$A, overlap_sets$B)
#' intersection_BC <- intersect(overlap_sets$B, overlap_sets$C)
#' intersection_AC <- intersect(overlap_sets$A, overlap_sets$C)
#' intersection_ABC <- intersect(intersection_AB, overlap_sets$C)
#'
#' print(paste("A ∩ B:", paste(intersection_AB, collapse = ", ")))
#' print(paste("B ∩ C:", paste(intersection_BC, collapse = ", ")))
#' print(paste("A ∩ C:", paste(intersection_AC, collapse = ", ")))
#' print(paste("A ∩ B ∩ C:", paste(intersection_ABC, collapse = ", ")))
#'
#' # Example 3: Real biological application - DEG analysis
#' deg_results <- data.frame(
#'   comparison = c(rep("TreatmentA_vs_Control", 6),
#'     rep("TreatmentB_vs_Control", 5),
#'     rep("TreatmentA_vs_TreatmentB", 4)),
#'   gene = c("ATG1", "ATG2", "PTEN", "TP53", "MYC", "EGFR",
#'     "ATG1", "PTEN", "BRCA1", "BRCA2", "MYC",
#'     "TP53", "MYC", "EGFR", "KRAS")
#' )
#'
#' deg_lists <- df_to_list(
#'   data = deg_results,
#'   group = "comparison",
#'   value = "gene"
#' )
#'
#' print("DEG lists for Venn diagram:")
#' print(deg_lists)
#'
#' # Example 4: Handling missing values
#' data_with_na <- data.frame(
#'   category = c("Group1", "Group1", "Group2", "Group2", "Group1"),
#'   item = c("A", NA, "B", "C", "D")
#' )
#'
#' clean_list <- df_to_list(
#'   data = data_with_na,
#'   group = "category",
#'   value = "item"
#' )
#' print("Data with NA values handled:")
#' print(clean_list)
#'
df_to_list <- function(data, group, value) {
  # Input validation
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }

  if (nrow(data) == 0) {
    stop("'data' cannot be empty")
  }

  if (!is.character(group) || length(group) != 1) {
    stop("'group' must be a single character string specifying the column name")
  }

  if (!is.character(value) || length(value) != 1) {
    stop("'value' must be a single character string specifying the column name")
  }

  if (!group %in% names(data)) {
    stop(paste("Column '", group, "' not found in data"))
  }

  if (!value %in% names(data)) {
    stop(paste("Column '", value, "' not found in data"))
  }

  # Check for completely missing group or value columns
  if (all(is.na(data[[group]]))) {
    stop("The group column contains only missing values")
  }

  # Warn about missing values that will be handled
  na_groups <- sum(is.na(data[[group]]))
  na_values <- sum(is.na(data[[value]]))

  if (na_groups > 0) {
    warning(paste("Found", na_groups, "missing values in group column. These rows will be excluded."))
  }

  if (na_values > 0) {
    warning(paste("Found", na_values, "missing values in value column. These will be excluded from groups."))
  }

  # Clean and prepare data
  clean_data <- data %>%
    dplyr::select(!!group, !!value) %>%
    # Remove rows where group is NA (can't assign to a group)
    dplyr::filter(!is.na(!!rlang::sym(group))) %>%
    # Create a row identifier for pivot operations
    dplyr::mutate(row_id = dplyr::row_number()) %>%
    # Rename columns to standard names for easier processing
    dplyr::rename(group_col = !!group, value_col = !!value)

  if (nrow(clean_data) == 0) {
    warning("No valid data remaining after removing missing group values")
    return(list())
  }

  # Get unique groups in order of first appearance
  unique_groups <- unique(clean_data$group_col)

  # Create the list using split function (more efficient than pivot)
  result_list <- clean_data %>%
    # Remove rows with missing values (but keep the group)
    dplyr::filter(!is.na(value_col)) %>%
    # Split values by group
    dplyr::pull(value_col) %>%
    split(clean_data$group_col[!is.na(clean_data$value_col)])

  # Ensure all groups are represented (even if empty after NA removal)
  missing_groups <- setdiff(unique_groups, names(result_list))
  for (grp in missing_groups) {
    result_list[[grp]] <- character(0)
  }

  # Reorder list to match original group order
  result_list <- result_list[unique_groups]

  # Add attributes for metadata
  attr(result_list, "source_data_rows") <- nrow(data)
  attr(result_list, "processed_rows") <- nrow(clean_data)
  attr(result_list, "group_column") <- group
  attr(result_list, "value_column") <- value
  attr(result_list, "total_unique_values") <- length(unique(unlist(result_list)))

  # Print summary in interactive mode
  if (interactive()) {
    cat("Data frame to list conversion completed:\n")
    cat("Groups found:", length(result_list), "\n")
    cat("Group names:", paste(names(result_list), collapse = ", "), "\n")
    cat("Group sizes:", paste(lengths(result_list), collapse = ", "), "\n")
    cat("Total unique values:", attr(result_list, "total_unique_values"), "\n")
    if (na_groups > 0 || na_values > 0) {
      cat("Note: Missing values were handled during conversion\n")
    }
    cat("\n")
  }

  return(result_list)
}
