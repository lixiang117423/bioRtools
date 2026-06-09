#' Convert Data Frame to Named List for Set Analysis and Venn Diagrams
#'
#' This function converts a data frame with group-value pairs into a named list
#' where each element contains all values belonging to a specific group. This
#' format is commonly used for Venn diagram plotting, set analysis, and overlap
#' studies in biological data analysis.
#'
#' @param data A data frame containing at least two columns: one for grouping
#'   and one for values.
#' @param group Either a column name (as string) for group identifiers, or a
#'   formula like \code{value ~ group}. When a formula is passed, the
#'   \code{value} parameter is ignored.
#' @param value Column name (as string) for the values to be grouped.
#'   Ignored when \code{group} is a formula.
#' @param unique Logical. If TRUE (default), remove duplicate values within
#'   each group. Set to FALSE to keep duplicates.
#' @param verbose Logical. Print summary information. Default is TRUE.
#'
#' @return A named list where:
#'   \itemize{
#'     \item Each element name corresponds to a unique group
#'     \item Each element contains a vector of values belonging to that group
#'     \item Missing values (NA) are automatically removed
#'   }
#'
#' @details
#' Commonly used to prepare data for:
#' \itemize{
#'   \item Venn diagram creation (ggVennDiagram, VennDiagram)
#'   \item Set overlap and intersection analysis
#'   \item Gene/feature comparison between conditions
#' }
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' library(bioRtools)
#'
#' # Formula interface
#' gene_data <- data.frame(
#'   condition = c(rep("Treatment_A", 5), rep("Treatment_B", 5), rep("Control", 4)),
#'   gene_id = c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5",
#'     "GENE3", "GENE4", "GENE6", "GENE7", "GENE8",
#'     "GENE1", "GENE9", "GENE10", "GENE11")
#' )
#'
#' gene_list <- df2list(gene_data, gene_id ~ condition)
#' print(gene_list)
#'
#' # String interface (still supported)
#' df2list(gene_data, group = "condition", value = "gene_id")
#'
#' # Find intersections
#' intersect(gene_list$Treatment_A, gene_list$Treatment_B)
#'
df2list <- function(data, group, value = NULL, unique = TRUE, verbose = TRUE) {

  # Formula interface
  if (inherits(group, "formula")) {
    if (length(group) != 3) {
      stop("Formula must be in the form 'value ~ group'")
    }
    value <- as.character(group[[2]])
    group <- as.character(group[[3]])
  }

  # Input validation
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  if (nrow(data) == 0) {
    stop("'data' cannot be empty")
  }
  if (!group %in% names(data)) {
    stop("Column '", group, "' not found in data")
  }
  if (!value %in% names(data)) {
    stop("Column '", value, "' not found in data")
  }

  # Drop rows with NA in group or value
  df <- data[!is.na(data[[group]]) & !is.na(data[[value]]),
             c(group, value), drop = FALSE]

  if (nrow(df) == 0) {
    warning("No valid data remaining after removing missing values")
    return(list())
  }

  # Split by group
  result <- split(df[[value]], df[[group]])

  # Remove duplicates within each group
  if (unique) {
    result <- lapply(result, base::unique)
  }

  if (verbose) {
    message("Groups: ", length(result),
            " | Values per group: ", paste(lengths(result), collapse = ", "))
  }

  result
}
