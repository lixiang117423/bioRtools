#' Extract hierarchical information from phylogenetic tree
#'
#' @description
#' This function extracts hierarchical node information from a phylogenetic tree,
#' creating a comprehensive table showing the path from each tip to the root node.
#' The resulting data frame contains parent-child relationships at multiple levels,
#' making it useful for downstream phylogenetic analyses and taxonomic classifications.
#'
#' @param tree_file Character string specifying the path to the phylogenetic tree
#'   file. Supported formats include Newick (.nwk, .tre, .tree) and other formats
#'   supported by \code{ggtree::read.tree()}.
#' @param include_internal_nodes Logical indicating whether to include internal
#'   (non-tip) nodes in the output. Default is FALSE (tips only).
#' @param max_path_length Numeric value specifying the maximum path length to prevent
#'   infinite loops in case of circular references. Default is 1000.
#' @param column_prefix Character string specifying the prefix for hierarchy level
#'   columns. Default is "parent_".
#' @param return_format Character string specifying the output format. Options are
#'   "wide" (default) for wide format with separate columns for each level, or
#'   "long" for long format with path information in a single column.
#' @param verbose Logical indicating whether to print progress messages. Default is TRUE.
#'
#' @return A data frame containing hierarchical information. In wide format,
#'   each row represents a tip node with columns showing the node IDs at each
#'   hierarchical level from root to tip. In long format, each row represents
#'   a single node-level combination.
#'
#' @details
#' This function performs the following steps:
#' \itemize{
#'   \item Reads the phylogenetic tree from the specified file
#'   \item Extracts node-parent relationships using ggtree
#'   \item Traces the path from each node to the root
#'   \item Converts the hierarchical information to a structured data frame
#' }
#'
#' The function uses an internal helper function \code{find_path_to_root} to
#' trace the complete path from any given node to the root node. This ensures
#' that all hierarchical relationships are captured accurately.
#'
#' @note
#' \itemize{
#'   \item Large trees may take significant time to process
#'   \item The function assumes the tree is properly rooted
#'   \item Circular references in the tree will trigger a warning
#'   \item Missing node labels will be preserved as NA values
#' }
#'
#' @importFrom ggtree read.tree ggtree
#' @importFrom dplyr rowwise mutate ungroup select filter group_by
#' @importFrom stringr str_split
#' @importFrom tidyr unnest pivot_wider
#'
#' @export
#'
#' @examples
#' library(bioRtools)
#' library(writexl)
#' \dontrun{
#' # Basic usage - extract hierarchy information
#' tree_hierarchy <- extract_tree_hierarchy("path/to/your/tree.nwk")
#'
#' # Save to Excel file
#' writexl::write_xlsx(tree_hierarchy, "hierarchy_output.xlsx")
#'
#' # Include internal nodes and use long format
#' tree_long <- extract_tree_hierarchy(
#'   tree_file = "path/to/your/tree.nwk",
#'   include_internal_nodes = TRUE,
#'   return_format = "long"
#' )
#'
#' # Custom column prefix and maximum path length
#' tree_custom <- extract_tree_hierarchy(
#'   tree_file = "path/to/your/tree.nwk",
#'   column_prefix = "level_",
#'   max_path_length = 500,
#'   verbose = FALSE
#' )
#' }
#'
extract_tree_hierarchy <- function(tree_file,
                                   include_internal_nodes = FALSE,
                                   max_path_length = 1000,
                                   column_prefix = "parent_",
                                   return_format = c("wide", "long"),
                                   verbose = TRUE) {
  # Input validation
  if (!is.character(tree_file) || length(tree_file) != 1) {
    stop("'tree_file' must be a single character string")
  }

  if (!file.exists(tree_file)) {
    stop(paste("Tree file does not exist:", tree_file))
  }

  if (!is.logical(include_internal_nodes) || length(include_internal_nodes) != 1) {
    stop("'include_internal_nodes' must be a single logical value")
  }

  if (!is.numeric(max_path_length) || max_path_length <= 0) {
    stop("'max_path_length' must be a positive number")
  }

  if (!is.character(column_prefix) || length(column_prefix) != 1) {
    stop("'column_prefix' must be a single character string")
  }

  return_format <- match.arg(return_format)

  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("'verbose' must be a single logical value")
  }

  # Check required packages
  required_packages <- c("ggtree", "dplyr", "stringr", "tidyr")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

  if (length(missing_packages) > 0) {
    stop(paste("Required packages not available:", paste(missing_packages, collapse = ", ")))
  }

  # Progress message
  if (verbose) {
    message("Reading phylogenetic tree from: ", tree_file)
  }

  # Read tree and extract data
  tryCatch(
    {
      tree <- ggtree::read.tree(tree_file) %>%
        ggtree::ggtree()

      df <- tree$data
    },
    error = function(e) {
      stop(paste("Failed to read tree file:", e$message))
    })

  if (verbose) {
    message("Tree loaded successfully. Processing ", nrow(df), " nodes...")
  }

  # Internal function to find path to root
  find_path_to_root <- function(node_id, parent_child_map, max_length = max_path_length) {
    path <- c()
    current_node <- node_id

    # Trace path until reaching root node
    while (TRUE) {
      path <- c(path, current_node)

      # Find parent node
      parent_node <- parent_child_map[[as.character(current_node)]]

      # Check if reached root or if parent equals current node
      if (is.null(parent_node) || parent_node == current_node) {
        break
      }

      current_node <- parent_node

      # Prevent infinite loops
      if (length(path) > max_length) {
        warning(paste("Possible circular reference detected for node:", node_id))
        break
      }
    }

    return(path)
  }

  # Create node-to-parent mapping
  parent_child_map <- setNames(df$parent, df$node)

  if (verbose) {
    message("Calculating hierarchical paths...")
  }

  # Calculate paths for each node
  df_processed <- df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      path_to_root = list(find_path_to_root(node, parent_child_map)),
      path_length = length(path_to_root),
      path_string = paste(rev(path_to_root), collapse = " -> ")
    ) %>%
    dplyr::ungroup()

  # Filter nodes based on include_internal_nodes parameter
  if (!include_internal_nodes) {
    df_processed <- df_processed %>%
      dplyr::filter(isTip == TRUE)
  }

  if (verbose) {
    message("Converting to structured format...")
  }

  # Process based on return format
  if (return_format == "wide") {
    # Convert to wide format
    df_final <- df_processed %>%
      dplyr::select(label, isTip, path_string) %>%
      dplyr::mutate(path = stringr::str_split(path_string, " -> ")) %>%
      dplyr::select(-path_string) %>%
      tidyr::unnest(path) %>%
      dplyr::group_by(label) %>%
      dplyr::mutate(levels = paste0(column_prefix, 1:dplyr::n())) %>%
      dplyr::ungroup() %>%
      tidyr::pivot_wider(names_from = levels, values_from = path) %>%
      dplyr::select(dplyr::starts_with(column_prefix), dplyr::everything()) %>%
      dplyr::filter(!is.na(label))

    if (!include_internal_nodes) {
      df_final <- df_final %>%
        dplyr::select(-isTip)
    }
  } else {
    # Return long format
    df_final <- df_processed %>%
      dplyr::select(label, isTip, path_string, path_length) %>%
      dplyr::filter(!is.na(label))

    if (!include_internal_nodes) {
      df_final <- df_final %>%
        dplyr::select(-isTip)
    }
  }

  if (verbose) {
    message("Processing completed. Returning ", nrow(df_final), " rows.")
  }

  return(df_final)
}
