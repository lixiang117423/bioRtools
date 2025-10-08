#' Plot Motif Locations on Sequences
#'
#' @description
#' Visualize the location of motifs identified by MEME on sequences. Can
#' optionally display sequences aligned with a phylogenetic tree.
#'
#' @param data A data frame from \code{get_motif_from_meme()} containing motif
#'   information with columns: input_seq_id, length, motif_id, start_position,
#'   and end_position.
#' @param tree_path Optional character string. Path to a phylogenetic tree file
#'   in Newick format. Tree tip labels must match the sequence IDs used in MEME
#'   analysis. Default is NULL.
#' @param tree_annotation Optional data frame for annotating tree tips. Must
#'   contain columns "label" (matching tree tip labels) and "Group" (for coloring).
#'   Default is NULL.
#'
#' @return A ggplot2 object showing motif locations. If a tree is provided,
#'   returns a combined plot with the tree and motif locations side-by-side.
#'
#' @author Xiang LI \email{lixiang117423@@foxmail.com}
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select mutate rename filter arrange left_join
#' @importFrom ggplot2 ggplot aes geom_segment geom_rect scale_x_continuous
#'   labs theme theme_classic element_blank
#' @importFrom ggtree ggtree geom_tippoint
#' @importFrom ape read.tree
#' @importFrom grid unit
#' @importFrom patchwork plot_layout
#' @importFrom treeio as_tibble as.treedata
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' # Without phylogenetic tree
#' filepath <- system.file("examples", "meme.xml", package = "ggmotif")
#' motif_data <- get_motif_from_meme(data = filepath, format = "xml")
#' plot <- plot_motif_location(data = motif_data)
#'
#' # With phylogenetic tree
#' filepath <- system.file("examples", "meme.xml", package = "ggmotif")
#' tree_path <- system.file("examples", "ara.nwk", package = "ggmotif")
#' motif_data <- get_motif_from_meme(data = filepath, format = "xml")
#' plot <- plot_motif_location(data = motif_data, tree_path = tree_path)
#'
#' # With tree annotation
#' tree_anno <- data.frame(
#'   label = c("seq1", "seq2", "seq3"),
#'   Group = c("A", "A", "B")
#' )
#' plot <- plot_motif_location(
#'   data = motif_data,
#'   tree_path = tree_path,
#'   tree_annotation = tree_anno
#' )
#' }
#'
#' @export
plot_motif_location <- function(data,
                                 tree_path = NULL,
                                 tree_annotation = NULL) {
  # Input validation
  validate_motif_data(data)
  
  if (!is.null(tree_path) && !file.exists(tree_path)) {
    stop("Tree file not found: ", tree_path, call. = FALSE)
  }
  
  if (!is.null(tree_annotation)) {
    validate_tree_annotation(tree_annotation)
  }
  
  # Create plot based on whether tree is provided
  if (is.null(tree_path)) {
    create_motif_plot_standalone(data)
  } else {
    create_motif_plot_with_tree(data, tree_path, tree_annotation)
  }
}


#' Validate Motif Data
#'
#' @param data Data frame to validate
#' @keywords internal
#' @noRd
validate_motif_data <- function(data) {
  required_cols <- c(
    "input_seq_id", "length", "motif_id",
    "start_position", "end_position"
  )
  
  missing_cols <- setdiff(required_cols, names(data))
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in data: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}


#' Validate Tree Annotation
#'
#' @param tree_annotation Data frame to validate
#' @keywords internal
#' @noRd
validate_tree_annotation <- function(tree_annotation) {
  if (!is.data.frame(tree_annotation)) {
    stop("tree_annotation must be a data frame", call. = FALSE)
  }
  
  if (!"label" %in% names(tree_annotation)) {
    stop("tree_annotation must contain a 'label' column", call. = FALSE)
  }
  
  if (!"Group" %in% names(tree_annotation)) {
    stop("tree_annotation must contain a 'Group' column", call. = FALSE)
  }
  
  invisible(TRUE)
}


#' Prepare Motif Data for Plotting
#'
#' @param data Raw motif data
#' @param gene_order Optional data frame with Genes and y columns for ordering
#' @return Prepared data frame ready for plotting
#' @keywords internal
#' @noRd
prepare_motif_data <- function(data, gene_order = NULL) {
  # Rename columns to standardized names
  motif_data <- data %>%
    dplyr::select(
      .data$input_seq_id,
      .data$length,
      .data$motif_id,
      .data$start_position,
      .data$end_position
    ) %>%
    dplyr::rename(
      genes = .data$input_seq_id,
      motif = .data$motif_id,
      start = .data$start_position,
      end = .data$end_position
    )
  
  # Add y-coordinates for plotting
  if (is.null(gene_order)) {
    unique_genes <- unique(motif_data$genes)
    gene_order <- data.frame(
      genes = unique_genes,
      y = seq_along(unique_genes),
      stringsAsFactors = FALSE
    )
  }
  
  # Merge and calculate plot coordinates
  motif_data %>%
    dplyr::left_join(gene_order, by = "genes") %>%
    dplyr::mutate(
      length = as.numeric(.data$length),
      x_max_length = max(.data$length),
      x_min = .data$start,
      x_max = .data$end,
      y_min = .data$y - 0.4,
      y_max = .data$y + 0.4,
      group_id = paste0(.data$genes, .data$motif)
    ) %>%
    dplyr::arrange(.data$y, .data$length) %>%
    dplyr::mutate(genes = factor(.data$genes, levels = unique(.data$genes)))
}


#' Create Standalone Motif Plot
#'
#' @param data Motif data frame
#' @return A ggplot2 object
#' @keywords internal
#' @noRd
create_motif_plot_standalone <- function(data) {
  prepared_data <- prepare_motif_data(data)
  
  create_motif_ggplot(prepared_data, show_y_axis = TRUE)
}


#' Create Base Motif ggplot
#'
#' @param data Prepared motif data
#' @param show_y_axis Logical, whether to show y-axis elements
#' @return A ggplot2 object
#' @keywords internal
#' @noRd
create_motif_ggplot <- function(data, show_y_axis = TRUE) {
  max_length <- max(data$length)
  breaks <- seq(0, max_length, by = 100)
  
  p <- ggplot2::ggplot(data) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = 0,
        xend = .data$length,
        y = .data$genes,
        yend = .data$genes
      ),
      linewidth = 0.8
    ) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = .data$x_min,
        xmax = .data$x_max,
        ymin = .data$y_min,
        ymax = .data$y_max,
        fill = .data$motif
      )
    ) +
    ggplot2::scale_x_continuous(
      expand = c(0, 0),
      breaks = breaks
    ) +
    ggplot2::labs(x = "", y = "Name") +
    ggplot2::theme_classic()
  
  # Customize y-axis based on whether tree is present
  if (!show_y_axis) {
    p <- p +
      ggplot2::theme(
        axis.title.y = ggplot2::element_blank(),
        plot.margin = grid::unit(c(0, 0, 0, -3), "cm")
      )
  }
  
  return(p)
}


#' Create Phylogenetic Tree Plot
#'
#' @param tree_path Path to tree file
#' @param tree_annotation Optional annotation data frame
#' @return A ggtree object
#' @keywords internal
#' @noRd
create_tree_plot <- function(tree_path, tree_annotation = NULL) {
  tree <- ape::read.tree(tree_path)
  
  if (is.null(tree_annotation)) {
    p_tree <- ggtree::ggtree(tree, branch.length = "none") +
      ggplot2::theme(
        plot.margin = grid::unit(c(0, -3, 0, 0), "cm")
      )
  } else {
    p_tree <- tree %>%
      treeio::as_tibble() %>%
      dplyr::left_join(tree_annotation, by = "label") %>%
      treeio::as.treedata() %>%
      ggtree::ggtree(branch.length = "none") +
      ggtree::geom_tippoint(
        ggplot2::aes(color = .data$Group),
        size = 2
      ) +
      ggplot2::theme(
        plot.margin = grid::unit(c(0, -3, 0, 0), "cm")
      )
  }
  
  return(p_tree)
}


#' Extract Tree Tip Coordinates
#'
#' @param tree_plot A ggtree object
#' @return A data frame with gene names and y-coordinates
#' @keywords internal
#' @noRd
extract_tree_coordinates <- function(tree_plot) {
  tree_plot[["data"]] %>%
    dplyr::filter(.data$isTip == TRUE) %>%
    dplyr::select(.data$label, .data$y) %>%
    dplyr::rename(genes = .data$label) %>%
    dplyr::arrange(.data$y)
}


#' Create Motif Plot with Phylogenetic Tree
#'
#' @param data Motif data frame
#' @param tree_path Path to tree file
#' @param tree_annotation Optional annotation data frame
#' @return A combined plot with tree and motifs
#' @keywords internal
#' @noRd
create_motif_plot_with_tree <- function(data,
                                         tree_path,
                                         tree_annotation = NULL) {
  # Create tree plot
  tree_plot <- create_tree_plot(tree_path, tree_annotation)
  
  # Extract tree coordinates for gene ordering
  tree_coords <- extract_tree_coordinates(tree_plot)
  
  # Prepare motif data with tree-based ordering
  prepared_data <- data %>%
    prepare_motif_data(gene_order = tree_coords) %>%
    dplyr::mutate(
      genes = factor(.data$genes, levels = tree_coords$genes)
    )
  
  # Create motif plot
  motif_plot <- create_motif_ggplot(prepared_data, show_y_axis = FALSE)
  
  # Combine plots
  combined_plot <- tree_plot + motif_plot +
    patchwork::plot_layout(guides = "collect")
  
  return(combined_plot)
}


# Global variables (for R CMD check)
utils::globalVariables(c(
  ".", ".data"
))