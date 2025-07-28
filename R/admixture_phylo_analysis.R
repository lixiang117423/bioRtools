#' Analyze Admixture results with phylogenetic ordering and visualization
#'
#' @description
#' Process Admixture .Q files and create population structure plots ordered 
#' according to phylogenetic tree topology. This function integrates phylogenetic 
#' relationships with population genetic structure analysis for comprehensive 
#' evolutionary insights.
#'
#' @param admixture_path Character string specifying the directory path containing 
#'   Admixture .Q files. All .Q files in this directory will be processed.
#' @param tree_file Character string specifying the path to the phylogenetic tree 
#'   file (Newick format). Sample names in the tree must match those in the 
#'   Admixture files.
#' @param output_dir Character string specifying the output directory for plots. 
#'   Default is the same as admixture_path.
#' @param population_info Optional data frame containing population information 
#'   with columns 'sample' and 'population'. If provided, samples will be 
#'   grouped by population before phylogenetic ordering. Default is NULL.
#' @param plot_width Numeric value specifying the width of the admixture plot 
#'   in inches. Default is 15.
#' @param plot_height Numeric value specifying the height of the admixture plot 
#'   in inches. Default is 1.
#' @param tree_width Numeric value specifying the width of the phylogenetic tree 
#'   plot in inches. Default is 3.
#' @param tree_height Numeric value specifying the height of the phylogenetic 
#'   tree plot in inches. Default is 15.
#' @param dpi Numeric value specifying the resolution for output plots. 
#'   Default is 500.
#' @param k_range Integer vector specifying which K values to include in the 
#'   analysis. Can be specified as a range (e.g., 2:5) or specific values 
#'   (e.g., c(2,3,5)). If NULL, all available K values will be used. Default is NULL.
#' @param color_palette Character vector of colors for different ancestry 
#'   components. If NULL, uses a default color scheme. Colors will be 
#'   dynamically adjusted based on the maximum K value in the analysis.
#' @param show_tree_branch_length Logical indicating whether to show branch 
#'   lengths in the phylogenetic tree. Default is FALSE.
#' @param output_prefix Character string for output file naming prefix. 
#'   Default is "admixture_phylo_analysis".
#' @param align_clusters Logical indicating whether to align clusters across 
#'   different K values. Default is TRUE.
#' @param verbose Logical indicating whether to print progress information. 
#'   Default is TRUE.
#'
#' @return A list containing six components:
#' \describe{
#'   \item{qlist.original}{Original Q matrices read from files.}
#'   \item{qlist.ordered}{Q matrices reordered according to phylogenetic topology.}
#'   \item{tree.plot}{ggplot object of the phylogenetic tree.}
#'   \item{sample.order}{Data frame showing the final sample ordering with 
#'     phylogenetic and population information.}
#'   \item{summary.statistics}{Summary statistics from pophelper analysis.}
#'   \item{output.files}{Character vector of created output file paths.}
#' }
#'
#' @details
#' The function workflow includes:
#' \enumerate{
#'   \item Reading all .Q files from the specified directory
#'   \item Processing and summarizing admixture results
#'   \item Aligning clusters across different K values (optional)
#'   \item Reading phylogenetic tree and extracting tip order
#'   \item Reordering samples based on phylogenetic relationships
#'   \item Incorporating population structure if provided
#'   \item Creating publication-ready plots with matching dimensions
#' }
#'
#' The phylogenetic ordering ensures that:
#' \itemize{
#'   \item Closely related samples appear adjacent in plots
#'   \item Population structure patterns can be interpreted in evolutionary context
#'   \item Tree and admixture plots have matching sample orders for easy comparison
#' }
#'
#' @note
#' \itemize{
#'   \item Sample names must be consistent across .Q files and phylogenetic tree
#'   \item Requires pophelper, ape, ggtree, and dplyr packages
#'   \item Output plots are designed for easy horizontal concatenation
#'   \item Tree plot height should match admixture plot height for alignment
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic analysis with specific K range
#' result <- admixture_phylo_analysis(
#'   admixture_path = "path/to/admixture/results/",
#'   tree_file = "path/to/phylogenetic_tree.nwk",
#'   output_dir = "path/to/output/",
#'   k_range = 2:5  # Only analyze K=2,3,4,5
#' )
#'
#' # View the plots
#' result$tree.plot
#' 
#' # Access reordered data
#' head(result$sample.order)
#' 
#' # With population information and specific K values
#' pop_info <- data.frame(
#'   sample = c("sample1", "sample2", "sample3"),
#'   population = c("Pop1", "Pop1", "Pop2")
#' )
#' 
#' result_with_pop <- admixture_phylo_analysis(
#'   admixture_path = "path/to/admixture/results/",
#'   tree_file = "path/to/phylogenetic_tree.nwk",
#'   population_info = pop_info,
#'   k_range = c(2, 3, 5),  # Only K=2,3,5
#'   output_dir = "path/to/output/",
#'   plot_width = 20,
#'   tree_width = 4
#' )
#' 
#' # Custom colors with automatic adjustment
#' custom_colors <- c("#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4", 
#'                    "#FFEAA7", "#DDA0DD", "#98D8C8", "#F7DC6F")
#' 
#' result_custom <- admixture_phylo_analysis(
#'   admixture_path = "path/to/admixture/results/",
#'   tree_file = "path/to/phylogenetic_tree.nwk",
#'   k_range = 2:6,
#'   color_palette = custom_colors,  # Will extend if needed
#'   show_tree_branch_length = TRUE,
#'   dpi = 300
#' )
#' }
#'
admixture_phylo_analysis <- function(admixture_path,
                                     tree_file,
                                     output_dir = NULL,
                                     population_info = NULL,
                                     k_range = NULL,
                                     plot_width = 15,
                                     plot_height = 1,
                                     tree_width = 3,
                                     tree_height = 15,
                                     dpi = 500,
                                     color_palette = NULL,
                                     show_tree_branch_length = FALSE,
                                     output_prefix = "admixture_phylo_analysis",
                                     align_clusters = TRUE,
                                     verbose = TRUE) {
  
  # Check required packages
  required_packages <- c("pophelper", "ape", "ggtree", "dplyr", "ggplot2")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  
  if (length(missing_packages) > 0) {
    stop("Required packages missing: ", paste(missing_packages, collapse = ", "))
  }
  
  # Input validation
  if (!dir.exists(admixture_path)) {
    stop("Admixture directory not found: ", admixture_path)
  }
  
  if (!file.exists(tree_file)) {
    stop("Tree file not found: ", tree_file)
  }
  
  if (is.null(output_dir)) {
    output_dir <- admixture_path
  }
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    if (verbose) message("Created output directory: ", output_dir)
  }
  
  if (verbose) {
    message("Starting admixture phylogenetic analysis...")
    message("Admixture path: ", admixture_path)
    message("Tree file: ", tree_file)
    message("Output directory: ", output_dir)
  }
  
  # Step 1: Read admixture files
  if (verbose) message("Reading admixture .Q files...")
  
  file_list <- list.files(
    path = admixture_path,
    pattern = "*.Q$",
    full.names = TRUE
  )
  
  if (length(file_list) == 0) {
    stop("No .Q files found in: ", admixture_path)
  }
  
  if (verbose) message("Found ", length(file_list), " .Q files")
  
  tryCatch({
    qlist_original <- pophelper::readQ(files = file_list)
  }, error = function(e) {
    stop("Failed to read .Q files: ", e$message)
  })
  
  # Step 2: Process admixture data
  if (verbose) message("Processing admixture data...")
  
  tab_result <- pophelper::tabulateQ(qlist_original)
  sum_result <- pophelper::summariseQ(tab_result)
  
  # Step 3: Filter K values if specified and align clusters
  if (!is.null(k_range)) {
    if (verbose) message("Filtering K values: ", paste(k_range, collapse = ", "))
    
    # Get K values from Q matrices
    available_k <- sapply(qlist_original, ncol)
    names(qlist_original) <- paste0("K", available_k)
    
    # Filter based on specified range
    selected_k_names <- paste0("K", k_range)
    available_k_names <- names(qlist_original)
    
    # Check which K values are available
    missing_k <- setdiff(selected_k_names, available_k_names)
    if (length(missing_k) > 0) {
      warning("Requested K values not found: ", paste(gsub("K", "", missing_k), collapse = ", "))
    }
    
    # Filter to available K values
    valid_k_names <- intersect(selected_k_names, available_k_names)
    if (length(valid_k_names) == 0) {
      stop("None of the requested K values are available in the data")
    }
    
    qlist_filtered <- qlist_original[valid_k_names]
    if (verbose) message("Selected K values: ", paste(gsub("K", "", valid_k_names), collapse = ", "))
  } else {
    qlist_filtered <- qlist_original
    if (verbose) message("Using all available K values")
  }
  
  # Align clusters if requested
  if (align_clusters) {
    if (verbose) message("Aligning clusters across K values...")
    qlist_aligned <- pophelper::alignK(qlist_filtered)
  } else {
    qlist_aligned <- qlist_filtered
  }
  
  # Step 4: Read phylogenetic tree
  if (verbose) message("Reading phylogenetic tree...")
  
  tryCatch({
    tree <- ape::read.tree(tree_file)
  }, error = function(e) {
    stop("Failed to read tree file: ", e$message)
  })
  
  tree_order <- tree$tip.label
  if (verbose) message("Tree contains ", length(tree_order), " tips")
  
  # Step 5: Create sample ordering data
  if (verbose) message("Creating phylogenetic ordering...")
  
  # Get sample names from Q matrices
  sample_names <- rownames(qlist_aligned[[1]])
  
  # Create base ordering data
  phylo_data <- data.frame(
    sample = sample_names,
    phylo_order = match(sample_names, tree_order),
    stringsAsFactors = FALSE
  )
  
  # Add population information if provided
  if (!is.null(population_info)) {
    if (verbose) message("Incorporating population information...")
    
    if (!all(c("sample", "population") %in% colnames(population_info))) {
      stop("population_info must contain 'sample' and 'population' columns")
    }
    
    phylo_data <- phylo_data %>%
      dplyr::left_join(population_info, by = "sample")
    
    # Order by population first, then by phylogenetic order within population
    ordered_data <- phylo_data %>%
      dplyr::arrange(population, phylo_order) %>%
      dplyr::filter(!is.na(phylo_order))
  } else {
    # Order only by phylogenetic relationships
    ordered_data <- phylo_data %>%
      dplyr::arrange(phylo_order) %>%
      dplyr::filter(!is.na(phylo_order))
  }
  
  if (nrow(ordered_data) == 0) {
    stop("No samples match between admixture files and phylogenetic tree")
  }
  
  if (verbose) {
    message("Successfully matched ", nrow(ordered_data), " samples")
    missing_samples <- setdiff(sample_names, ordered_data$sample)
    if (length(missing_samples) > 0) {
      message("Warning: ", length(missing_samples), " samples not found in tree: ", 
              paste(head(missing_samples, 5), collapse = ", "))
    }
  }
  
  # Step 6: Reorder Q matrices
  if (verbose) message("Reordering admixture data...")
  
  final_order_indices <- match(ordered_data$sample, sample_names)
  qlist_ordered <- lapply(qlist_aligned, function(x) {
    x[final_order_indices, , drop = FALSE]
  })
  
  # Step 7: Set up colors dynamically based on maximum K
  max_k <- max(sapply(qlist_ordered, ncol))
  if (verbose) message("Maximum K value: ", max_k)
  
  if (is.null(color_palette)) {
    # Create a comprehensive color palette
    base_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", 
                     "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
                     "#bcbd22", "#17becf", "#aec7e8", "#ffbb78",
                     "#98df8a", "#ff9896", "#c5b0d5", "#c49c94",
                     "#f7b6d3", "#c7c7c7", "#dbdb8d", "#9edae5")
    
    # If we need more colors than available, generate additional colors
    if (max_k > length(base_colors)) {
      additional_colors <- rainbow(max_k - length(base_colors))
      color_palette <- c(base_colors, additional_colors)
    } else {
      color_palette <- base_colors[1:max_k]
    }
  } else {
    # Ensure we have enough colors for the maximum K
    if (length(color_palette) < max_k) {
      warning("Color palette has fewer colors (", length(color_palette), 
              ") than maximum K (", max_k, "). Extending with default colors.")
      additional_colors <- rainbow(max_k - length(color_palette))
      color_palette <- c(color_palette, additional_colors)
    } else {
      # Use only the needed number of colors
      color_palette <- color_palette[1:max_k]
    }
  }
  
  if (verbose) message("Using ", length(color_palette), " colors for visualization")
  
  # Step 8: Create admixture plot
  if (verbose) message("Creating admixture structure plot...")
  
  admixture_output_file <- file.path(output_dir, paste0(output_prefix, "_structure"))
  
  tryCatch({
    pophelper::plotQ(
      qlist_ordered,
      imgoutput = "join",
      exportpath = output_dir,
      height = plot_height,
      width = plot_width,
      imgtype = "png",
      dpi = dpi,
      clustercol = color_palette,
      splab = paste0("K=", sapply(qlist_ordered, ncol)),
      showyaxis = FALSE,
      showlegend = FALSE,
      outputfilename = basename(admixture_output_file)
    )
  }, error = function(e) {
    warning("Failed to create admixture plot: ", e$message)
  })
  
  # Step 9: Create phylogenetic tree plot
  if (verbose) message("Creating phylogenetic tree plot...")
  
  # Reorder tree to match sample order
  tree_subset <- ape::keep.tip(tree, ordered_data$sample)
  
  tree_plot <- ggtree::ggtree(tree_subset, branch.length = if(show_tree_branch_length) NULL else "none") +
    ggplot2::theme(
      plot.margin = ggplot2::margin(5, 5, 5, 5, "pt")
    )
  
  tree_output_file <- file.path(output_dir, paste0(output_prefix, "_phylogeny.png"))
  
  tryCatch({
    ggplot2::ggsave(
      tree_plot,
      filename = tree_output_file,
      width = tree_width,
      height = tree_height,
      dpi = dpi,
      bg = "white"
    )
  }, error = function(e) {
    warning("Failed to save tree plot: ", e$message)
  })
  
  # Step 10: Create output file list
  output_files <- c(
    paste0(admixture_output_file, ".png"),
    tree_output_file
  )
  
  # Filter existing files
  output_files <- output_files[file.exists(output_files)]
  
  if (verbose) {
    message("Analysis completed successfully!")
    message("Output files created:")
    for (file in output_files) {
      message("  ", file)
    }
  }
  
  # Return comprehensive results
  return(list(
    qlist.original = qlist_original,
    qlist.filtered = if(!is.null(k_range)) qlist_filtered else NULL,
    qlist.ordered = qlist_ordered,
    tree.plot = tree_plot,
    sample.order = ordered_data,
    summary.statistics = sum_result,
    k.values.used = sapply(qlist_ordered, ncol),
    color.palette.used = color_palette,
    output.files = output_files
  ))
}