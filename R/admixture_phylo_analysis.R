#' Analyze Admixture results with phylogenetic ordering and ggplot2 visualization
#'
#' @description
#' Process Admixture .Q files and create population structure plots ordered 
#' according to phylogenetic tree topology using ggplot2. This function integrates 
#' phylogenetic relationships with population genetic structure analysis and 
#' returns customizable ggplot2 objects.
#'
#' @param admixture_path Character string specifying the directory path containing 
#'   Admixture .Q files and .fam file. All .Q files in this directory will be processed.
#' @param tree_file Character string specifying the path to the phylogenetic tree 
#'   file (Newick format). Sample names in the tree must match those in the 
#'   .fam file.
#' @param population_info Optional data frame containing population information 
#'   with columns 'sample' and 'population'. If provided, samples will be 
#'   grouped by population before phylogenetic ordering. Default is NULL.
#' @param k_range Integer vector specifying which K values to include in the 
#'   analysis. Can be specified as a range (e.g., 2:5) or specific values 
#'   (e.g., c(2,3,5)). If NULL, all available K values will be used. Default is NULL.
#' @param color_scale Character string specifying the color scale to use. 
#'   Options: "default" (uses optimized base colors), "aaas", "npg", "lancet", 
#'   "jco", "ucscgb", "uchicago", "simpsons", "rickandmorty", or "manual". 
#'   Default is "default".
#' @param manual_colors Character vector of colors for manual color scale. 
#'   Only used when color_scale = "manual". Default is NULL.
#' @param verbose Logical indicating whether to print progress information. 
#'   Default is TRUE.
#'
#' @return A ggplot2 object of the admixture structure plot, ready for display 
#' or further customization.
#'
#' @details
#' The function workflow includes:
#' \enumerate{
#'   \item Reading phylogenetic tree and extracting tip order with proper label cleaning
#'   \item Reading .fam file to get sample information
#'   \item Reading all .Q files and combining with sample names
#'   \item Converting data to long format suitable for ggplot2
#'   \item Creating publication-ready ggplot2 visualization
#'   \item Applying phylogenetic ordering to ensure related samples are adjacent
#' }
#'
#' The default color palette intelligently selects colors based on the maximum K value:
#' \itemize{
#'   \item Uses a curated set of 20 distinct colors
#'   \item Automatically extends palette if more colors are needed
#'   \item Optimized for scientific publication and color-blind accessibility
#' }
#'
#' @note
#' \itemize{
#'   \item Requires tidyverse, ape, treeio, and ggsci packages
#'   \item Returns ggplot2 object for easy customization
#'   \item Default styling removes axes, legends, and panel spacing for clean look
#'   \item Sample names are automatically cleaned (quotes removed from tree labels)
#'   \item No files are saved automatically - use ggsave() if needed
#' }
#'
#' @export
#' @examples
#' \dontrun{
#' # Basic analysis with default colors
#' plot <- admixture_phylo_analysis(
#'   admixture_path = "path/to/admixture/results/",
#'   tree_file = "path/to/phylogenetic_tree.nwk"
#' )
#' print(plot)
#' 
#' # Customize the plot further
#' plot + 
#'   theme(strip.text.y = element_text(size = 12)) +
#'   labs(title = "Population Structure Analysis")
#' 
#' # Use ggsci color scale
#' plot_npg <- admixture_phylo_analysis(
#'   admixture_path = "path/to/admixture/results/",
#'   tree_file = "path/to/phylogenetic_tree.nwk",
#'   color_scale = "npg",
#'   k_range = 2:6
#' )
#' 
#' # Custom colors
#' plot_custom <- admixture_phylo_analysis(
#'   admixture_path = "path/to/admixture/results/",
#'   tree_file = "path/to/phylogenetic_tree.nwk",
#'   color_scale = "manual",
#'   manual_colors = c("#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4")
#' )
#' 
#' # Save if needed
#' ggsave("admixture_plot.png", plot, width = 15, height = 8, dpi = 500)
#' }
#'
admixture_phylo_analysis <- function(admixture_path,
                                     tree_file,
                                     population_info = NULL,
                                     k_range = NULL,
                                     color_scale = "default",
                                     manual_colors = NULL,
                                     verbose = TRUE) {
  
  # Check required packages
  required_packages <- c("dplyr", "tidyr", "ggplot2", "readr", "stringr", 
                        "purrr", "magrittr", "ape", "treeio", "ggsci")
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
  
  if (verbose) {
    message("Starting admixture phylogenetic analysis...")
    message("Admixture path: ", admixture_path)
    message("Tree file: ", tree_file)
    message("Color scale: ", color_scale)
  }
  
  # Step 1: Read phylogenetic tree and get sample order
  if (verbose) message("Reading phylogenetic tree and extracting sample order...")
  
  tryCatch({
    df.order.tree <- ape::read.tree(tree_file) %>% 
      treeio::as_tibble() %>% 
      dplyr::filter(!is.na(label)) %>% 
      dplyr::select(node, label) %>% 
      dplyr::mutate(sample = stringr::str_remove_all(label, "\\'")) %>% 
      dplyr::select(sample, node) %>% 
      dplyr::rename(order = node)
    
    if (verbose) message("Found ", nrow(df.order.tree), " samples in tree")
  }, error = function(e) {
    stop("Failed to read tree file: ", e$message)
  })
  
  # Step 2: Read .fam file to get sample information
  if (verbose) message("Reading .fam file...")
  
  tryCatch({
    df.sample <- dir(admixture_path, pattern = "fam") %>% 
      as.data.frame() %>% 
      magrittr::set_names("file") %>% 
      dplyr::mutate(path = file.path(admixture_path, file)) %>% 
      dplyr::mutate(fam = purrr::map(path, ~ readr::read_delim(.x, col_names = FALSE, show_col_types = FALSE) %>% 
                                           dplyr::select(1))) %>% 
      tidyr::unnest(fam) %>% 
      dplyr::select(3) %>% 
      magrittr::set_names("sample")
    
    if (verbose) message("Read ", nrow(df.sample), " samples from .fam file")
  }, error = function(e) {
    stop("Failed to read .fam file: ", e$message)
  })
  
  # Step 3: Read .Q files and create admixture data
  if (verbose) message("Reading .Q files...")
  
  tryCatch({
    q_files <- dir(admixture_path, pattern = "Q$")
    
    if (length(q_files) == 0) {
      stop("No .Q files found in: ", admixture_path)
    }
    
    df.admixture <- q_files %>% 
      as.data.frame() %>% 
      magrittr::set_names("file") %>% 
      dplyr::mutate(k_value = stringr::str_split(file, "\\.") %>% sapply("[", 2),
                    k_value = as.numeric(k_value),
                    k = paste0("K = ", k_value),
                    path = file.path(admixture_path, file)) %>% 
      dplyr::mutate(q_matrix = purrr::map(path, ~ readr::read_delim(.x, col_names = FALSE, show_col_types = FALSE) %>% 
                                                magrittr::set_names(paste0("Cluster", 1:ncol(.))) %>% 
                                                dplyr::mutate(sample = df.sample$sample[1:nrow(.)]) %>% 
                                                dplyr::select(sample, 1:(ncol(.)-1)) %>% 
                                                tidyr::pivot_longer(cols = 2:ncol(.)))) %>% 
      dplyr::select(k, k_value, q_matrix) %>% 
      tidyr::unnest(q_matrix)
    
    if (verbose) message("Found K values: ", paste(sort(unique(df.admixture$k_value)), collapse = ", "))
  }, error = function(e) {
    stop("Failed to read .Q files: ", e$message)
  })
  
  # Step 4: Filter K values if specified
  if (!is.null(k_range)) {
    if (verbose) message("Filtering K values: ", paste(k_range, collapse = ", "))
    
    df.admixture <- df.admixture %>% 
      dplyr::filter(k_value %in% k_range)
    
    if (nrow(df.admixture) == 0) {
      stop("No data found for specified K values: ", paste(k_range, collapse = ", "))
    }
    
    if (verbose) message("Using K values: ", paste(sort(unique(df.admixture$k_value)), collapse = ", "))
  }
  
  # Step 5: Add population information if provided
  if (!is.null(population_info)) {
    if (verbose) message("Incorporating population information...")
    
    if (!all(c("sample", "population") %in% colnames(population_info))) {
      stop("population_info must contain 'sample' and 'population' columns")
    }
    
    # Add population info to tree order
    df.order.tree <- df.order.tree %>% 
      dplyr::left_join(population_info, by = "sample") %>% 
      dplyr::arrange(population, order)
  }
  
  # Step 6: Check sample matching
  matched_samples <- intersect(df.admixture$sample, df.order.tree$sample)
  if (length(matched_samples) == 0) {
    cat("\n=== DEBUGGING INFORMATION ===\n")
    cat("Admixture samples (first 10):\n")
    print(head(unique(df.admixture$sample), 10))
    cat("\nTree samples (first 10):\n")
    print(head(df.order.tree$sample, 10))
    cat("=============================\n")
    stop("No samples match between .Q files and phylogenetic tree")
  }
  
  if (verbose) message("Successfully matched ", length(matched_samples), " samples")
  
  # Step 7: Determine colors needed and prepare palette
  max_clusters <- max(table(df.admixture$k, df.admixture$name))
  max_k <- max(df.admixture$k_value)
  
  # Calculate total unique clusters across all K values
  total_clusters <- length(unique(df.admixture$name))
  
  if (verbose) message("Maximum K value: ", max_k, ", Total clusters needed: ", total_clusters)
  
  # Define base color palette (your provided colors)
  base_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", 
                   "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
                   "#bcbd22", "#17becf", "#aec7e8", "#ffbb78",
                   "#98df8a", "#ff9896", "#c5b0d5", "#c49c94",
                   "#f7b6d3", "#c7c7c7", "#dbdb8d", "#9edae5")
  
  # Step 8: Create the plot
  if (verbose) message("Creating ggplot2 visualization...")
  
  p <- df.admixture %>% 
    ggplot2::ggplot(ggplot2::aes(sample, value, fill = name)) +
    ggplot2::geom_bar(stat = "identity", width = 1, position = ggplot2::position_fill()) +
    ggplot2::scale_x_discrete(limits = df.order.tree$sample) +
    ggplot2::facet_wrap(. ~ k, ncol = 1, strip.position = "left") +
    ggplot2::theme(
      # Clean appearance - no axes, ticks, or legends
      legend.position = "none",
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      
      # Facet styling
      strip.text.y = ggplot2::element_text(angle = 90),
      strip.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      
      # Minimal spacing between facets
      panel.spacing.y = ggplot2::unit(0, "lines")
    )
  
  # Step 9: Apply color scale
  if (color_scale == "default") {
    # Use base colors, extending if necessary
    if (total_clusters <= length(base_colors)) {
      selected_colors <- base_colors[1:total_clusters]
    } else {
      # Extend with additional colors if needed
      additional_colors <- rainbow(total_clusters - length(base_colors))
      selected_colors <- c(base_colors, additional_colors)
      if (verbose) message("Extended color palette with ", length(additional_colors), " additional colors")
    }
    p <- p + ggplot2::scale_fill_manual(values = selected_colors)
    
  } else if (color_scale == "manual") {
    if (is.null(manual_colors)) {
      warning("manual_colors not provided, using default colors")
      if (total_clusters <= length(base_colors)) {
        selected_colors <- base_colors[1:total_clusters]
      } else {
        additional_colors <- rainbow(total_clusters - length(base_colors))
        selected_colors <- c(base_colors, additional_colors)
      }
      p <- p + ggplot2::scale_fill_manual(values = selected_colors)
    } else {
      # Check if enough manual colors provided
      if (length(manual_colors) < total_clusters) {
        warning("Not enough manual colors provided (", length(manual_colors), 
                ") for clusters needed (", total_clusters, "). Extending with default colors.")
        extended_colors <- c(manual_colors, base_colors[1:(total_clusters - length(manual_colors))])
        p <- p + ggplot2::scale_fill_manual(values = extended_colors)
      } else {
        p <- p + ggplot2::scale_fill_manual(values = manual_colors[1:total_clusters])
      }
    }
    
  } else {
    # Use ggsci color scales
    p <- switch(color_scale,
      "aaas" = p + ggsci::scale_fill_aaas(),
      "npg" = p + ggsci::scale_fill_npg(),
      "lancet" = p + ggsci::scale_fill_lancet(),
      "jco" = p + ggsci::scale_fill_jco(),
      "ucscgb" = p + ggsci::scale_fill_ucscgb(),
      "uchicago" = p + ggsci::scale_fill_uchicago(),
      "simpsons" = p + ggsci::scale_fill_simpsons(),
      "rickandmorty" = p + ggsci::scale_fill_rickandmorty(),
      {
        warning("Unknown color scale: ", color_scale, ". Using default colors.")
        if (total_clusters <= length(base_colors)) {
          selected_colors <- base_colors[1:total_clusters]
        } else {
          additional_colors <- rainbow(total_clusters - length(base_colors))
          selected_colors <- c(base_colors, additional_colors)
        }
        p + ggplot2::scale_fill_manual(values = selected_colors)
      }
    )
  }
  
  if (verbose) message("Analysis completed successfully! Returning ggplot2 object.")
  
  # Return the ggplot2 object directly
  return(p)
}