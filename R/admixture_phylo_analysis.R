#' Analyze Admixture results with phylogenetic tree and visualization
#'
#' @description
#' Process Admixture .Q files and create population structure plots with integrated 
#' phylogenetic tree visualization using ggplot2 and ggtree. This function creates 
#' a rotated layout with samples on Y-axis and values on X-axis, with K values 
#' arranged horizontally, and optionally combines with a phylogenetic tree.
#'
#' @param admixture_path Character string specifying the directory path containing 
#'   Admixture .Q files and .fam file. All .Q files in this directory will be processed.
#' @param tree_file Optional character string specifying the path to the phylogenetic tree 
#'   file (Newick format). If provided, a phylogenetic tree will be drawn and combined 
#'   with the admixture plot. If NULL, only admixture plot will be returned. Default is NULL.
#' @param sample_group Optional data frame containing sample grouping information 
#'   with columns 'label' (sample names), 'group', and optionally 'color'. 
#'   Can also be a file path to read the grouping information.
#'   If provided when tree_file is given, tree tips will be colored by groups. Default is NULL.
#' @param k_range Integer vector specifying which K values to include in the 
#'   analysis. Can be specified as a range (e.g., 2:5) or specific values 
#'   (e.g., c(2,3,5)). If NULL, all available K values will be used. Default is NULL.
#' @param k_ncol Integer specifying number of columns for K value facets. 
#'   Default is 9.
#' @param color_scale Character string specifying the color scale to use for admixture. 
#'   Options: "default" (uses optimized base colors), "aaas", "npg", "lancet", 
#'   "jco", "ucscgb", "uchicago", "simpsons", "rickandmorty", or "manual". 
#'   Default is "aaas".
#' @param tree_color_scale Character string specifying the color scale for tree groups.
#'   Options: "d3", "aaas", "npg", "lancet", etc. Default is "d3".
#' @param manual_colors Character vector of colors for manual admixture color scale. 
#'   Only used when color_scale = "manual". Default is NULL.
#' @param tree_width Numeric value specifying the relative width of the tree plot. 
#'   Default is 0.2 (20% of total width).
#' @param tree_linewidth Numeric value specifying the line width for tree branches. 
#'   Default is 0.2.
#' @param tip_point_size Numeric value specifying the size of tip points in tree. 
#'   Default is 0.1.
#' @param show_branch_length Logical indicating whether to show branch lengths. 
#'   Default is FALSE.
#' @param verbose Logical indicating whether to print progress information. 
#'   Default is TRUE.
#'
#' @return A list containing:
#' \describe{
#'   \item{combined_plot}{Combined plot object (aplot) when tree_file provided, otherwise NULL.}
#'   \item{admixture_plot}{ggplot2 object of the admixture structure plot.}
#'   \item{tree_plot}{ggplot2/ggtree object when tree_file provided, otherwise NULL.}
#' }
#'
#' @details
#' The function has two distinct modes based on whether tree_file is provided:
#' 
#' **Mode 1: Without tree_file (Traditional Layout)**
#' \enumerate{
#'   \item Samples arranged horizontally (X-axis)
#'   \item K values arranged vertically as facets (ncol = 1)
#'   \item Strip labels on the left with 90-degree rotation
#'   \item Returns only admixture plot
#' }
#' 
#' **Mode 2: With tree_file (Rotated Layout)**
#' \enumerate{
#'   \item Reading phylogenetic tree and extracting tip order with proper label cleaning
#'   \item Samples arranged vertically (Y-axis)
#'   \item K values arranged horizontally as facets (ncol = auto-calculated from k_range)
#'   \item Strip labels on top with 90-degree rotation
#'   \item Optional sample grouping with colored tree tips
#'   \item Returns combined plot, admixture plot, and tree plot
#' }
#' 
#' The sample_group parameter can be:
#' \itemize{
#'   \item A data.frame with columns: 'label', 'group', and optionally 'color'
#'   \item A file path to a tab-delimited file with the same columns
#'   \item NULL for no grouping (default)
#' }
#'
#' @note
#' \itemize{
#'   \item Requires tidyverse and ggplot2 packages
#'   \item When using phylogenetic features: requires ape, treeio, ggtree, aplot packages
#'   \item When using ggsci color scales: requires ggsci package  
#'   \item Returns aplot object when tree provided, ggplot2 object otherwise
#'   \item No files are saved automatically - use ggsave() if needed
#' }
#'
#' @export
#' @examples
#' \dontrun{
#' # Mode 1: Analysis without tree (traditional layout)
#' result <- admixture_phylo_analysis(
#'   admixture_path = "path/to/admixture/results/"
#' )
#' print(result$admixture_plot)  # ggplot2 object
#' # result$combined_plot and result$tree_plot will be NULL
#' 
#' # Mode 2: Analysis with tree (rotated layout)
#' result <- admixture_phylo_analysis(
#'   admixture_path = "path/to/admixture/results/",
#'   tree_file = "path/to/phylogenetic_tree.nwk"
#' )
#' print(result$combined_plot)   # aplot combined object
#' print(result$admixture_plot)  # ggplot2 admixture plot
#' print(result$tree_plot)       # ggtree plot
#' 
#' # Mode 2 with sample grouping (data frame)
#' sample_group <- data.frame(
#'   label = c("sample1", "sample2", "sample3"),
#'   group = c("A", "A", "B"),
#'   color = c("A", "A", "B")
#' )
#' 
#' result <- admixture_phylo_analysis(
#'   admixture_path = "path/to/admixture/results/",
#'   tree_file = "path/to/phylogenetic_tree.nwk",
#'   sample_group = sample_group,
#'   k_range = 2:6
#' )
#' 
#' # Mode 2 with sample grouping (file path)
#' result <- admixture_phylo_analysis(
#'   admixture_path = "path/to/admixture/results/",
#'   tree_file = "path/to/phylogenetic_tree.nwk",
#'   sample_group = "path/to/sample_groups.txt"  # Tab-delimited file
#' )
#' 
#' # Custom styling
#' result <- admixture_phylo_analysis(
#'   admixture_path = "path/to/admixture/results/",
#'   tree_file = "path/to/phylogenetic_tree.nwk",
#'   color_scale = "npg",
#'   tree_color_scale = "aaas",
#'   tree_width = 0.3,
#'   tip_point_size = 0.2
#' )
#' }
#'
admixture_phylo_analysis <- function(admixture_path,
                                     tree_file = NULL,
                                     sample_group = NULL,
                                     k_range = NULL,
                                     color_scale = "aaas",
                                     tree_color_scale = "d3",
                                     manual_colors = NULL,
                                     tree_width = 0.2,
                                     tree_linewidth = 0.2,
                                     tip_point_size = 0.1,
                                     show_branch_length = FALSE,
                                     verbose = TRUE) {
  
  # Check required packages
  required_packages <- c("dplyr", "tidyr", "ggplot2", "readr", "stringr", 
                        "purrr", "magrittr")
  
  # Add tree-related packages if tree file is provided
  if (!is.null(tree_file)) {
    required_packages <- c(required_packages, "ape", "treeio", "ggtree", "aplot")
  }
  
  # Add ggsci if using ggsci color scales
  ggsci_scales <- c("aaas", "npg", "lancet", "jco", "ucscgb", "uchicago", 
                    "simpsons", "rickandmorty", "d3")
  if (color_scale %in% ggsci_scales || tree_color_scale %in% ggsci_scales) {
    required_packages <- c(required_packages, "ggsci")
  }
  
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  
  if (length(missing_packages) > 0) {
    stop("Required packages missing: ", paste(missing_packages, collapse = ", "))
  }
  
  # Input validation
  if (!dir.exists(admixture_path)) {
    stop("Admixture directory not found: ", admixture_path)
  }
  
  if (!is.null(tree_file) && !file.exists(tree_file)) {
    stop("Tree file not found: ", tree_file)
  }
  
  if (verbose) {
    message("Starting admixture phylogenetic analysis...")
    message("Admixture path: ", admixture_path)
    if (!is.null(tree_file)) {
      message("Tree file: ", tree_file)
    } else {
      message("No tree file provided - creating admixture plot only")
    }
    message("Color scale: ", color_scale)
  }
  
  # Step 1: Read phylogenetic tree and get sample order (if provided)
  if (!is.null(tree_file)) {
    if (verbose) message("Reading phylogenetic tree and extracting sample order...")
    
    tryCatch({
      df.order.tree <- ape::read.tree(tree_file) %>% 
        treeio::as_tibble() %>% 
        dplyr::filter(!is.na(label)) %>% 
        dplyr::select(node, label) %>% 
        dplyr::mutate(sample = stringr::str_remove_all(label, "\\'")) %>% 
        dplyr::select(sample, node) %>% 
        as.data.frame() %>% 
        dplyr::rename(order = node)
      
      if (verbose) message("Found ", nrow(df.order.tree), " samples in tree")
    }, error = function(e) {
      stop("Failed to read tree file: ", e$message)
    })
  } else {
    if (verbose) message("No tree file provided - will use default sample order from .fam file")
    df.order.tree <- NULL  # Will be created after reading .fam file
  }
  
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
    
    # If no tree file provided, create default ordering based on .fam file order
    if (is.null(df.order.tree)) {
      df.order.tree <- data.frame(
        sample = df.sample$sample,
        order = seq_len(nrow(df.sample)),
        stringsAsFactors = FALSE
      )
      if (verbose) message("Created default sample ordering based on .fam file order")
    }
    
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
  
  # Fix K value ordering by converting to factor with proper levels
  df.admixture <- df.admixture %>% 
    dplyr::mutate(
      k = factor(k, levels = paste0("K = ", sort(unique(k_value))))
    )
  
  # Step 5: Check sample matching
  if (!is.null(tree_file)) {
    # Only check matching if tree file was provided
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
    if (verbose) message("Successfully matched ", length(matched_samples), " samples with tree")
  } else {
    # When no tree file, all samples should match since we created df.order.tree from .fam file
    matched_samples <- intersect(df.admixture$sample, df.order.tree$sample)
    if (verbose) message("Using ", length(matched_samples), " samples in default order")
  }
  
  # Step 6: Determine colors needed and prepare palette
  total_clusters <- length(unique(df.admixture$name))
  
  if (verbose) message("Total clusters needed: ", total_clusters)
  
  # Define base color palette
  base_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", 
                   "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
                   "#bcbd22", "#17becf", "#aec7e8", "#ffbb78",
                   "#98df8a", "#ff9896", "#c5b0d5", "#c49c94",
                   "#f7b6d3", "#c7c7c7", "#dbdb8d", "#9edae5")
  
  # Step 7: Create the admixture plot with different layouts based on tree_file
  if (verbose) message("Creating admixture plot...")
  
  if (is.null(tree_file)) {
    # Mode 1: No tree file - Traditional layout (sample on X-axis, vertical facets)
    if (verbose) message("Using traditional layout (no tree): samples on X-axis, K values vertically")
    
    p.admixture <- df.admixture %>% 
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
        
        # Remove all background elements
        panel.background = ggplot2::element_blank(),
        plot.background = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        
        # Facet styling (vertical layout)
        strip.text.y = ggplot2::element_text(angle = 90),
        strip.background = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        
        # Minimal spacing between facets (vertical)
        panel.spacing.y = ggplot2::unit(0, "lines")
      )
    
  } else {
    # Mode 2: With tree file - Rotated layout (sample on Y-axis, horizontal facets)
    if (verbose) message("Using rotated layout (with tree): samples on Y-axis, K values horizontally")
    
    # Calculate ncol based on k_range
    if (!is.null(k_range)) {
      k_ncol <- length(k_range)
    } else {
      k_ncol <- length(unique(df.admixture$k_value))
    }
    
    if (verbose) message("Setting facet columns to: ", k_ncol)
    
    p.admixture <- df.admixture %>% 
      ggplot2::ggplot(ggplot2::aes(value, sample, fill = name)) +
      ggplot2::geom_bar(stat = "identity", width = 1, position = ggplot2::position_fill()) +
      ggplot2::scale_y_discrete(limits = df.order.tree$sample) +
      ggplot2::facet_wrap(k ~ ., ncol = k_ncol, strip.position = "top") +
      ggplot2::theme(
        # Clean appearance - no axes, ticks, or legends
        legend.position = "none",
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        
        # Remove all background elements
        panel.background = ggplot2::element_blank(),
        plot.background = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        
        # Facet styling (horizontal layout)
        strip.text.x = ggplot2::element_text(angle = 90),
        strip.background = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        
        # Minimal spacing between facets (horizontal)
        panel.spacing.x = ggplot2::unit(0, "lines")
      )
  }
  
  # Step 8: Apply color scale to admixture plot
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
    p.admixture <- p.admixture + ggplot2::scale_fill_manual(values = selected_colors)
    
  } else if (color_scale == "manual") {
    if (is.null(manual_colors)) {
      warning("manual_colors not provided, using default colors")
      if (total_clusters <= length(base_colors)) {
        selected_colors <- base_colors[1:total_clusters]
      } else {
        additional_colors <- rainbow(total_clusters - length(base_colors))
        selected_colors <- c(base_colors, additional_colors)
      }
      p.admixture <- p.admixture + ggplot2::scale_fill_manual(values = selected_colors)
    } else {
      # Check if enough manual colors provided
      if (length(manual_colors) < total_clusters) {
        warning("Not enough manual colors provided (", length(manual_colors), 
                ") for clusters needed (", total_clusters, "). Extending with default colors.")
        extended_colors <- c(manual_colors, base_colors[1:(total_clusters - length(manual_colors))])
        p.admixture <- p.admixture + ggplot2::scale_fill_manual(values = extended_colors)
      } else {
        p.admixture <- p.admixture + ggplot2::scale_fill_manual(values = manual_colors[1:total_clusters])
      }
    }
    
  } else {
    # Use ggsci color scales
    p.admixture <- switch(color_scale,
      "aaas" = p.admixture + ggsci::scale_fill_aaas(),
      "npg" = p.admixture + ggsci::scale_fill_npg(),
      "lancet" = p.admixture + ggsci::scale_fill_lancet(),
      "jco" = p.admixture + ggsci::scale_fill_jco(),
      "ucscgb" = p.admixture + ggsci::scale_fill_ucscgb(),
      "uchicago" = p.admixture + ggsci::scale_fill_uchicago(),
      "simpsons" = p.admixture + ggsci::scale_fill_simpsons(),
      "rickandmorty" = p.admixture + ggsci::scale_fill_rickandmorty(),
      {
        warning("Unknown color scale: ", color_scale, ". Using default colors.")
        if (total_clusters <= length(base_colors)) {
          selected_colors <- base_colors[1:total_clusters]
        } else {
          additional_colors <- rainbow(total_clusters - length(base_colors))
          selected_colors <- c(base_colors, additional_colors)
        }
        p.admixture + ggplot2::scale_fill_manual(values = selected_colors)
      }
    )
  }
  
  # Step 9: Create phylogenetic tree plot and combine (if tree file provided)
  p.tree <- NULL
  combined_plot <- NULL
  
  if (!is.null(tree_file)) {
    if (verbose) message("Creating phylogenetic tree plot...")
    
    # Read sample_group information if provided
    sample_group_data <- NULL
    if (!is.null(sample_group)) {
      if (is.character(sample_group) && length(sample_group) == 1) {
        # It's a file path
        if (verbose) message("Reading sample group information from file: ", sample_group)
        tryCatch({
          sample_group_data <- readr::read_delim(sample_group, delim = "\t", show_col_types = FALSE)
        }, error = function(e) {
          warning("Failed to read sample group file: ", e$message)
          sample_group_data <- NULL
        })
      } else if (is.data.frame(sample_group)) {
        # It's a data frame
        sample_group_data <- sample_group
      } else {
        warning("sample_group must be a file path or data frame")
      }
    }
    
    # Create base tree plot
    p.tree <- ape::read.tree(tree_file) %>% 
      treeio::as_tibble() %>% 
      dplyr::mutate(label = stringr::str_remove_all(label, "\\'")) %>% 
      treeio::as.phylo() %>% 
      ggtree::ggtree(branch.length = if(show_branch_length) NULL else "none", 
                    ladderize = FALSE,
                     linewidth = tree_linewidth)
    
    # Add group information if provided
    if (!is.null(sample_group_data)) {
      if (verbose) message("Adding sample group information to tree...")
      
      # Validate sample_group_data
      required_cols <- c("label", "group")
      if (!all(required_cols %in% colnames(sample_group_data))) {
        warning("sample_group must contain 'label' and 'group' columns. Skipping group coloring.")
      } else {
        # Use 'color' column if available, otherwise use 'group'
        if (!"color" %in% colnames(sample_group_data)) {
          sample_group_data$color <- sample_group_data$group
        }
        
        # Add group info to tree
        p.tree <- p.tree %<+% sample_group_data +
          ggtree::geom_tippoint(ggplot2::aes(color = color), size = tip_point_size)
        
        # Apply color scale to tree
        p.tree <- switch(tree_color_scale,
          "d3" = p.tree + ggsci::scale_color_d3(),
          "aaas" = p.tree + ggsci::scale_color_aaas(),
          "npg" = p.tree + ggsci::scale_color_npg(),
          "lancet" = p.tree + ggsci::scale_color_lancet(),
          "jco" = p.tree + ggsci::scale_color_jco(),
          "ucscgb" = p.tree + ggsci::scale_color_ucscgb(),
          "uchicago" = p.tree + ggsci::scale_color_uchicago(),
          "simpsons" = p.tree + ggsci::scale_color_simpsons(),
          "rickandmorty" = p.tree + ggsci::scale_color_rickandmorty(),
          {
            warning("Unknown tree color scale: ", tree_color_scale, ". Using default.")
            p.tree + ggplot2::scale_color_discrete()
          }
        )
      }
    }
    
    # Combine plots using aplot
    if (verbose) message("Combining tree and admixture plots...")
    combined_plot <- p.admixture %>% 
      aplot::insert_left(p.tree, width = tree_width)
  }
  
  # Return list with all three plot objects
  if (verbose) message("Analysis completed! Returning plot objects.")
  
  return(list(
    combined_plot = combined_plot,    # aplot object (NULL if no tree)
    admixture_plot = p.admixture,    # ggplot2 object
    tree_plot = p.tree               # ggtree object (NULL if no tree)
  ))
}