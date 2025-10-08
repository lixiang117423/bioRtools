#' Plot protein domain (from Pfam and other tools) distribution with optional phylogenetic tree
#'
#' @param data Character string, dataframe of input data with four columns, id, length, domain, start, end
#' @param tree_file Character string, path to phylogenetic tree file (optional, default is NULL)
#' @param tile_height Numeric, height of geom_tile, default is 0.8
#' @param fill_scale Fill color setting, can be a function or color vector, default is bioRtools::scale_fill_cell()
#' @param tree_xlim_factor Numeric, factor to adjust phylogenetic tree x-axis range, default is 1.5
#' @param tree_width Numeric, width ratio of phylogenetic tree in combined plot, default is 1
#' @param tree_layout Character string, phylogenetic tree layout type, default is "ellipse"
#' @param rounded_corners Logical, whether to use rounded corners for domain rectangles, default is FALSE
#' @param corner_radius Numeric, radius of rounded corners (only used when rounded_corners = TRUE), default is 0.1
#'
#' @return ggplot object or aplot combined plot object
#' @note When rounded_corners = TRUE, the ggforce package is required
#' @export
#'
#' @examples
#' # Without phylogenetic tree
#' plot_pfam(data = "path/to/pfam.tsv")
#' 
#' # With phylogenetic tree
#' plot_pfam(data = "path/to/pfam.tsv", 
#'           tree_file = "path/to/tree.nwk")
#' 
#' # Custom parameters with rounded corners
#' plot_pfam(data = "path/to/pfam.tsv", 
#'           tree_file = "path/to/tree.nwk",
#'           tile_height = 0.6,
#'           tree_xlim_factor = 2.0,
#'           tree_width = 1.5,
#'           rounded_corners = TRUE,
#'           corner_radius = 0.15)

plot_pfam <- function(data, 
                     tree_file = NULL,
                     tile_height = 0.8,
                     fill_scale = bioRtools::scale_fill_cell(),
                     tree_xlim_factor = 1.5,
                     tree_width = 1,
                     tree_layout = "ellipse",
                     rounded_corners = FALSE,
                     corner_radius = 0.1) {
  
  # Load required packages
  library(readr)
  library(dplyr)
  library(magrittr)
  library(ggplot2)
  
  # Load ggforce if rounded corners are requested
  if (rounded_corners) {
    library(ggforce)
  }
  
  # Read pfam data
  # pfam_data <- readr::read_delim(data, col_names = FALSE) %>% 
  #   dplyr::select(1, 3, 6:8) %>% 
  #   magrittr::set_names(c("id", "length", "domain", "start", "end"))

  pfam_data <- data

  # If no phylogenetic tree file provided, plot domain diagram only
  if (is.null(tree_file)) {
    p <- pfam_data %>% 
      dplyr::mutate(y = as.factor(id) %>% as.numeric()) %>% 
      ggplot(aes(x = 0, xend = length, y = id, yend = id)) +
      geom_segment(size = 0.6)
    
    # Add rectangles (rounded or regular)
    if (rounded_corners) {
      p <- p + geom_rounded_rect(aes(xmin = start, xmax = end, 
                                    ymin = as.numeric(as.factor(id)) - tile_height/2, 
                                    ymax = as.numeric(as.factor(id)) + tile_height/2, 
                                    fill = domain), 
                                radius = unit(corner_radius, "npc"))
    } else {
      p <- p + geom_tile(aes(x = (start + end)/2, y = id, width = end - start, height = tile_height, fill = domain))
    }
    
    p <- p + labs(x = "Position") +
      fill_scale +
      theme(legend.title = element_blank(),
            axis.title.y = element_blank())
    
    return(p)
  }
  
  # If phylogenetic tree file is provided
  else {
    # Load phylogenetic tree related packages
    library(ggtree)
    library(aplot)
    
    # Read phylogenetic tree
    tree <- ggtree::read.tree(tree_file)
    
    # Get y coordinate information from tree
    tree_y <- ggtree::ggtree(tree)[["data"]] %>% 
      dplyr::filter(isTip == "TRUE") %>%
      dplyr::select("label", "y") %>% 
      dplyr::rename(id = label)
    
    # Plot protein domain diagram
    p_domain_base <- pfam_data %>% 
      dplyr::left_join(tree_y, by = "id") %>% 
      ggplot(aes(x = 0, xend = length, y = id, yend = id)) +
      geom_segment(size = 0.6)
    
    # Add rectangles (rounded or regular)
    if (rounded_corners) {
      p_domain <- p_domain_base + 
        geom_rounded_rect(aes(xmin = start, xmax = end, 
                             ymin = y - tile_height/2, 
                             ymax = y + tile_height/2, 
                             fill = domain), 
                         radius = unit(corner_radius, "npc"))
    } else {
      p_domain <- p_domain_base + 
        geom_tile(aes(x = (start + end)/2, y = y, width = end - start, height = tile_height, fill = domain))
    }
    
    p_domain <- p_domain +
      labs(x = "Position") +
      fill_scale +
      theme(legend.title = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank())
    
    # Plot phylogenetic tree
    p_tree <- ggtree::ggtree(tree, layout = tree_layout) +
      ggtree::geom_tiplab(align = TRUE)
    
    # Adjust phylogenetic tree x-axis range
    p_tree_final <- p_tree + xlim(NA, max(p_tree$data$x) * tree_xlim_factor)
    
    # Combine plots
    p_final <- aplot::insert_left(p_domain, p_tree_final, width = tree_width)
    
    return(p_final)
  }
}

# Optional convenience function
#' Quick function for plotting protein domains
#'
#' @param data Path to pfam file
#' @param tree_file Path to tree file (optional)
#' @param ... Other parameters passed to plot_pfam
#'
#' @return ggplot object
#' @export

quick_pfam_plot <- function(data, tree_file = NULL, ...) {
  plot_pfam(data = data, 
            tree_file = tree_file, 
            ...)
}