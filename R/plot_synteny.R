#' Create Gene Synteny Plot Using ggplot2
#'
#' @description
#' Create a publication-quality gene synteny plot to visualize conserved
#' gene order and collinearity between multiple species/genomes. This function
#' uses gggenes and ggforce packages to draw genes as arrows and connect
#' homologous genes with Bezier curves.
#'
#' @param gene_data Data frame containing gene information. Must include columns:
#'   \itemize{
#'     \item start: Gene start position (numeric)
#'     \item end: Gene end position (numeric)
#'     \item y: Y-axis position for the gene (numeric, typically 1, 2, 3, etc.)
#'     \item gene_name: Gene name/identifier (character)
#'     \item species: Species name (character)
#'     \item direction: Gene strand direction (character, "+" or "-")
#'     \item gene_type: Gene category for coloring (character)
#'   }
#' @param syntenic_data Data frame containing synteny relationships between genes.
#'   Must include columns:
#'   \itemize{
#'     \item x: X-coordinate (numeric)
#'     \item y: Y-coordinate (numeric)
#'     \item group: Grouping variable for connecting homologous genes (character/numeric)
#'   }
#'   For ribbon-style links (link_type = "ribbon"), each group must have 4 rows
#'   defining the corners of the link region. For curve-style links (link_type = "curve"),
#'   each group must have 2 rows defining the start and end points.
#'   If NULL, syntenic_data will be automatically generated from synteny_groups parameter.
#'   Default is NULL.
#' @param start_col Character string specifying the column name for gene start position.
#'   Default is "start".
#' @param end_col Character string specifying the column name for gene end position.
#'   Default is "end".
#' @param y_col Character string specifying the column name for y-axis position.
#'   Default is "y".
#' @param gene_name_col Character string specifying the column name for gene names.
#'   Default is "gene_name".
#' @param species_col Character string specifying the column name for species.
#'   Default is "species".
#' @param direction_col Character string specifying the column name for gene direction.
#'   Default is "direction".
#' @param gene_type_col Character string specifying the column name for gene categories.
#'   Default is "gene_type".
#' @param x_col Character string specifying the column name for x-coordinate in syntenic_data.
#'   Default is "x".
#' @param y_synt_col Character string specifying the column name for y-coordinate in syntenic_data.
#'   Default is "y".
#' @param group_col Character string specifying the column name for grouping in syntenic_data.
#'   Default is "group".
#' @param synteny_groups Named list or vector defining homologous gene groups.
#'   Each element contains gene names that should be connected as syntenic.
#'   For example: list(group1 = c("Gene1", "Gene2", "Gene3"), group2 = c("Gene4", "Gene5")).
#'   If provided, syntenic_data will be automatically generated. Ignored if syntenic_data is provided.
#'   Default is NULL.
#' @param syn_group_col Character string specifying the column name in gene_data that
#'   contains synteny group identifiers. If provided and syntenic_data is NULL,
#'   syntenic_data will be automatically generated from this column. This is more convenient
#'   than synteny_groups when group information is already in the gene_data (e.g., from Excel).
#'   Default is "syn_group".
#' @param species_labels Character vector of labels for species on y-axis.
#'   Default is NULL (automatically uses species column values from gene_data).
#' @param fill_colors Character vector of colors for gene types.
#'   Default is c("#e59f01", "#56b4e8", "#009f73", "#0072b1").
#' @param link_color Character string specifying the color for synteny links.
#'   Default is "grey".
#' @param link_alpha Numeric value between 0 and 1 for link transparency.
#'   Default is 0.5.
#' @param link_strength Numeric value controlling the curvature of links.
#'   For ribbon links, controls the "bulge" of the diagonal. For curve links,
#'   controls the curvature. Default is 0.5.
#' @param link_type Character string specifying the style of synteny links.
#'   Options: "ribbon" for filled ribbon-style links (requires 4 points per group)
#'   or "curve" for simple curved lines (requires 2 points per group). Default is "ribbon".
#' @param gene_arrow_height Numeric value for gene arrow head height in mm.
#'   Default is 3.
#' @param gene_arrow_width Numeric value for gene arrow head width in mm.
#'   Default is 1.5.
#' @param show_gene_labels Logical indicating whether to display gene names.
#'   Default is TRUE.
#' @param gene_label_fontface Character string specifying font face for gene labels.
#'   Options: "plain", "bold", "italic", "bold.italic". Default is "italic".
#' @param axis_text_italic Logical indicating whether y-axis species labels should be italic.
#'   Default is TRUE.
#' @param title Character string for plot title. Default is NULL.
#' @param subtitle Character string for plot subtitle. Default is NULL.
#'
#' @return A list containing two components:
#' \describe{
#'   \item{plot.synteny}{A ggplot object displaying the gene synteny plot.}
#'   \item{data.summary}{A summary table showing the number of genes and links
#'     per species.}
#' }
#'
#' @details
#' This function creates a gene synteny plot commonly used in comparative genomics
#' to visualize conserved gene order and collinear relationships between species.
#'
#' The plot consists of:
#' \itemize{
#'   \item Genes displayed as arrows pointing in the direction of transcription
#'   \item Species arranged on separate horizontal tracks
#'   \item Bezier curves connecting homologous genes between species
#'   \item Color coding by gene type/category
#' }
#'
#' The function automatically adjusts gene positions to prevent overlap between
#' species tracks. For custom positioning, users should pre-process the gene_data
#' to manually set y-coordinates.
#'
#' @note
#' This function requires the gggenes and ggforce packages. Install them using:
#' \code{install.packages("gggenes")} and \code{install.packages("ggforce")}.
#'
#' The syntenic_data must be manually prepared to specify which genes should be
#' connected. Each row represents one connection between two genes.
#'
#' @references
#' This function is inspired by the synteny plots commonly seen in comparative
#' genomics papers, particularly in CELL, Nature, and Science journals.
#'
#' @examples
#' library(bioRtools)
#'
#' # Load example data
#' data(df.synteny.gene)
#' data(df.synteny.link)
#'
#' # Create basic synteny plot
#' synteny_result <- plot_synteny(
#'   gene_data = df.synteny.gene,
#'   syntenic_data = df.synteny.link,
#'   species_labels = c("A", "B", "C", "D")
#' )
#'
#' # Display the plot
#' synteny_result$plot.synteny
#'
#' # View summary statistics
#' synteny_result$data.summary
#'
#' # Customized plot with different colors
#' plot_synteny(
#'   gene_data = df.synteny.gene,
#'   syntenic_data = df.synteny.link,
#'   species_labels = c("Sp1", "Sp2", "Sp3", "Sp4"),
#'   fill_colors = c("#FF6B6B", "#4ECDC4", "#45B7D1", "#FFA07A"),
#'   link_color = "grey50",
#'   link_alpha = 0.3,
#'   title = "Gene Synteny Analysis",
#'   subtitle = "Cross-species comparison"
#' )$plot.synteny
#'
#' # Automatic syntenic data generation using synteny_groups
#' plot_synteny(
#'   gene_data = df.synteny.gene,
#'   synteny_groups = list(
#'     group1 = c("Gene03", "Gene06", "Gene09", "Gene12"),
#'     group2 = c("Gene04", "Gene07", "Gene10", "Gene13"),
#'     group3 = c("Gene02", "Gene05", "Gene11", "Gene14")
#'   ),
#'   species_labels = c("A", "B", "C", "D")
#' )$plot.synteny
#'
#' # Easiest method: use syn_group column from gene_data
#' # If your data has a "syn_group" column, simply call:
#' plot_synteny(
#'   gene_data = df.synteny.gene
#'   # syntenic_data is auto-generated from syn_group column
#'   # species_labels are auto-generated from species column
#' )$plot.synteny
#'
#' @export
#' @author Xiang LI \email{lixiang117423@@foxmail.com}
#'
plot_synteny <- function(gene_data,
                         syntenic_data = NULL,
                         start_col = "start",
                         end_col = "end",
                         y_col = "y",
                         gene_name_col = "gene_name",
                         species_col = "species",
                         direction_col = "direction",
                         gene_type_col = "gene_type",
                         x_col = "x",
                         y_synt_col = "y",
                         group_col = "group",
                         synteny_groups = NULL,
                         syn_group_col = "syn_group",
                         species_labels = NULL,
                         fill_colors = c("#e59f01", "#56b4e8", "#009f73", "#0072b1"),
                         link_color = "grey",
                         link_alpha = 0.5,
                         link_strength = 0.5,
                         link_type = "ribbon",
                         gene_arrow_height = 3,
                         gene_arrow_width = 1.5,
                         show_gene_labels = TRUE,
                         gene_label_fontface = "italic",
                         axis_text_italic = TRUE,
                         title = NULL,
                         subtitle = NULL) {
  # Auto-generate syntenic_data from synteny_groups or syn_group_col if syntenic_data is NULL
  if (is.null(syntenic_data)) {
    # First check if syn_group_col exists in gene_data
    if (syn_group_col %in% colnames(gene_data) && !all(is.na(gene_data[[syn_group_col]]))) {
      # Generate from syn_group column
      syntenic_data <- auto_generate_synteny_links_from_col(
        gene_data = gene_data,
        syn_group_col = syn_group_col,
        start_col = start_col,
        end_col = end_col,
        y_col = y_col
      )
    } else if (!is.null(synteny_groups)) {
      # Generate from synteny_groups list
      syntenic_data <- auto_generate_synteny_links(
        gene_data = gene_data,
        synteny_groups = synteny_groups,
        gene_name_col = gene_name_col,
        start_col = start_col,
        end_col = end_col,
        y_col = y_col
      )
    } else {
      stop("Either 'syntenic_data', 'synteny_groups', or a gene_data column specified by 'syn_group_col' must be provided",
        call. = FALSE)
    }
  }

  # Input validation
  validate_synteny_inputs(
    gene_data = gene_data,
    syntenic_data = syntenic_data,
    start_col = start_col,
    end_col = end_col,
    y_col = y_col,
    gene_name_col = gene_name_col,
    species_col = species_col,
    direction_col = direction_col,
    gene_type_col = gene_type_col,
    x_col = x_col,
    y_synt_col = y_synt_col,
    group_col = group_col,
    link_type = link_type,
    fill_colors = fill_colors,
    link_alpha = link_alpha,
    link_strength = link_strength,
    gene_arrow_height = gene_arrow_height,
    gene_arrow_width = gene_arrow_width,
    show_gene_labels = show_gene_labels,
    axis_text_italic = axis_text_italic
  )

  # Process gene data
  gene_data_processed <- process_gene_data(
    gene_data = gene_data,
    start_col = start_col,
    end_col = end_col,
    y_col = y_col,
    species_col = species_col
  )

  # Process synteny link data
  syntenic_data_processed <- process_syntenic_data(
    syntenic_data = syntenic_data,
    x_col = x_col,
    y_synt_col = y_synt_col,
    group_col = group_col
  )

  # Prepare segment data for background lines
  segment_data <- prepare_segment_data(
    gene_data = gene_data_processed,
    start_col = start_col,
    end_col = end_col,
    y_col = y_col
  )

  # Generate species labels if not provided
  if (is.null(species_labels)) {
    # Use species column values as labels (ordered by first appearance in y)
    species_labels <- unique(gene_data_processed[[species_col]][order(gene_data_processed[[y_col]])])
  }

  # Create the synteny plot
  plot_synteny <- create_synteny_plot(
    gene_data = gene_data_processed,
    syntenic_data = syntenic_data_processed,
    segment_data = segment_data,
    start_col = start_col,
    end_col = end_col,
    y_col = y_col,
    gene_name_col = gene_name_col,
    species_col = species_col,
    direction_col = direction_col,
    gene_type_col = gene_type_col,
    x_col = x_col,
    y_synt_col = y_synt_col,
    group_col = group_col,
    link_type = link_type,
    species_labels = species_labels,
    fill_colors = fill_colors,
    link_color = link_color,
    link_alpha = link_alpha,
    link_strength = link_strength,
    gene_arrow_height = gene_arrow_height,
    gene_arrow_width = gene_arrow_width,
    show_gene_labels = show_gene_labels,
    gene_label_fontface = gene_label_fontface,
    axis_text_italic = axis_text_italic,
    title = title,
    subtitle = subtitle
  )

  # Create summary statistics
  data_summary <- create_synteny_summary(
    gene_data = gene_data_processed,
    syntenic_data = syntenic_data_processed,
    species_col = species_col,
    group_col = group_col
  )

  # Return results following project conventions
  return(list(
    plot.synteny = plot_synteny,
    data.summary = data_summary
  ))
}


#' Validate Synteny Plot Inputs
#'
#' @param gene_data Gene data frame
#' @param syntenic_data Synteny link data frame
#' @param ... Various column name and parameter specifications
#' @keywords internal
#' @noRd
validate_synteny_inputs <- function(gene_data, syntenic_data, ...) {
  # Check data frames
  if (!is.data.frame(gene_data)) {
    stop("'gene_data' must be a data frame", call. = FALSE)
  }

  if (!is.data.frame(syntenic_data)) {
    stop("'syntenic_data' must be a data frame", call. = FALSE)
  }

  # Extract parameters
  args <- list(...)

  # Check required columns in gene_data
  required_gene_cols <- c(
    args$start_col, args$end_col, args$y_col,
    args$gene_name_col, args$species_col,
    args$direction_col, args$gene_type_col
  )

  missing_gene_cols <- required_gene_cols[!required_gene_cols %in% colnames(gene_data)]
  if (length(missing_gene_cols) > 0) {
    stop(
      "Required columns missing from gene_data: ",
      paste(missing_gene_cols, collapse = ", "),
      call. = FALSE
    )
  }

  # Check required columns in syntenic_data
  required_synt_cols <- c(args$x_col, args$y_synt_col, args$group_col)
  missing_synt_cols <- required_synt_cols[!required_synt_cols %in% colnames(syntenic_data)]
  if (length(missing_synt_cols) > 0) {
    stop(
      "Required columns missing from syntenic_data: ",
      paste(missing_synt_cols, collapse = ", "),
      call. = FALSE
    )
  }

  # Validate numeric parameters
  if (!is.numeric(args$link_alpha) || args$link_alpha < 0 || args$link_alpha > 1) {
    stop("'link_alpha' must be a numeric value between 0 and 1", call. = FALSE)
  }

  if (!is.numeric(args$link_strength) || length(args$link_strength) != 1) {
    stop("'link_strength' must be a single numeric value", call. = FALSE)
  }

  if (!is.numeric(args$gene_arrow_height) || args$gene_arrow_height <= 0) {
    stop("'gene_arrow_height' must be a positive numeric value", call. = FALSE)
  }

  if (!is.numeric(args$gene_arrow_width) || args$gene_arrow_width <= 0) {
    stop("'gene_arrow_width' must be a positive numeric value", call. = FALSE)
  }

  if (!is.logical(args$show_gene_labels) || length(args$show_gene_labels) != 1) {
    stop("'show_gene_labels' must be a single logical value (TRUE/FALSE)", call. = FALSE)
  }

  if (!is.logical(args$axis_text_italic) || length(args$axis_text_italic) != 1) {
    stop("'axis_text_italic' must be a single logical value (TRUE/FALSE)", call. = FALSE)
  }

  # Validate link_type
  if (!args$link_type %in% c("ribbon", "curve")) {
    stop("'link_type' must be either 'ribbon' or 'curve'", call. = FALSE)
  }

  # Validate syntenic_data structure based on link_type
  n_rows_per_group <- syntenic_data %>%
    dplyr::group_by(!!rlang::sym(args$group_col)) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop")

  if (args$link_type == "ribbon") {
    if (any(n_rows_per_group$n != 4)) {
      stop("For link_type = 'ribbon', each group in syntenic_data must have exactly 4 rows",
        call. = FALSE)
    }
  } else if (args$link_type == "curve") {
    if (any(n_rows_per_group$n != 2)) {
      stop("For link_type = 'curve', each group in syntenic_data must have exactly 2 rows",
        call. = FALSE)
    }
  }

  # Check if required packages are available
  if (!requireNamespace("gggenes", quietly = TRUE)) {
    stop("Package 'gggenes' is required. Please install it using: install.packages('gggenes')",
      call. = FALSE)
  }

  if (!requireNamespace("ggforce", quietly = TRUE)) {
    stop("Package 'ggforce' is required. Please install it using: install.packages('ggforce')",
      call. = FALSE)
  }

  invisible(TRUE)
}


#' Auto-generate Synteny Links from Gene Groups
#'
#' @param gene_data Gene data frame
#' @param synteny_groups Named list of gene groups
#' @param gene_name_col Column name for gene names
#' @param start_col Column name for start position
#' @param end_col Column name for end position
#' @param y_col Column name for y position
#' @return Syntenic data frame in ribbon format
#' @keywords internal
#' @noRd
auto_generate_synteny_links <- function(gene_data, synteny_groups,
                                        gene_name_col, start_col, end_col, y_col) {
  # Initialize result
  syntenic_data <- NULL

  # Process each synteny group
  for (group_name in names(synteny_groups)) {
    gene_names <- synteny_groups[[group_name]]

    # Filter genes in this group
    group_genes <- gene_data[gene_data[[gene_name_col]] %in% gene_names, ]

    # Skip if not enough genes
    if (nrow(group_genes) < 2) next

    # Sort by y position
    group_genes <- group_genes[order(group_genes[[y_col]]), ]

    # Create connections between adjacent genes
    for (i in 1:(nrow(group_genes) - 1)) {
      j <- i + 1

      # Extract gene i and gene j
      genes_ij <- group_genes[c(i, j), ]

      # Create ribbon data (4 points)
      ribbon_data <- genes_ij[, c(start_col, end_col, y_col)]
      colnames(ribbon_data)[3] <- "y"

      # Reshape to long format
      ribbon_long <- ribbon_data %>%
        tidyr::pivot_longer(
          cols = c(start_col, end_col),
          names_to = "coord_type",
          values_to = "x"
        ) %>%
        dplyr::select(y, x) %>%
        dplyr::mutate(group = paste0(group_name, "_", i, "_", j))

      # Bind to result
      syntenic_data <- dplyr::bind_rows(syntenic_data, ribbon_long)
    }
  }

  return(syntenic_data)
}


#' Auto-generate Synteny Links from syn_group Column
#'
#' @param gene_data Gene data frame with syn_group column
#' @param syn_group_col Column name containing synteny group identifiers
#' @param start_col Column name for start position
#' @param end_col Column name for end position
#' @param y_col Column name for y position
#' @return Syntenic data frame in ribbon format
#' @keywords internal
#' @noRd
auto_generate_synteny_links_from_col <- function(gene_data, syn_group_col,
                                                 start_col, end_col, y_col) {
  # Get unique synteny groups
  syn_groups <- unique(gene_data[[syn_group_col]])
  syn_groups <- syn_groups[!is.na(syn_groups)]

  # Initialize result
  syntenic_data <- NULL

  # Process each synteny group
  for (group_name in syn_groups) {
    # Filter genes in this group
    group_genes <- gene_data[!is.na(gene_data[[syn_group_col]]) &
      gene_data[[syn_group_col]] == group_name, ]

    # Skip if not enough genes
    if (nrow(group_genes) < 2) next

    # Sort by y position
    group_genes <- group_genes[order(group_genes[[y_col]]), ]

    # Create connections between adjacent genes
    for (i in 1:(nrow(group_genes) - 1)) {
      j <- i + 1

      # Extract gene i and gene j
      genes_ij <- group_genes[c(i, j), ]

      # Create ribbon data (4 points)
      ribbon_data <- genes_ij[, c(start_col, end_col, y_col)]
      colnames(ribbon_data)[3] <- "y"

      # Reshape to long format
      ribbon_long <- ribbon_data %>%
        tidyr::pivot_longer(
          cols = c(start_col, end_col),
          names_to = "coord_type",
          values_to = "x"
        ) %>%
        dplyr::select(y, x) %>%
        dplyr::mutate(group = paste0(group_name, "_", i, "_", j))

      # Bind to result
      syntenic_data <- dplyr::bind_rows(syntenic_data, ribbon_long)
    }
  }

  return(syntenic_data)
}


#' Process Gene Data
#'
#' @param gene_data Raw gene data frame
#' @param start_col Column name for start position
#' @param end_col Column name for end position
#' @param y_col Column name for y position
#' @param species_col Column name for species
#' @return Processed gene data frame
#' @keywords internal
#' @noRd
process_gene_data <- function(gene_data, start_col, end_col, y_col, species_col) {
  # Remove rows with NA values in critical columns
  gene_data_clean <- gene_data %>%
    dplyr::filter(
      !is.na(!!rlang::sym(start_col)),
      !is.na(!!rlang::sym(end_col)),
      !is.na(!!rlang::sym(y_col))
    )

  if (nrow(gene_data_clean) == 0) {
    stop("No valid gene data remaining after filtering NA values", call. = FALSE)
  }

  return(gene_data_clean)
}


#' Process Syntenic Data
#'
#' @param syntenic_data Raw synteny data frame
#' @param x_col Column name for x coordinate
#' @param y_synt_col Column name for y coordinate
#' @param group_col Column name for group
#' @return Processed synteny data frame
#' @keywords internal
#' @noRd
process_syntenic_data <- function(syntenic_data, x_col, y_synt_col, group_col) {
  # Remove rows with NA values
  syntenic_data_clean <- syntenic_data %>%
    dplyr::filter(
      !is.na(!!rlang::sym(x_col)),
      !is.na(!!rlang::sym(y_synt_col)),
      !is.na(!!rlang::sym(group_col))
    )

  if (nrow(syntenic_data_clean) == 0) {
    stop("No valid synteny data remaining after filtering NA values", call. = FALSE)
  }

  return(syntenic_data_clean)
}


#' Prepare Segment Data for Background Lines
#'
#' @param gene_data Gene data frame
#' @param start_col Column name for start position
#' @param end_col Column name for end position
#' @param y_col Column name for y position
#' @return Data frame with segment data
#' @keywords internal
#' @noRd
prepare_segment_data <- function(gene_data, start_col, end_col, y_col) {
  segment_df <- gene_data %>%
    dplyr::group_by(!!rlang::sym(y_col)) %>%
    dplyr::summarise(
      min = min(!!rlang::sym(start_col)),
      max = max(!!rlang::sym(end_col)),
      y = min(!!rlang::sym(y_col)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      min = min - 3,
      max = max + 3
    )

  return(segment_df)
}


#' Create Synteny Plot
#'
#' @param gene_data Processed gene data frame
#' @param syntenic_data Processed synteny data frame
#' @param segment_data Segment data for background lines
#' @param ... Various plotting parameters
#' @return A ggplot object
#' @keywords internal
#' @noRd
create_synteny_plot <- function(gene_data, syntenic_data, segment_data, ...) {
  args <- list(...)

  # Create base plot
  p <- ggplot2::ggplot() +
    # Add background segment lines
    ggplot2::geom_segment(
      data = segment_data,
      ggplot2::aes(
        x = min,
        xend = max,
        y = y,
        yend = y
      ),
      color = "grey",
      linewidth = 0.5
    )

  # Add synteny links based on link_type
  if (args$link_type == "ribbon") {
    # Use geom_diagonal_wide for ribbon-style links
    p <- p +
      ggforce::geom_diagonal_wide(
        data = syntenic_data,
        ggplot2::aes(
          x = !!rlang::sym(args$x_col),
          y = !!rlang::sym(args$y_synt_col),
          group = !!rlang::sym(args$group_col)
        ),
        fill = args$link_color,
        alpha = args$link_alpha,
        strength = args$link_strength,
        radius = grid::unit(0, "mm"),
        orientation = "y"
      )
  } else if (args$link_type == "curve") {
    # Use geom_curve for simple curved lines
    p <- p +
      ggplot2::geom_curve(
        data = syntenic_data %>%
          dplyr::group_by(!!rlang::sym(args$group_col)) %>%
          dplyr::summarise(
            x_start = .data[[args$x_col]][1],
            y_start = .data[[args$y_synt_col]][1],
            x_end = .data[[args$x_col]][2],
            y_end = .data[[args$y_synt_col]][2],
            .groups = "drop"
          ),
        ggplot2::aes(
          x = x_start,
          y = y_start,
          xend = x_end,
          yend = y_end
        ),
        color = args$link_color,
        alpha = args$link_alpha,
        curvature = args$link_strength,
        linewidth = 0.5,
        ncp = 20
      )
  }

  # Continue building the plot
  p <- p +
    # Add gene arrows
    gggenes::geom_gene_arrow(
      data = gene_data,
      ggplot2::aes(
        xmin = !!rlang::sym(args$start_col),
        xmax = !!rlang::sym(args$end_col),
        y = !!rlang::sym(args$y_col),
        fill = !!rlang::sym(args$gene_type_col),
        forward = !!rlang::sym(args$direction_col) == "+"
      ),
      arrowhead_height = grid::unit(args$gene_arrow_height, "mm"),
      arrowhead_width = grid::unit(args$gene_arrow_width, "mm"),
      alpha = 1
    )

  # Add gene labels if requested
  if (args$show_gene_labels) {
    p <- p +
      ggplot2::geom_text(
        data = gene_data,
        ggplot2::aes(
          x = !!rlang::sym(args$start_col),
          y = !!rlang::sym(args$y_col),
          label = !!rlang::sym(args$gene_name_col)
        ),
        hjust = 0,
        vjust = 2,
        fontface = args$gene_label_fontface,
        size = 3
      )
  }

  # Apply theme and scales
  p <- p +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      legend.position = "none",
      axis.text.y = ggplot2::element_text(
        face = ifelse(args$axis_text_italic, "italic", "plain")
      ),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 12)
    ) +
    ggplot2::scale_y_continuous(
      breaks = unique(gene_data[[args$y_col]]),
      labels = args$species_labels
    ) +
    ggplot2::scale_fill_manual(values = args$fill_colors) +
    ggplot2::labs(
      title = args$title,
      subtitle = args$subtitle
    )

  return(p)
}


#' Create Synteny Summary Statistics
#'
#' @param gene_data Gene data frame
#' @param syntenic_data Synteny data frame
#' @param species_col Species column name
#' @param group_col Group column name
#' @return Summary data frame
#' @keywords internal
#' @noRd
create_synteny_summary <- function(gene_data, syntenic_data, species_col, group_col) {
  summary <- gene_data %>%
    dplyr::group_by(!!rlang::sym(species_col)) %>%
    dplyr::summarise(
      n_genes = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::arrange(!!rlang::sym(species_col))

  # Add total number of links
  total_links <- length(unique(syntenic_data[[group_col]]))

  summary$total_links <- total_links

  return(summary)
}
