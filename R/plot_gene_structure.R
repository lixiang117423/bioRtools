#' Plot Gene Structure from GFF/GFF3 File
#'
#' Visualize gene structure including exons, introns, and other genomic features
#' from GFF/GFF3 format files. This function parses gene annotation files and
#' creates publication-ready gene structure plots with optional phylogenetic tree.
#'
#' @param gff Path to GFF/GFF3 file containing gene annotation data
#' @param tree Either "none" (default) or path to a Newick tree file for
#'   phylogenetic ordering of genes
#' @param exon_color Color for exon features (default: "#4DAF4A")
#' @param cds_color Color for CDS features (default: "#377EB8")
#' @param utr_color Color for UTR features (default: "#FF7F00")
#' @param gene_line_size Size of gene line (default: 0.8)
#' @param gene_height Height of gene features (default: 0.4)
#' @param x_label X-axis label (default: "")
#' @param base_size Base font size for the plot (default: 12)
#'
#' @return A ggplot object displaying gene structure with:
#'   - Horizontal lines representing gene spans
#'   - Rectangles showing exons, CDS, and UTR regions
#'   - Color-coded feature types (exon, CDS, UTR)
#'   - Optional phylogenetic tree layout
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(bioRtools)
#'
#' # Basic usage without tree
#' gff_file <- system.file("extdata", "genes.gff3", package = "bioRtools")
#' p1 <- plot_gene_structure(gff = gff_file)
#' print(p1)
#'
#' # With phylogenetic tree
#' tree_file <- system.file("extdata", "gene_tree.nwk", package = "bioRtools")
#' p2 <- plot_gene_structure(
#'   gff = gff_file,
#'   tree = tree_file,
#'   exon_color = "#E41A1C",
#'   cds_color = "#377EB8"
#' )
#' print(p2)
#'
#' # Save plot
#' ggplot2::ggsave("gene_structure.pdf", p1, width = 10, height = 6)
#' }
plot_gene_structure <- function(gff,
                                 tree = "none",
                                 exon_color = "#4DAF4A",
                                 cds_color = "#377EB8",
                                 utr_color = "#FF7F00",
                                 gene_line_size = 0.8,
                                 gene_height = 0.4,
                                 x_label = "",
                                 base_size = 12) {

  # Input validation
  if (missing(gff) || is.null(gff)) {
    stop("'gff' file path is required")
  }

  if (!is.character(gff)) {
    stop("'gff' must be a character string specifying the file path")
  }

  if (!file.exists(gff)) {
    stop(paste("GFF file not found:", gff))
  }

  if (!is.character(tree)) {
    stop("'tree' must be either 'none' or a file path")
  }

  if (tree != "none" && !file.exists(tree)) {
    stop(paste("Tree file not found:", tree))
  }

  # Validate color parameters
  if (!all(sapply(list(exon_color, cds_color, utr_color), is.character))) {
    stop("Color parameters must be character strings (hex codes or color names)")
  }

  # Validate numeric parameters
  if (!is.numeric(gene_line_size) || gene_line_size <= 0) {
    stop("'gene_line_size' must be a positive number")
  }

  if (!is.numeric(gene_height) || gene_height <= 0) {
    stop("'gene_height' must be a positive number")
  }

  # Read GFF file
  df_gff <- tryCatch({
    data.table::fread(gff, header = FALSE, verbose = FALSE)
  }, error = function(e) {
    stop(paste("Failed to read GFF file:", e$message))
  })

  # Check if GFF file has required columns
  if (ncol(df_gff) < 9) {
    stop("GFF file must have at least 9 columns (standard GFF format)")
  }

  # Extract gene IDs from column 9 (attributes)
  df_gff$gene_id <- ""

  for (i in seq_len(nrow(df_gff))) {
    attributes <- stringr::str_split(df_gff$V9[i], ";")[[1]]

    for (attr in attributes) {
      if (stringr::str_sub(attr, 1, 6) == "Parent") {
        trans_id <- stringr::str_split(attr, "=")[[1]][2]
        df_gff$gene_id[i] <- stringr::str_split(trans_id, "\\.")[[1]][1]
      } else if (stringr::str_sub(attr, 1, 3) == "ID=") {
        # Some GFF files use ID= instead of Parent=
        trans_id <- stringr::str_split(attr, "=")[[1]][2]
        df_gff$gene_id[i] <- stringr::str_split(trans_id, "\\.")[[1]][1]
      }
    }
  }

  # Remove rows without gene IDs
  if (all(df_gff$gene_id == "")) {
    stop("No valid gene IDs found in GFF file. Check the format of column 9.")
  }

  df_gff <- df_gff[df_gff$gene_id != "", ]

  # Calculate gene coordinates
  df_gff <- df_gff %>%
    dplyr::group_by(gene_id) %>%
    dplyr::mutate(
      min_coord = min(min(V4), min(V5)),
      max_coord = max(max(V4), max(V5)),
      .before = 1
    ) %>%
    dplyr::mutate(
      start = ifelse(V7 == "+", V4 - min_coord, abs(V4 - max_coord)),
      end = ifelse(V7 == "+", V5 - min_coord, abs(V5 - max_coord)),
      gene_length = max_coord - min_coord
    ) %>%
    dplyr::ungroup()

  # Determine gene ordering
  if (tree == "none") {
    # Simple alphabetical ordering
    df_order <- data.frame(
      gene_id = unique(df_gff$gene_id),
      y = seq_along(unique(df_gff$gene_id))
    )
  } else {
    # Use phylogenetic tree for ordering
    phylo_tree <- tryCatch({
      ape::read.tree(file = tree)
    }, error = function(e) {
      stop(paste("Failed to read tree file:", e$message))
    })

    tree_plot <- ggtree::ggtree(phylo_tree)

    df_order <- tree_plot$data %>%
      dplyr::filter(isTip == TRUE) %>%
      dplyr::select(label, y) %>%
      dplyr::rename(gene_id = label)
  }

  # Merge gene order with annotation
  df_gff <- df_gff %>%
    dplyr::left_join(df_order, by = "gene_id") %>%
    dplyr::mutate(
      y_min = y - gene_height / 2,
      y_max = y + gene_height / 2
    ) %>%
    dplyr::arrange(y) %>%
    dplyr::mutate(gene_id = factor(gene_id, levels = unique(gene_id)))

  # Create plot
  p <- ggplot2::ggplot()

  # Add gene lines
  if (tree == "none") {
    p <- p +
      ggplot2::geom_segment(
        data = df_gff %>%
          dplyr::distinct(gene_id, start, gene_length, y),
        ggplot2::aes(x = start, xend = gene_length, y = gene_id, yend = gene_id),
        size = gene_line_size,
        color = "black"
      )
  } else {
    p <- p +
      ggplot2::geom_segment(
        data = df_gff %>%
          dplyr::distinct(gene_id, start, gene_length, y),
        ggplot2::aes(x = start, xend = gene_length, y = y, yend = y),
        size = gene_line_size,
        color = "black"
      )
  }

  # Define feature colors
  feature_colors <- c(
    "exon" = exon_color,
    "CDS" = cds_color,
    "five_prime_UTR" = utr_color,
    "three_prime_UTR" = utr_color
  )

  # Add feature rectangles (excluding mRNA)
  df_features <- df_gff %>%
    dplyr::filter(V3 != "mRNA")

  p <- p +
    ggplot2::geom_rect(
      data = df_features,
      ggplot2::aes(
        xmin = start,
        xmax = end,
        ymin = y_min,
        ymax = y_max,
        fill = V3
      ),
      color = NA
    ) +
    ggplot2::scale_fill_manual(
      values = feature_colors,
      na.value = "grey70"
    ) +
    ggplot2::scale_x_continuous(
      expand = c(0, 0),
      breaks = scales::pretty_breaks()(c(0, max(df_gff$gene_length)))
    ) +
    ggplot2::labs(
      x = x_label,
      y = ifelse(tree == "none", "Gene", "")
    ) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      legend.title = ggplot2::element_blank()
    )

  # Adjust theme based on tree presence
  if (tree != "none") {
    p <- p +
      ggplot2::theme(
        axis.ticks.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.line.y = ggplot2::element_blank()
      )
  }

  p
}
