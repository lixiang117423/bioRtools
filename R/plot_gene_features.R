#' Plot Gene Features from Simple Format
#'
#' Visualize gene structure from a simple tab-delimited format containing gene
#' IDs, feature types, start and end positions. Each gene has a 'gene' row
#' defining its span, and multiple feature rows (exon, CDS, etc.) showing
#' features within the gene.
#'
#' @param data Path to tab-delimited file with columns: id, type, start, end,
#'   or a data frame with the same structure
#' @param gene_color Color for gene lines (default: "black")
#' @param base_size Base font size for the plot (default: 12)
#' @param feature_colors Named vector for feature type colors.
#'   Names should match feature types (e.g., exon, CDS, UTR).
#'   Default is c(exon = "#4DAF4A", CDS = "#377EB8")
#' @param max_genes Maximum number of genes to display (default: NULL for all)
#'
#' @return A ggplot object showing gene structure with features
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export?
#'
#' @examples
#' \dontrun{
#' library(bioRtools)
#'
#' # From file
#' p1 <- plot_gene_features("genes.txt")
#' print(p1)
#'
#' # From data frame
#' p2 <- plot_gene_features(df.gene, max_genes = 20)
#' print(p2)
#'
#' # Custom colors
#' p3 <- plot_gene_features(
#'   data = df.gene,
#'   feature_colors = c(exon = "#E41A1C", CDS = "#377EB8")
#' )
#' print(p3)
#' }
plot_gene_features <- function(data,
                               gene_color = "black",
                               base_size = 12,
                               feature_colors = c(exon = "#4DAF4A",
                                                  CDS = "#377EB8"),
                               max_genes = NULL) {

  # Input validation and reading
  if (missing(data) || is.null(data)) {
    stop("'data' is required (file path or data frame)")
  }

  # Read data from file or use provided data frame
  if (is.character(data)) {
    if (!file.exists(data)) {
      stop(paste("File not found:", data))
    }
    df <- tryCatch({
      data.table::fread(data, header = TRUE, verbose = FALSE)
    }, error = function(e) {
      stop(paste("Failed to read file:", e$message))
    })
  } else if (is.data.frame(data)) {
    df <- data
  } else {
    stop("'data' must be a file path (character) or a data frame")
  }

  # Validate and detect column names
  cols <- names(df)

  # Detect ID column
  id_col <- intersect(cols, c("id", "gene_id", "gene", "ID"))[1]
  if (is.na(id_col)) {
    stop(
      "Could not find ID column. Expected: id, gene_id, gene, or ID. ",
      "Found: ", paste(cols, collapse = ", ")
    )
  }

  # Detect type/feature column
  type_col <- intersect(cols, c("type", "feature", "X3", "feature_type"))[1]
  if (is.na(type_col)) {
    stop(
      "Could not find type/feature column. Expected: type, feature, X3, or feature_type. ",
      "Found: ", paste(cols, collapse = ", ")
    )
  }

  # Detect start column
  start_col <- intersect(cols, c("start", "X4", "Start"))[1]
  if (is.na(start_col)) {
    stop(
      "Could not find start column. Expected: start, X4, or Start. ",
      "Found: ", paste(cols, collapse = ", ")
    )
  }

  # Detect end column
  end_col <- intersect(cols, c("end", "X5", "End"))[1]
  if (is.na(end_col)) {
    stop(
      "Could not find end column. Expected: end, X5, or End. ",
      "Found: ", paste(cols, collapse = ", ")
    )
  }

  # Rename columns for clarity
  df <- df %>%
    dplyr::rename(
      gene_id = !!rlang::sym(id_col),
      feature_type = !!rlang::sym(type_col),
      start = !!rlang::sym(start_col),
      end = !!rlang::sym(end_col)
    )

  # Validate that coordinates are numeric
  if (!is.numeric(df$start) || !is.numeric(df$end)) {
    stop("Coordinate columns (X4, X5) must be numeric")
  }

  # Get gene spans from 'gene' rows, or calculate from feature coordinates
  gene_rows <- df %>%
    dplyr::filter(.data$feature_type == "gene")

  if (nrow(gene_rows) > 0) {
    # Use explicit gene rows if available
    gene_spans <- gene_rows %>%
      dplyr::select(.data$gene_id, .data$start, .data$end) %>%
      dplyr::mutate(gene_length = .data$end - .data$start)
  } else {
    # Calculate gene spans from feature coordinates
    gene_spans <- df %>%
      dplyr::group_by(.data$gene_id) %>%
      dplyr::summarise(
        start = min(.data$start, na.rm = TRUE),
        end = max(.data$end, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(gene_length = .data$end - .data$start)
  }

  # Limit genes if max_genes specified
  if (!is.null(max_genes) && nrow(gene_spans) > max_genes) {
    warning(
      "Data contains ", nrow(gene_spans), " genes. ",
      "Showing first ", max_genes, " genes.",
      call. = FALSE
    )
    genes_to_keep <- gene_spans$gene_id[1:max_genes]
    gene_spans <- gene_spans[gene_spans$gene_id %in% genes_to_keep, ]
  }

  # Assign y-coordinates to genes
  gene_spans <- gene_spans %>%
    dplyr::arrange(.data$gene_id) %>%
    dplyr::mutate(y_num = row_number()) %>%
    dplyr::mutate(y = factor(.data$gene_id, levels = .data$gene_id))

  # Calculate relative coordinates for features
  df_plot <- df %>%
    dplyr::filter(.data$feature_type != "gene") %>%
    dplyr::left_join(
      gene_spans %>% dplyr::select(gene_id, gene_start = start, y_num),
      by = "gene_id"
    ) %>%
    dplyr::mutate(
      rel_start = .data$start - .data$gene_start,
      rel_end = .data$end - .data$gene_start
    ) %>%
    dplyr::arrange(.data$y_num, .data$rel_start)

  # Create plot
  p <- ggplot2::ggplot()

  # Add gene lines
  p <- p +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = 0,
        xend = .data$gene_length,
        y = .data$y_num,
        yend = .data$y_num
      ),
      data = gene_spans,
      size = 0.8,
      color = gene_color,
      inherit.aes = FALSE
    )

  # Add feature rectangles
  p <- p +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = .data$rel_start,
        xmax = .data$rel_end,
        ymin = .data$y_num - 0.4,
        ymax = .data$y_num + 0.4,
        fill = .data$feature_type
      ),
      data = df_plot,
      color = NA,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = feature_colors,
      na.value = "grey70"
    ) +
    ggplot2::scale_x_continuous(
      expand = c(0, 0),
      labels = scales::comma_format()
    ) +
    ggplot2::scale_y_continuous(
      breaks = gene_spans$y_num,
      labels = gene_spans$y,
      expand = c(0.02, 0)
    ) +
    ggplot2::labs(
      x = "Position (bp)",
      y = "Gene"
    ) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      legend.title = ggplot2::element_text(face = "bold"),
      legend.position = "top"
    )

  p
}


#' Plot Gene Features with Gene Labels (Alias)
#'
#' Alias for plot_gene_features for backwards compatibility
#'
#' @param data Path to tab-delimited file with columns: id, type, start, end,
#'   or a data frame with the same structure
#' @param gene_color Color for gene lines (default: "black")
#' @param base_size Base font size for the plot (default: 12)
#' @param feature_colors Named vector for feature type colors
#' @param max_genes Maximum number of genes to display (prevents overcrowding)
#'
#' @return A ggplot object showing gene structure with gene labels
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export?
#'
#' @examples
#' \dontrun{
#' library(bioRtools)
#'
#' # Show first 20 genes with labels
#' p1 <- plot_gene_features_labeled("genes.txt", max_genes = 20)
#' print(p1)
#' }
plot_gene_features_labeled <- function(data,
                                      gene_color = "black",
                                      base_size = 12,
                                      feature_colors = c(exon = "#4DAF4A",
                                                         CDS = "#377EB8"),
                                      max_genes = NULL) {
  plot_gene_features(
    data = data,
    gene_color = gene_color,
    base_size = base_size,
    feature_colors = feature_colors,
    max_genes = max_genes
  )
}
