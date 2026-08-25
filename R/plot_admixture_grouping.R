#' Plot Admixture Grouping Bar Ordered by Tree
#'
#' Draw a single-K population grouping bar chart: each sample gets a
#' height-1 bar coloured by its dominant Admixture cluster at a given K,
#' with samples arranged left to right following the tree tip order (top
#' to bottom of the drawn tree). Uses \code{scale_fill_bioinfo()} for the
#' cluster palette and \code{theme_void()} for a clean structure-plot
#' look.
#'
#' @param tree Path to a Newick tree file, or a \code{phylo} object. The
#'   tree is drawn with \code{ggtree::ggtree()} and the sample order is
#'   taken from the drawn tips via \code{ggtree::get_taxa_name()}.
#' @param data Path to a delimited file readable by \code{\link{read_data}},
#'   or a data frame containing an ID column and one cluster-assignment
#'   column per K value (e.g. \code{K2}, \code{K3}, ...). Each \code{Kx}
#'   column holds the dominant cluster number for each sample at K = x.
#' @param k K value to plot (default: 8). A numeric \code{k} selects the
#'   column \code{K8}; a character string is used as the column name
#'   directly.
#' @param id_col Column of \code{data} holding the sample IDs
#'   (default: \code{"FID"}).
#'
#' @param show_legend Logical; whether to draw the cluster legend
#'   (default: \code{FALSE}).
#' @return A ggplot2 object.
#'
#' @details
#' Samples run left to right in the same order as tree tips top to bottom
#' (the reversed tip order of the drawn tree). Samples present in the tree
#' but missing from \code{data} keep an empty slot on the axis (with a
#' warning); samples present in \code{data} but missing from the tree are
#' dropped from the plot (with a warning).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' p <- plot_admixture_grouping(
#'   tree = "path/to/tree.nwk",
#'   data = "path/to/admixture_grouping.csv",
#'   k = 8
#' )
#' p
#' }
plot_admixture_grouping <- function(tree, data, k = 8, id_col = "FID",
                                    show_legend = FALSE) {
  # 1. Input validation -------------------------------------------------
  if (inherits(tree, "phylo")) {
    phylo_tree <- tree
  } else if (is.character(tree) && length(tree) == 1L) {
    phylo_tree <- read_data(tree)
    if (!inherits(phylo_tree, "phylo")) {
      stop("'tree' file did not parse as a phylo object: ", tree)
    }
  } else {
    stop("'tree' must be a phylo object or a path to a Newick file")
  }

  if (is.data.frame(data)) {
    df <- data
  } else if (is.character(data) && length(data) == 1L) {
    df <- read_data(data)
  } else {
    stop("'data' must be a data frame or a path to a delimited file")
  }

  if (!is.character(id_col) || length(id_col) != 1L || !nzchar(id_col)) {
    stop("'id_col' must be a single column name")
  }
  if (!id_col %in% colnames(df)) {
    stop("'id_col' (", id_col, ") not found in data columns: ",
         paste(colnames(df), collapse = ", "))
  }

  if (is.numeric(k) && length(k) == 1L && !is.na(k) && k >= 1 && k == floor(k)) {
    cluster_col <- paste0("K", k)
  } else if (is.character(k) && length(k) == 1L && nzchar(k)) {
    cluster_col <- k
  } else {
    stop("'k' must be a single positive integer (selects the column K<k>) ",
         "or a column name")
  }
  if (!cluster_col %in% colnames(df)) {
    stop("Cluster column '", cluster_col, "' not found in data. Available ",
         "K columns: ",
         paste(grep("^K[0-9]+$", colnames(df), value = TRUE), collapse = ", "))
  }

  if (!is.logical(show_legend) || length(show_legend) != 1L || is.na(show_legend)) {
    stop("'show_legend' must be a single logical value")
  }

  # 2. Tip order of the drawn tree (bottom to top) -----------------------
  label_order <- phylo_tree %>%
    ggtree::ggtree() %>%
    ggtree::get_taxa_name()

  # 3. Build the plotting data -------------------------------------------
  plot_df <- df %>%
    dplyr::select(sample = dplyr::all_of(id_col),
                  cluster = dplyr::all_of(cluster_col))

  if (anyNA(plot_df$cluster)) {
    warning(sum(is.na(plot_df$cluster)),
            " sample(s) have no cluster value in column '", cluster_col,
            "'; dropped from the plot")
    plot_df <- plot_df %>%
      dplyr::filter(!is.na(.data$cluster))
  }

  # cluster labels like "1".."20" must sort numerically, not lexicographically
  lv <- unique(as.character(plot_df$cluster))
  lv_num <- suppressWarnings(as.numeric(lv))
  if (all(!is.na(lv_num))) {
    lv <- as.character(sort(lv_num))
  } else {
    lv <- sort(lv)
  }
  plot_df$cluster <- factor(as.character(plot_df$cluster), levels = lv)

  # 4. Match samples with the tree ---------------------------------------
  matched <- intersect(label_order, plot_df$sample)
  only_tree <- setdiff(label_order, plot_df$sample)
  only_data <- setdiff(plot_df$sample, label_order)

  if (length(matched) == 0L) {
    stop("No samples match between tree and data.\n",
         "  Tree samples (first 10): ",
         paste(utils::head(label_order, 10), collapse = ", "), "\n",
         "  Data samples (first 10): ",
         paste(utils::head(unique(plot_df$sample), 10), collapse = ", "))
  }
  if (length(only_data) > 0L) {
    warning(length(only_data),
            " sample(s) in data but not in the tree; dropped from the plot: ",
            paste(utils::head(only_data, 10), collapse = ", "))
  }
  if (length(only_tree) > 0L) {
    warning(length(only_tree),
            " sample(s) in the tree but not in data; shown as empty slots: ",
            paste(utils::head(only_tree, 10), collapse = ", "))
  }

  # 5. Plot: samples left to right = tree tips top to bottom --------------
  plot_df %>%
    ggplot2::ggplot(ggplot2::aes(.data$sample, y = 1, fill = .data$cluster)) +
    ggplot2::geom_col(show.legend = show_legend) +
    ggplot2::scale_x_discrete(limits = rev(label_order)) +
    scale_fill_bioinfo() +
    ggplot2::theme_void()
}
