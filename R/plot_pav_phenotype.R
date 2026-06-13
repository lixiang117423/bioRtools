#' Plot PAV Heatmap with Phenotype Annotation
#'
#' Creates a side-by-side heatmap: gene/kmer presence-absence variation (PAV)
#' on the left (clustered by sample) and a phenotype (e.g., disease level)
#' annotation heatmap on the right. Samples are ordered by hierarchical
#' clustering on the PAV matrix.
#'
#' @param data_pav Long-format data frame with gene, sample, and PAV value
#'   columns (0 = absent, 1 = present).
#' @param data_phenotype Data frame with sample, group, and phenotype level
#'   columns.
#' @param pav_gene_col Column in \code{data_pav} for gene/kmer identifier
#'   (default: "Start").
#' @param pav_sample_col Column in \code{data_pav} for sample identifier.
#' @param pav_value_col Column in \code{data_pav} for PAV value
#'   (default: "value").
#' @param pheno_sample_col Column in \code{data_phenotype} for sample
#'   identifier (default: "sample").
#' @param pheno_group_col Column in \code{data_phenotype} for group
#'   (displayed on x-axis; default: "group").
#' @param pheno_level_col Column in \code{data_phenotype} for phenotype
#'   level (fill color; default: "level").
#' @param title Plot title for the PAV panel (default: NULL).
#' @param pav_colors Named vector of two colors for PAV values 0 and 1
#'   (default: \code{c("0" = "#4A90D9", "1" = "#D9544A")}).
#' @param cluster_method Hierarchical clustering method passed to
#'   \code{hclust} (default: "complete").
#' @param show_gene_axis Logical; show gene labels on PAV x-axis
#'   (default: FALSE).
#'
#' @return A named list containing:
#'   \describe{
#'     \item{\code{plot}}{Combined plot object (via \code{aplot}).}
#'     \item{\code{sample_order}}{Character vector of samples in clustered
#'       order.}
#'     \item{\code{p_pav}}{PAV heatmap as standalone ggplot object.}
#'     \item{\code{p_phenotype}}{Phenotype heatmap as standalone ggplot
#'       object.}
#'   }
#'
#' @details
#' Samples are clustered using Euclidean distance + hierarchical clustering
#' on the binary PAV matrix. The resulting dendrogram order is used to align
#' both the PAV and phenotype panels. The two panels are combined using
#' \code{aplot::insert_right()}.
#'
#' @examples
#' \dontrun{
#' library(bioRtools)
#'
#' # PAV data in long format
#' df.pav <- data.frame(
#'   Start = rep(paste0("kmer_", 1:10), each = 20),
#'   sample = rep(paste0("S", 1:20), 10),
#'   value = sample(0:1, 200, replace = TRUE)
#' )
#'
#' # Phenotype data
#' df.pheno <- data.frame(
#'   sample = paste0("S", 1:20),
#'   level = sample(1:5, 20, replace = TRUE),
#'   group = sample(c("G1", "G2", "G3"), 20, replace = TRUE)
#' )
#'
#' res <- bioRtools::plot_pav_phenotype(df.pav, df.pheno,
#'   pav_sample_col = "sample",
#'   title = "Rps1-k"
#' )
#' print(res$plot)
#' }
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
plot_pav_phenotype <- function(data_pav, data_phenotype,
                               pav_gene_col = "Start",
                               pav_sample_col,
                               pav_value_col = "value",
                               pheno_sample_col = "sample",
                               pheno_group_col = "group",
                               pheno_level_col = "level",
                               title = NULL,
                               pav_colors = c("0" = "#4A90D9", "1" = "#D9544A"),
                               cluster_method = "complete",
                               show_gene_axis = FALSE) {

  # --- Input validation ---
  if (!is.data.frame(data_pav)) stop("'data_pav' must be a data frame")
  if (!is.data.frame(data_phenotype)) stop("'data_phenotype' must be a data frame")
  if (missing(pav_sample_col)) stop("'pav_sample_col' is required")

  for (col in c(pav_gene_col, pav_sample_col, pav_value_col)) {
    if (!col %in% names(data_pav)) {
      stop("Column '", col, "' not found in data_pav")
    }
  }
  for (col in c(pheno_sample_col, pheno_group_col, pheno_level_col)) {
    if (!col %in% names(data_phenotype)) {
      stop("Column '", col, "' not found in data_phenotype")
    }
  }

  # --- Prepare PAV matrix for clustering ---
  df_pav <- data_pav %>%
    dplyr::select(
      gene = !!rlang::sym(pav_gene_col),
      sample = !!rlang::sym(pav_sample_col),
      value = !!rlang::sym(pav_value_col)
    )

  mat_pav <- df_pav %>%
    tidyr::pivot_wider(names_from = gene, values_from = value, id_cols = sample) %>%
    tibble::column_to_rownames("sample") %>%
    as.matrix()

  # Fill NA with 0
  mat_pav[is.na(mat_pav)] <- 0

  # Cluster samples
  hc <- stats::hclust(stats::dist(mat_pav), method = cluster_method)
  sample_order <- rownames(mat_pav)[hc$order]

  # Filter phenotype data to matching samples
  common_samples <- intersect(sample_order, data_phenotype[[pheno_sample_col]])
  if (length(common_samples) == 0) {
    stop("No matching samples between PAV and phenotype data")
  }
  sample_order <- sample_order[sample_order %in% common_samples]

  # --- PAV heatmap ---
  p_pav <- df_pav %>%
    dplyr::filter(sample %in% common_samples) %>%
    ggplot2::ggplot(ggplot2::aes(
      x = factor(gene),
      y = sample,
      fill = factor(value)
    )) +
    ggplot2::geom_tile() +
    ggplot2::scale_y_discrete(limits = sample_order) +
    ggplot2::scale_fill_manual(values = pav_colors) +
    ggplot2::labs(title = title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      legend.position = "none",
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
    )

  if (!show_gene_axis) {
    p_pav <- p_pav + ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
  }

  # --- Phenotype heatmap ---
  df_pheno <- data_phenotype %>%
    dplyr::filter(!!rlang::sym(pheno_sample_col) %in% common_samples) %>%
    dplyr::mutate(across(dplyr::all_of(pheno_level_col), as.character))

  p_phenotype <- df_pheno %>%
    ggplot2::ggplot(ggplot2::aes(
      x = !!rlang::sym(pheno_group_col),
      y = !!rlang::sym(pheno_sample_col),
      fill = !!rlang::sym(pheno_level_col)
    )) +
    ggplot2::geom_tile() +
    ggplot2::scale_y_discrete(limits = sample_order) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    )

  # Use scale_fill_research if available, otherwise scale_fill_brewer
  tryCatch({
    p_phenotype <- p_phenotype + bioRtools::scale_fill_research()
  }, error = function(e) {
    p_phenotype <<- p_phenotype + ggplot2::scale_fill_brewer(palette = "Set1")
  })

  # --- Combine ---
  p_combined <- aplot::insert_right(p_pav, p_phenotype)

  list(
    plot = p_combined,
    sample_order = sample_order,
    p_pav = p_pav,
    p_phenotype = p_phenotype
  )
}
