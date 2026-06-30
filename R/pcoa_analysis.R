#' Perform principal coordinate analysis (PCoA).
#'
#' @param data Feature table, with rows representing samples and columns representing feature values, such as OTUs. The transpose (features as rows, samples as columns) is also accepted — see \code{feature_as_row}.
#' @param sample Sample table, with the first column containing sample names. There are no specific requirements for the names of the subsequent columns, but the sample names must match those in the feature table.
#' @param method The method for calculating distances in vegan::vegdist() is defaulted to Bray-Curtis.
#' @param x The principal coordinate used for the X-axis in the plot is defaulted to PCo1.
#' @param y The principal coordinate used for the y-axis in the plot is defaulted to PCo2.
#' @param size The size of point.
#' @param color Column names in the sample data frame used for coloring the points.
#' @param alpha The alpha of point.
#' @param feature_as_row Logical or \code{NA}. \code{NA} (default) auto-detects
#'   the orientation by matching sample IDs from \code{sample} against the row
#'   and column names of \code{data}; \code{TRUE} forces features-as-rows;
#'   \code{FALSE} forces samples-as-rows. When detected or forced, the matrix is
#'   transposed internally so a manual \code{t()} is not needed.
#'
#' @return A list containing four components:
#' \describe{
#'   \item{result_pcoa}{The output from ape::pcoa(), which can be called directly.}
#'   \item{plot_pcoa}{The plotting results, defaulting to PCo1 and PCo2. The output is a ggplot object, which can be fine-tuned using ggplot2.}
#'   \item{point_data}{The data used for plotting, which users can call for their own plots or export for use in other software.}
#'   \item{eigenvalue_pcoa}{The explained variance of the principal coordinates, starting from PCo1 by default.}
#' }
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' library(dplyr)
#' library(bioRtools)
#'
#' data(df.pcoa.otu)
#' data(df.pcoa.sample)
#'
#' pcoa_analysis(data = df.pcoa.otu, sample = df.pcoa.sample) -> pcoa_res
#'
pcoa_analysis <- function(data, sample, method = "bray", x = "PCo1", y = "PCo2", size = 2, color = "group", alpha = 1, feature_as_row = NA) {
  # Resolve orientation before any row-name-based sample matching / vegdist()
  data <- orient_to_sample_row(data, sample, NULL, feature_as_row, FALSE)

  # Auto-detect sample ID column: if no "sample" column, find one matching data rownames
  if (!"sample" %in% colnames(sample)) {
    data_rn <- rownames(data)
    if (!is.null(data_rn)) {
      for (col in names(sample)) {
        if (is.character(sample[[col]]) && all(data_rn %in% sample[[col]])) {
          names(sample)[names(sample) == col] <- "sample"
          break
        }
      }
    }
  }

  # Step 1: Calculate distance matrix and perform PCoA
  data %>%
    vegan::vegdist(method = method) %>%
    ape::pcoa() -> pcoa_res

  # Step 2: Extract eigenvalues and calculate explained variance
  pcoa_res$values %>%
    as.data.frame() %>%
    dplyr::select(2) %>%
    dplyr::mutate(pcoa = paste0("PCo", seq_len(nrow(.)))) %>%
    dplyr::select(pcoa, Relative_eig) %>%
    dplyr::mutate(
      relative_eig = round(Relative_eig * 100, 2),
      percentage = paste0(relative_eig, "%"),
      label = paste0(pcoa, " (", percentage, ")")
    ) -> eigenvalue_pcoa

  # Step 3: Extract coordinate data and merge with sample information
  pcoa_res$vectors %>%
    as.data.frame() %>%
    dplyr::select(1:5) %>%
    {
      df <- .
      pco_labels <- paste0("PCo", seq_len(ncol(df)), " (",
        eigenvalue_pcoa$percentage[seq_len(ncol(df))], ")")
      magrittr::set_names(df, pco_labels)
    } %>%
    tibble::rownames_to_column(var = "sample") %>%
    dplyr::left_join(sample, by = "sample") -> point_data

  # Step 4: Create axis labels with explained variance
  x_label <- eigenvalue_pcoa %>%
    dplyr::filter(pcoa == x) %>%
    dplyr::pull(label)

  y_label <- eigenvalue_pcoa %>%
    dplyr::filter(pcoa == y) %>%
    dplyr::pull(label)

  # Step 5: Create PCoA plot
  plot_pcoa <- point_data %>%
    ggplot2::ggplot(ggplot2::aes(
      x = !!rlang::sym(x),
      y = !!rlang::sym(y),
      color = !!rlang::sym(color)
    )) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "#222222",
      linetype = "dashed",
      linewidth = 0.5
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      color = "#222222",
      linetype = "dashed",
      linewidth = 0.5
    ) +
    ggplot2::geom_point(size = size, alpha = alpha) +
    ggplot2::labs(x = x_label, y = y_label) +
    ggsci::scale_color_d3() +
    ggplot2::theme_bw()

  # Step 6: Return results as a list following project conventions
  list(
    result_pcoa = pcoa_res,
    plot_pcoa = plot_pcoa,
    point_data = point_data,
    eigenvalue_pcoa = eigenvalue_pcoa
  )
}
