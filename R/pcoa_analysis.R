#' Perform principal coordinate analysis (PCoA).
#'
#' @param data Feature table, with rows representing samples and columns representing feature values, such as OTUs.
#' @param sample Sample table, with the first column containing sample names. There are no specific requirements for the names of the subsequent columns, but the sample names must match those in the feature table.
#' @param method The method for calculating distances in vegan::vegdist() is defaulted to Bray-Curtis.
#' @param x The principal coordinate used for the X-axis in the plot is defaulted to PCo1.
#' @param y The principal coordinate used for the y-axis in the plot is defaulted to PCo2.
#' @param size The size of point.
#' @param color Column names in the sample data frame used for coloring the points.
#' @param alpha The alpha of point.
#'
#' @return A list containing four components:
#' \describe{
#'   \item{result.pcoa}{The output from ape::pcoa(), which can be called directly.}
#'   \item{plot.pcoa}{The plotting results, defaulting to PCo1 and PCo2. The output is a ggplot object, which can be fine-tuned using ggplot2.}
#'   \item{point.data}{The data used for plotting, which users can call for their own plots or export for use in other software.}
#'   \item{eigenvalue.pcoa}{The explained variance of the principal coordinates, starting from PCo1 by default.}
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
#' pcoa_analysis(data = df.pcoa.otu, sample = df.pcoa.sample) -> pcoa.res
#'
pcoa_analysis <- function(data, sample, method = "bray", x = "PCo1", y = "PCo2", size = 2, color = "group", alpha = 1) {
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
    ape::pcoa() -> pcoa.res

  # Step 2: Extract eigenvalues and calculate explained variance
  pcoa.res$values %>%
    as.data.frame() %>%
    dplyr::select(2) %>%
    dplyr::mutate(pcoa = paste0("PCo", seq_len(nrow(.)))) %>%
    dplyr::select(pcoa, Relative_eig) %>%
    dplyr::mutate(
      relative_eig = round(Relative_eig * 100, 2),
      percentage = paste0(relative_eig, "%"),
      label = paste0(pcoa, " (", percentage, ")")
    ) -> eigenvalue.pcoa

  # Step 3: Extract coordinate data and merge with sample information
  pcoa.res$vectors %>%
    as.data.frame() %>%
    dplyr::select(1:5) %>%
    {
      df <- .
      pco_labels <- paste0("PCo", seq_len(ncol(df)), " (",
        eigenvalue.pcoa$percentage[seq_len(ncol(df))], ")")
      magrittr::set_names(df, pco_labels)
    } %>%
    tibble::rownames_to_column(var = "sample") %>%
    dplyr::left_join(sample, by = "sample") -> point.data

  # Step 4: Create axis labels with explained variance
  x_label <- eigenvalue.pcoa %>%
    dplyr::filter(pcoa == x) %>%
    dplyr::pull(label)

  y_label <- eigenvalue.pcoa %>%
    dplyr::filter(pcoa == y) %>%
    dplyr::pull(label)

  # Step 5: Create PCoA plot
  plot.pcoa <- point.data %>%
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
    result.pcoa = pcoa.res,
    plot.pcoa = plot.pcoa,
    point.data = point.data,
    eigenvalue.pcoa = eigenvalue.pcoa
  )
}
