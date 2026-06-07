#' UpSet Plot from a Group-Feature Data Frame
#'
#' Creates an UpSet plot from a simple two-column data frame where one column
#' is the group label and the other is the feature identifier. Useful for
#' visualizing overlap of differential features across pairwise comparisons.
#'
#' Works directly with \code{\link{find_dams_deseq2}} output:
#' \preformatted{
#' res <- find_dams_deseq2(mat, meta, groupCol = "treatment")
#' plot_upset(res, group = "comparison", feature = "feature_id")
#' }
#'
#' @param data A data frame with at least two columns: one for group labels
#'   and one for feature identifiers.
#' @param group Column name for group labels (default: auto-detect
#'   "comparison" or "cluster").
#' @param feature Column name for feature identifiers (default: auto-detect
#'   "feature_id" or "gene").
#' @param n_intersections Maximum number of intersections to show (default: 30).
#' @param order_by How to order intersections: "freq" (default) or "degree".
#' @param fill Bar fill color (default: "#2874A6").
#' @param show_counts Show count labels on bars (default: TRUE).
#'
#' @return A list:
#'   \describe{
#'     \item{plot}{ggplot object.}
#'     \item{data.pav}{Presence/absence matrix (feature × group).}
#'   }
#'
#' @examples
#' # With find_dams_deseq2 output
#' \dontrun{
#' res <- find_dams_deseq2(mat, meta, groupCol = "treatment")
#' result <- plot_upset(res, group = "comparison", feature = "feature_id")
#' result$plot
#' head(result$data.pav)
#' }
#'
#' @seealso \code{\link{find_dams_deseq2}}, \code{\link{plot_tax_upset}}
#'
#' @export
plot_upset <- function(data,
                       group = NULL,
                       feature = NULL,
                       n_intersections = 30,
                       order_by = "freq",
                       fill = "#2874A6",
                       show_counts = TRUE) {

  if (!requireNamespace("ggupset", quietly = TRUE)) {
    stop("Package 'ggupset' is required. Install with: install.packages('ggupset')")
  }

  data <- as.data.frame(data)

  # -- Auto-detect columns ----------------------------------------------------
  nms <- names(data)
  if (is.null(group)) {
    group <- intersect(c("comparison", "cluster", "group"), nms)
    if (length(group) == 0) stop("Cannot auto-detect group column. Specify 'group'.")
    group <- group[1]
  }
  if (is.null(feature)) {
    feature <- intersect(c("feature_id", "gene", "feature"), nms)
    if (length(feature) == 0) stop("Cannot auto-detect feature column. Specify 'feature'.")
    feature <- feature[1]
  }

  if (!group %in% nms) stop(sprintf("Column '%s' not found", group))
  if (!feature %in% nms) stop(sprintf("Column '%s' not found", feature))

  # -- Build presence/absence -------------------------------------------------
  df_pa <- data.frame(
    group   = as.character(data[[group]]),
    feature = as.character(data[[feature]]),
    stringsAsFactors = FALSE
  )
  df_pa <- unique(df_pa)

  # -- UpSet summary ----------------------------------------------------------
  df_upset <- dplyr::summarize(
    dplyr::group_by(df_pa, feature),
    groups = list(sort(unique(group))),
    .groups = "drop"
  )

  # -- PAV matrix -------------------------------------------------------------
  all_groups <- sort(unique(df_pa$group))
  df_pa$present <- 1L
  df_pav <- tidyr::pivot_wider(
    df_pa,
    names_from = group, values_from = present,
    values_fill = 0L
  )
  df_pav <- df_pav[, c("feature", all_groups), drop = FALSE]

  # -- Plot -------------------------------------------------------------------
  p <- ggplot2::ggplot(df_upset, ggplot2::aes(x = groups)) +
    ggplot2::geom_bar(fill = fill, width = 0.6) +
    ggupset::scale_x_upset(
      n_intersections = n_intersections,
      order_by = order_by
    ) +
    ggplot2::labs(y = "Number of features")

  if (show_counts) {
    p <- p +
      ggplot2::geom_text(
        stat = "count",
        ggplot2::aes(label = ggplot2::after_stat(count)),
        vjust = -0.3, size = 3.5
      )
  }

  p <- p +
    ggupset::theme_combmatrix() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "grey92", linewidth = 0.3)
    )

  list(plot = p, data.pav = df_pav)
}
