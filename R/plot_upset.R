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
#' @param relative_height Height ratio between upper and lower panels (default: 2).
#' @param relative_width Width ratio for the upset layout (default: 0.3).
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
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
plot_upset <- function(data,
                       group = NULL,
                       feature = NULL,
                       n_intersections = 30,
                       relative_height = 2,
                       relative_width = 0.3) {

  if (!requireNamespace("ggVennDiagram", quietly = TRUE)) {
    stop("Package 'ggVennDiagram' is required. Install with: install.packages('ggVennDiagram')")
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

  # -- PAV matrix -------------------------------------------------------------
  all_groups <- sort(unique(df_pa$group))
  df_pa$present <- 1L
  df_pav <- tidyr::pivot_wider(
    df_pa,
    names_from = group, values_from = present,
    values_fill = 0L
  )
  df_pav <- df_pav[, c("feature", all_groups), drop = FALSE]

  # -- Convert to named list for ggVennDiagram --------------------------------
  groups_list <- split(df_pa$feature[df_pa$present == 1L],
                       df_pa$group[df_pa$present == 1L])
  groups_list <- groups_list[all_groups]

  # -- Plot -------------------------------------------------------------------
  venn <- ggVennDiagram::Venn(groups_list)
  p <- ggVennDiagram::plot_upset(
    venn,
    nintersects = n_intersections,
    relative_height = relative_height,
    relative_width = relative_width
  )

  list(plot = p, data.pav = df_pav)
}
