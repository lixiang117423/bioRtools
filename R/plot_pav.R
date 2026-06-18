#' Plot Pangenome Presence/Absence Matrix with Category Annotation
#'
#' Creates a combined plot showing a feature category annotation bar on top
#' and a presence/absence heatmap below, commonly used for pangenome visualization.
#'
#' @param data A data frame with at least three columns: one for group/category,
#'   one for sample/genome identifiers, and one for feature names.
#' @param group_col Character string specifying the column containing
#'   group/category labels (e.g., "Core", "Accessory", "Unique").
#'   Default is "group".
#' @param sample_col Character string specifying the column containing
#'   sample/genome identifiers. Default is "sample".
#' @param feature_col Character string specifying the column containing
#'   feature identifiers (e.g., gene/orthogroup names). Default is "feature".
#' @param present_color Character string, fill color for "Present" cells.
#'   Default is "#A63446".
#' @param absent_color Character string, fill color for "Absent" cells.
#'   Default is "#0C6291".
#' @param category_height Numeric, height ratio for the category annotation bar.
#'   Default is 0.05.
#' @param annotate_axes Logical, whether to show axis text and ticks.
#'   Default is FALSE (clean plot for large datasets).
#'
#' @return A combined plot object (from \code{aplot::insert_top}).
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' \dontrun{
#' # Pangenome visualization
#' df <- data.frame(
#'   Category = c("Core","Core","Core","Shell","Shell","Unique"),
#'   Genome = c("G1","G2","G3","G1","G2","G3"),
#'   Gene = c("OG1","OG1","OG1","OG2","OG2","OG3")
#' )
#' plot_pav(df, group_col = "Category", sample_col = "Genome", feature_col = "Gene")
#' }
plot_pav <- function(data, group_col = "group", sample_col = "sample",
                     feature_col = "feature",
                     present_color = "#A63446", absent_color = "#0C6291",
                     category_height = 0.05, annotate_axes = FALSE) {

  if (!requireNamespace("aplot", quietly = TRUE)) {
    stop("Package 'aplot' is required. Install with install.packages('aplot')")
  }

  # Validate inputs
  if (!is.data.frame(data)) stop("'data' must be a data frame")
  for (col in c(group_col, sample_col, feature_col)) {
    if (!col %in% names(data)) stop("Column '", col, "' not found in data")
  }

  # Prepare long format: distinct rows
  df <- data[, c(group_col, sample_col, feature_col)] %>%
    dplyr::distinct() %>%
    dplyr::mutate(value = "Present")

  # Order features by occurrence (most common first)
  df_x_order <- df %>%
    dplyr::select(!!rlang::sym(feature_col), !!rlang::sym(sample_col)) %>%
    dplyr::distinct() %>%
    dplyr::group_by(!!rlang::sym(feature_col)) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(.data$n) %>%
    dplyr::mutate(!!feature_col := factor(
      !!rlang::sym(feature_col), levels = unique(!!rlang::sym(feature_col))))

  # Category annotation bar
  axis_theme <- if (annotate_axes) {
    ggplot2::theme(
      legend.position = "none",
      axis.title = ggplot2::element_blank())
  } else {
    ggplot2::theme(
      legend.position = "none",
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank())
  }

  p_type <- df %>%
    dplyr::select(!!rlang::sym(group_col), !!rlang::sym(feature_col)) %>%
    dplyr::distinct() %>%
    ggplot2::ggplot(ggplot2::aes(
      x = !!rlang::sym(feature_col), y = 1,
      fill = !!rlang::sym(group_col))) +
    ggplot2::geom_tile() +
    scale_fill_research() +
    scale_x_discrete(limits = df_x_order[[feature_col]]) +
    axis_theme

  # PAV heatmap
  df_pav <- df %>%
    dplyr::select(!!rlang::sym(sample_col), !!rlang::sym(feature_col)) %>%
    dplyr::distinct() %>%
    dplyr::mutate(value = "Present") %>%
    tidyr::pivot_wider(
      names_from = !!rlang::sym(feature_col),
      values_from = .data$value) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ replace(.x, is.na(.x), "Absent"))) %>%
    tidyr::pivot_longer(
      cols = -!!rlang::sym(sample_col),
      names_to = feature_col,
      values_to = "value") %>%
    dplyr::mutate(value = factor(.data$value, levels = c("Present", "Absent")))

  p_pav <- df_pav %>%
    ggplot2::ggplot(ggplot2::aes(
      x = !!rlang::sym(feature_col),
      y = !!rlang::sym(sample_col),
      fill = .data$value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_x_discrete(limits = df_x_order[[feature_col]]) +
    ggplot2::scale_fill_manual(values = c(present_color, absent_color)) +
    axis_theme

  # Combine
  p <- aplot::insert_top(p_pav, p_type, height = category_height)
  p
}
