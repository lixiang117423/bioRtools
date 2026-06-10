#' Plot Venn/Euler Diagram
#'
#' Creates Venn or Euler diagrams from a long-format data frame with group and
#' feature columns. Supports two modes: \code{"euler"} (proportional area,
#' default) and \code{"classic"} (equal-size circles).
#'
#' @param data A data frame with at least two columns: one for group membership
#'   and one for feature identifiers.
#' @param group_col Character string specifying the column containing group names.
#'   Default is "group".
#' @param feature_col Character string specifying the column containing feature
#'   identifiers. Default is "feature".
#' @param type Character string: \code{"euler"} (default) for proportional-area
#'   Euler diagrams via \pkg{eulerr}, or \code{"classic"} for equal-size Venn
#'   diagrams via \pkg{ggvenn}.
#' @param show_sets Optional character vector specifying which sets to display.
#'   If NULL (default), all sets are used. Useful when there are more than 4 groups.
#' @param fill_colors Optional character vector of fill colors. If NULL, a
#'   default palette is used.
#' @param fill_alpha Numeric, fill transparency (0-1). Default is 0.58.
#'   Only used when \code{type = "euler"}.
#' @param stroke_color Character string, border color. Default is "#3f3f46".
#' @param stroke_width Numeric, border line width. Default is 1.2.
#' @param label_size Numeric, set label font size. Default is 12.
#' @param label_color Character string, set label color. Default is "#202124".
#' @param show_counts Logical, whether to display counts. Default is TRUE.
#' @param show_percent Logical, whether to display percentages alongside counts.
#'   Default is FALSE. Only used when \code{type = "euler"}.
#' @param count_size Numeric, count label font size. Default is 10.
#' @param count_color Character string, count label color. Default is "#202124".
#' @param shape Character string, shape for Euler fitting. Default is "ellipse".
#'   Only used when \code{type = "euler"}.
#'
#' @return A plot object (ggplot for classic mode, grid object for euler mode).
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' \dontrun{
#' # Example data
#' df <- data.frame(
#'   feature = c("g1", "g1", "g2", "g3", "g3", "g4", "g5"),
#'   group = c("A", "B", "A", "B", "C", "C", "A")
#' )
#'
#' # Euler (proportional area)
#' plot_venn(df, feature_col = "feature", group_col = "group")
#'
#' # Classic (equal size)
#' plot_venn(df, type = "classic")
#'
#' # Custom colors
#' plot_venn(df, fill_colors = c("#4E79D9", "#E36C61", "#58A55C"))
#'
#' # Select specific sets
#' plot_venn(df, show_sets = c("A", "B"))
#' }
plot_venn <- function(data, group_col = "group", feature_col = "feature",
                      type = c("euler", "classic"),
                      show_sets = NULL,
                      fill_colors = NULL, fill_alpha = 0.58,
                      stroke_color = "#3f3f46", stroke_width = 1.2,
                      label_size = 12, label_color = "#202124",
                      show_counts = TRUE, show_percent = FALSE,
                      count_size = 10, count_color = "#202124",
                      shape = "ellipse") {

  type <- match.arg(type)

  # Validate inputs
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  if (!group_col %in% names(data)) {
    stop("'group_col' (", group_col, ") not found in data")
  }
  if (!feature_col %in% names(data)) {
    stop("'feature_col' (", feature_col, ") not found in data")
  }

  # Convert to named list: group -> vector of features
  feature_list <- split(as.character(data[[feature_col]]), data[[group_col]])
  feature_list <- lapply(feature_list, unique)

  # Filter to selected sets
  if (!is.null(show_sets)) {
    feature_list <- feature_list[show_sets]
  }

  n_sets <- length(feature_list)
  if (n_sets < 2) {
    stop("Need at least 2 sets for a Venn diagram")
  }
  if (type == "classic" && n_sets > 4) {
    stop("Classic mode (ggvenn) supports at most 4 sets. Use type = 'euler' or specify show_sets.")
  }

  # Default colors
  if (is.null(fill_colors)) {
    default_palettes <- list(
      c("#4E79D9", "#E36C61"),
      c("#4E79D9", "#E36C61", "#58A55C"),
      c("#4E79D9", "#E36C61", "#58A55C", "#F28E2B"),
      c("#4E79D9", "#E36C61", "#58A55C", "#F28E2B", "#B07AA1")
    )
    fill_colors <- default_palettes[[min(n_sets, 5)]]
  }

  if (type == "classic") {
    # Classic Venn using ggvenn
    if (!requireNamespace("ggvenn", quietly = TRUE)) {
      stop("Package 'ggvenn' is required. Install with install.packages('ggvenn')")
    }

    p <- ggvenn::ggvenn(
      feature_list,
      stroke_linetype = 1,
      stroke_size = stroke_width,
      stroke_color = stroke_color,
      set_name_color = label_color,
      set_name_size = label_size,
      fill_color = fill_colors,
      text_size = count_size,
      text_color = count_color,
      show_percentage = show_percent,
      show_counts = show_counts
    )

    return(p)

  } else {
    # Euler diagram using eulerr
    if (!requireNamespace("eulerr", quietly = TRUE)) {
      stop("Package 'eulerr' is required. Install with install.packages('eulerr')")
    }

    # Compute intersection counts
    all_combos <- list()
    set_names <- names(feature_list)

    for (k in seq_len(n_sets)) {
      combos <- combn(set_names, k, simplify = FALSE)
      for (combo in combos) {
        if (k == n_sets) {
          # Full intersection
          features_in_all <- Reduce(intersect, feature_list[combo])
          others <- setdiff(set_names, combo)
          # Features only in these sets (not in any other)
          if (length(others) > 0) {
            features_in_others <- Reduce(union, feature_list[others])
            combo_only <- setdiff(features_in_all, features_in_others)
          } else {
            combo_only <- features_in_all
          }
        } else {
          # Features in all these sets but not in any of the others
          features_in_all <- Reduce(intersect, feature_list[combo])
          others <- setdiff(set_names, combo)
          features_in_others <- Reduce(union, feature_list[others])
          combo_only <- setdiff(features_in_all, features_in_others)
        }
        combo_name <- paste(combo, collapse = "&")
        all_combos[[combo_name]] <- length(combo_only)
      }
    }

    set_counts <- unlist(all_combos)
    fit <- eulerr::euler(set_counts, shape = shape)

    quantity_type <- if (show_percent && show_counts) c("counts", "percent") else
                    if (show_percent) "percent" else "counts"

    # eulerr fills need colors for each set, not each subset
    plot_colors <- rep_len(fill_colors, n_sets)

    p <- plot(
      fit,
      fills = list(fill = plot_colors, alpha = fill_alpha),
      edges = list(col = stroke_color, lwd = stroke_width),
      labels = list(fontsize = label_size, fontface = "bold", col = label_color),
      quantities = list(
        type = quantity_type,
        fontsize = count_size,
        col = count_color
      )
    )

    return(p)
  }
}
