#' Generate iTOL Configuration Files
#'
#' A unified wrapper around `itol.toolkit` to generate iTOL annotation files
#' for tree visualization. Supports color strips (TREE_COLORS), symbols
#' (DATASET_SYMBOL), bar charts (DATASET_SIMPLEBAR), and heatmaps
#' (DATASET_HEATMAP).
#'
#' @param data Data frame where the first column contains node/tip IDs matching
#'   the tree, and remaining columns contain the values to visualize
#' @param tree A phylo object (from `ape::read.tree()` or similar)
#' @param key Character string label for the iTOL annotation unit
#' @param type Type of annotation to create. One of:
#'   \describe{
#'     \item{`"color"`}{Color strips for categorical variables (TREE_COLORS/range)}
#'     \item{`"symbol"`}{Symbols at tree nodes (DATASET_SYMBOL)}
#'     \item{`"bar"`}{Bar chart of numeric values (DATASET_SIMPLEBAR)}
#'     \item{`"heatmap"`}{Heatmap of numeric matrix (DATASET_HEATMAP)}
#'   }
#' @param output File path to write the iTOL annotation file. If NULL, the unit
#'   object is returned without writing to disk
#' @param color_min Minimum color for heatmap gradient (hex code, e.g. "#ffffff")
#' @param color_max Maximum color for heatmap gradient (hex code, e.g. "#8ccdd7")
#' @param ... Additional parameters passed to `itol.toolkit::create_unit()`
#'   (e.g., `subtype`, `position`)
#'
#' @return The itol.toolkit unit object (invisibly if `output` is specified)
#' @export
#'
#' @examples
#' \dontrun{
#' library(ape)
#' tree <- read.tree("tree.nwk")
#'
#' # Color strips for categorical groups
#' df_group <- data.frame(id = c("sp1", "sp2"), Phylum = c("Ascomycota", "Basidiomycota"))
#' itol_config(df_group, tree, key = "Phylum", type = "color", output = "phylum.txt")
#'
#' # Bar chart
#' df_bar <- data.frame(id = c("sp1", "sp2"), count = c(5, 12))
#' itol_config(df_bar, tree, key = "Gene count", type = "bar", output = "count.txt")
#'
#' # Heatmap with custom colors
#' df_heat <- data.frame(id = c("sp1", "sp2"), var1 = c(1, 2), var2 = c(3, 1))
#' itol_config(df_heat, tree, key = "Expression", type = "heatmap",
#'   color_min = "#ffffff", color_max = "#8ccdd7", output = "heatmap.txt")
#' }
itol_config <- function(data, tree, key, type = c("color", "symbol", "bar", "heatmap"),
                        output = NULL, color_min = NULL, color_max = NULL, ...) {
  type <- match.arg(type)

  if (!requireNamespace("itol.toolkit", quietly = TRUE)) {
    stop("Package 'itol.toolkit' is required. Install it from CRAN.")
  }

  type_map <- c(
    color = "TREE_COLORS",
    symbol = "DATASET_SYMBOL",
    bar = "DATASET_SIMPLEBAR",
    heatmap = "DATASET_HEATMAP"
  )

  itol_type <- unname(type_map[type])
  dots <- list(...)

  # Build create_unit arguments
  args <- list(
    data = data,
    key = key,
    type = itol_type,
    tree = tree
  )

  # Add type-specific defaults
  if (type == "color" && is.null(dots$subtype)) {
    args$subtype <- "range"
  }

  # Merge extra arguments
  args <- c(args, dots)

  unit <- do.call(itol.toolkit::create_unit, args)

  # Post-creation customization for heatmap colors
  if (type == "heatmap") {
    if (!is.null(color_min)) unit@specific_themes$heatmap$color$min <- color_min
    if (!is.null(color_max)) unit@specific_themes$heatmap$color$max <- color_max
  }

  if (!is.null(output)) {
    itol.toolkit::write_unit(unit, output)
  }

  invisible(unit)
}
