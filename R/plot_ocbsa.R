#' Genome-wide SNP-index Plot for OcBSA Results
#'
#' Replicates the OcBSA-style genome-wide SNP-index plot: chromosomes are laid
#' end to end on a cumulative x-axis, the binomial-distribution threshold curve
#' (\code{VirtualIndex}) and both pool index curves are drawn as lines, and
#' Δ(SNP-index) windows are shown as points colored by marker count.
#'
#' @param data A file path (character) or data frame from OcBSA. Expected
#'   columns: Chr, Pos, ChrLen, VirtualIndex, Pool1Index, Pool2Index,
#'   DeltaIndex, MarkerNum.
#' @param chr_sep Optional separator for cleaning chromosome names. When
#'   \code{"_"}, \code{"So120_Chr1"} becomes \code{"Chr1"} by keeping the
#'   \code{chr_index}-th segment. \code{NULL} (default) uses the Chr column
#'   as is.
#' @param chr_index 1-based segment to keep when \code{chr_sep} is set
#'   (default: 2).
#' @param line_colors Named character vector of colors for the three index
#'   lines; names must cover VirtualIndex, Pool1Index, Pool2Index.
#' @param point_size Size of the Δ(SNP-index) points (default: 1.5).
#' @param base_size Base font size passed to \code{theme_minimal()}
#'   (default: 13).
#'
#' @return A ggplot object.
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' \dontrun{
#' data(df.ocbsa.window, package = "bioRtools")
#' p <- plot_ocbsa(df.ocbsa.window, chr_sep = "_", chr_index = 2)
#' p
#' }
plot_ocbsa <- function(data,
                       chr_sep = NULL,
                       chr_index = 2,
                       line_colors = c(VirtualIndex = "red",
                                       Pool1Index = "green3",
                                       Pool2Index = "blue3"),
                       point_size = 1.5,
                       base_size = 13) {

  if (is.character(data)) {
    df <- read_data(data)
  } else if (is.data.frame(data)) {
    df <- data
  } else {
    stop("'data' must be a file path (character) or data frame")
  }

  required_cols <- c("Chr", "Pos", "ChrLen", "VirtualIndex",
                     "Pool1Index", "Pool2Index", "DeltaIndex", "MarkerNum")
  missing <- setdiff(required_cols, names(df))
  if (length(missing) > 0) {
    stop("Required columns missing: ", paste(missing, collapse = ", "))
  }

  numeric_cols <- setdiff(required_cols, "Chr")
  not_numeric <- numeric_cols[!vapply(df[numeric_cols], is.numeric, logical(1))]
  if (length(not_numeric) > 0) {
    stop("Columns must be numeric, got '", class(df[[not_numeric[1]]])[1],
         "' for: ", paste(not_numeric, collapse = ", "))
  }

  if (!is.null(chr_sep) &&
      (!is.character(chr_sep) || length(chr_sep) != 1L || !nzchar(chr_sep))) {
    stop("'chr_sep' must be NULL or a single non-empty string")
  }
  if (!is.numeric(chr_index) || length(chr_index) != 1L || is.na(chr_index) ||
      chr_index < 1 || chr_index != floor(chr_index)) {
    stop("'chr_index' must be a single positive integer")
  }

  index_cols <- c("VirtualIndex", "Pool1Index", "Pool2Index")
  missing_names <- setdiff(index_cols, names(line_colors))
  if (length(missing_names) > 0) {
    stop("'line_colors' must be named with: ",
         paste(index_cols, collapse = ", "))
  }

  if (!is.null(chr_sep)) {
    df$Chr <- split_pick(df$Chr, sep = chr_sep, index = chr_index)
  }

  # Natural sort keeps Chr2 before Chr10 regardless of prefix state
  chr_levels <- stringr::str_sort(unique(df$Chr), numeric = TRUE)

  chr_info <- df %>%
    dplyr::group_by(.data$Chr) %>%
    dplyr::summarise(chrlen = max(.data$ChrLen), .groups = "drop") %>%
    dplyr::mutate(chr_order = match(.data$Chr, chr_levels)) %>%
    dplyr::arrange(.data$chr_order) %>%
    dplyr::mutate(
      cum_start = c(0, utils::head(cumsum(.data$chrlen), -1)),
      mid = .data$cum_start + .data$chrlen / 2
    )

  # geom_line needs windows ordered along the cumulative axis
  plot_df <- df %>%
    dplyr::left_join(chr_info[, c("Chr", "cum_start")], by = "Chr") %>%
    dplyr::arrange(.data$cum_start, .data$Pos) %>%
    dplyr::mutate(x = .data$Pos + .data$cum_start)

  p <- ggplot2::ggplot(plot_df) +
    ggplot2::geom_vline(
      data = chr_info,
      ggplot2::aes(xintercept = .data$cum_start),
      linetype = "dashed", colour = "grey40", linewidth = 0.3
    ) +
    ggplot2::geom_line(
      ggplot2::aes(x = .data$x, y = .data$VirtualIndex),
      colour = line_colors[["VirtualIndex"]], linewidth = 0.5
    ) +
    ggplot2::geom_line(
      ggplot2::aes(x = .data$x, y = .data$Pool1Index),
      colour = line_colors[["Pool1Index"]], linewidth = 0.3
    ) +
    ggplot2::geom_line(
      ggplot2::aes(x = .data$x, y = .data$Pool2Index),
      colour = line_colors[["Pool2Index"]], linewidth = 0.3
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = .data$x, y = .data$DeltaIndex, colour = .data$MarkerNum),
      size = point_size
    ) +
    ggplot2::scale_colour_viridis_c(
      option = "plasma", direction = -1,
      name = "Marker\ncount"
    ) +
    ggplot2::scale_x_continuous(
      breaks = chr_info$mid, labels = chr_levels,
      expand = c(0.01, 0.01)
    ) +
    ggplot2::labs(x = NULL, y = "SNP-index") +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )

  p
}
