#' Create a dual Y-axis plot with area and scatter layers
#'
#' @description
#' Create a dual Y-axis plot combining an area chart (mapped to the right Y-axis)
#' with scatter points (mapped to the left Y-axis). Supports optional phase
#' annotations with colored rectangles and labels. Useful for visualizing
#' time-series data with different measurement scales, such as rates alongside
#' percentages.
#'
#' @param data A data frame in wide format containing the plotting data.
#' @param x Character string specifying the column name for the x-axis.
#' @param area_col Character string specifying the column name for the area chart,
#'   mapped to the right Y-axis. Default is NULL.
#' @param scatter_cols Character vector of column names for scatter points,
#'   mapped to the left Y-axis. Default is NULL.
#' @param left_range Numeric vector of length 2 for left Y-axis range (c(min, max)).
#'   Auto-calculated from scatter data if NULL.
#' @param right_range Numeric vector of length 2 for right Y-axis range (c(min, max)).
#'   Auto-calculated from area data if NULL.
#' @param phases Optional data frame with columns: \code{start}, \code{end},
#'   \code{label}, and optionally \code{fill_color}, \code{border_color}.
#'   Adds colored rectangles at the top of the plot to indicate phases.
#' @param area_fill Character string for area fill color (supports 8-digit hex
#'   with alpha, e.g., "#bad6f980"). Default is "#bad6f980".
#' @param area_fill_name Character string for the area legend label.
#'   Default uses the area column name.
#' @param area_border_color Character string for area border color.
#'   Default is "grey70".
#' @param scatter_colors Named character vector of colors for scatter series.
#'   Default uses RColorBrewer "Set1" palette.
#' @param scatter_shapes Named integer vector of shapes for scatter series.
#'   Default uses sequential integers starting from 1.
#' @param left_label Character string for left Y-axis label. Default is "Value".
#' @param right_label Character string for right Y-axis label.
#'   Default uses the area column name.
#' @param x_label Character string for x-axis label. Default uses the x column name.
#' @param x_breaks Numeric vector of breaks for the x-axis. Default is NULL (auto).
#' @param left_breaks Numeric vector of breaks for the left Y-axis. Default is NULL (auto).
#' @param right_breaks Numeric vector of breaks for the right Y-axis. Default is NULL (auto).
#' @param point_size Numeric value for scatter point size. Default is 2.5.
#' @param legend_position Character string or numeric vector for legend position.
#'   Default is "top".
#' @param base_size Numeric value for base font size. Default is 11.
#'
#' @return A ggplot object.
#'
#' @details
#' The dual Y-axis is implemented by scaling the area data to fit the left Y-axis
#' range, then using \code{sec_axis()} to display the original scale on the right.
#' The scaling factor is calculated as \code{left_range[2] / right_range[2]}.
#'
#' For best results, provide explicit \code{left_range} and \code{right_range}
#' to ensure clean axis labels.
#'
#' The \code{phases} data frame, if provided, must contain columns \code{start},
#' \code{end}, and \code{label}. Optional columns \code{fill_color} and
#' \code{border_color} control the appearance of phase rectangles.
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' \dontrun{
#' library(bioRtools)
#'
#' df <- data.frame(
#'   time = 1:100,
#'   rate_a = runif(100, 1, 5),
#'   rate_b = runif(100, 2, 6),
#'   efficiency = runif(100, 20, 90)
#' )
#'
#' phases <- data.frame(
#'   start = c(0, 50),
#'   end = c(50, 100),
#'   label = c("Phase A", "Phase B"),
#'   fill_color = c("#90719f90", "#FDDDA090")
#' )
#'
#' p <- plot_dual_axis(
#'   data = df,
#'   x = "time",
#'   area_col = "efficiency",
#'   scatter_cols = c("rate_a", "rate_b"),
#'   left_range = c(0, 6),
#'   right_range = c(0, 100),
#'   phases = phases,
#'   left_label = "Rate (mmol/L/day)",
#'   right_label = "Efficiency (%)"
#' )
#' print(p)
#' }
plot_dual_axis <- function(data,
                           x,
                           area_col = NULL,
                           scatter_cols = NULL,
                           left_range = NULL,
                           right_range = NULL,
                           phases = NULL,
                           area_fill = "#bad6f980",
                           area_fill_name = NULL,
                           area_border_color = "grey70",
                           scatter_colors = NULL,
                           scatter_shapes = NULL,
                           left_label = NULL,
                           right_label = NULL,
                           x_label = NULL,
                           x_breaks = NULL,
                           left_breaks = NULL,
                           right_breaks = NULL,
                           point_size = 2.5,
                           legend_position = "top",
                           base_size = 11) {

  # Input validation
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame", call. = FALSE)
  }

  if (is.null(area_col) && is.null(scatter_cols)) {
    stop("At least one of 'area_col' or 'scatter_cols' must be provided", call. = FALSE)
  }

  required_cols <- x
  if (!is.null(area_col)) required_cols <- c(required_cols, area_col)
  if (!is.null(scatter_cols)) required_cols <- c(required_cols, scatter_cols)
  missing_cols <- required_cols[!required_cols %in% colnames(data)]
  if (length(missing_cols) > 0) {
    stop(
      "Required columns missing: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  # Set default labels
  if (is.null(area_fill_name) && !is.null(area_col)) {
    area_fill_name <- area_col
  }
  if (is.null(left_label)) {
    left_label <- "Value"
  }
  if (is.null(right_label) && !is.null(area_col)) {
    right_label <- area_col
  }
  if (is.null(x_label)) {
    x_label <- x
  }

  # Auto-calculate ranges
  if (is.null(left_range) && !is.null(scatter_cols)) {
    scatter_max <- max(data[, scatter_cols, drop = FALSE], na.rm = TRUE)
    left_range <- c(0, ceiling(scatter_max))
  }
  if (is.null(right_range) && !is.null(area_col)) {
    area_max <- max(data[[area_col]], na.rm = TRUE)
    right_range <- c(0, ceiling(area_max))
  }
  if (is.null(left_range) && !is.null(area_col)) {
    left_range <- right_range
  }

  # Scaling factor: maps right Y values to left Y range
  scale_factor <- left_range[2] / right_range[2]

  # Prepare area data with scaled values
  if (!is.null(area_col)) {
    area_data <- data %>%
      dplyr::select(dplyr::all_of(c(x, area_col))) %>%
      dplyr::arrange(!!rlang::sym(x)) %>%
      dplyr::mutate(
        .scaled = !!rlang::sym(area_col) * scale_factor,
        .group = area_fill_name
      )
  }

  # Prepare scatter data in long format
  if (!is.null(scatter_cols)) {
    scatter_data <- data %>%
      dplyr::select(dplyr::all_of(c(x, scatter_cols))) %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(scatter_cols),
        names_to = ".series",
        values_to = ".value"
      )
  }

  # Set scatter aesthetic defaults
  n_scatter <- length(scatter_cols)
  if (is.null(scatter_colors) && !is.null(scatter_cols)) {
    pal <- RColorBrewer::brewer.pal(max(3, n_scatter), "Set1")
    scatter_colors <- stats::setNames(pal[seq_len(n_scatter)], scatter_cols)
  }
  if (is.null(scatter_shapes) && !is.null(scatter_cols)) {
    scatter_shapes <- stats::setNames(seq_len(n_scatter), scatter_cols)
  }

  # Build plot
  p <- ggplot2::ggplot()

  # Area layer (border color outside aes avoids color scale conflict with scatter)
  if (!is.null(area_col)) {
    p <- p +
      ggplot2::geom_area(
        data = area_data,
        ggplot2::aes(
          x = !!rlang::sym(x),
          y = .data$.scaled,
          fill = .data$.group
        ),
        color = area_border_color,
        key_glyph = "rect"
      ) +
      ggplot2::scale_fill_manual(
        values = stats::setNames(area_fill, area_fill_name),
        guide = ggplot2::guide_legend(
          theme = ggplot2::theme(
            legend.key.height = ggplot2::unit(0.3, "cm"),
            legend.key.width = ggplot2::unit(0.6, "cm")
          )
        )
      )
  }

  # Scatter layer
  if (!is.null(scatter_cols)) {
    p <- p +
      ggplot2::geom_point(
        data = scatter_data,
        ggplot2::aes(
          x = !!rlang::sym(x),
          y = .data$.value,
          color = .data$.series,
          shape = .data$.series
        ),
        size = point_size
      ) +
      ggplot2::scale_color_manual(values = scatter_colors) +
      ggplot2::scale_shape_manual(values = scatter_shapes)
  }

  # Phase annotations
  if (!is.null(phases)) {
    phase_required <- c("start", "end", "label")
    phase_missing <- phase_required[!phase_required %in% colnames(phases)]
    if (length(phase_missing) > 0) {
      stop(
        "'phases' must have columns: ",
        paste(phase_required, collapse = ", "),
        call. = FALSE
      )
    }
    if (!"fill_color" %in% colnames(phases)) {
      phases$fill_color <- "#90719f90"
    }
    if (!"border_color" %in% colnames(phases)) {
      phases$border_color <- "black"
    }

    # Phase bar position: top 7% of plot area
    ymin_bar <- left_range[2] * 0.93
    ymax_bar <- left_range[2]
    y_text <- left_range[2] * 0.965

    for (i in seq_len(nrow(phases))) {
      p <- p +
        ggplot2::annotate(
          "rect",
          xmin = phases$start[i],
          xmax = phases$end[i],
          ymin = ymin_bar,
          ymax = ymax_bar,
          fill = phases$fill_color[i],
          color = phases$border_color[i]
        ) +
        ggplot2::annotate(
          "text",
          x = (phases$start[i] + phases$end[i]) / 2,
          y = y_text,
          label = phases$label[i],
          fontface = "bold"
        )
    }

    # Vertical lines at phase boundaries
    boundaries <- phases$end[-nrow(phases)]
    for (b in boundaries) {
      p <- p + ggplot2::geom_vline(xintercept = b)
    }
  }

  # Coordinate system
  p <- p + ggplot2::coord_cartesian(clip = "off")

  # X-axis scale
  x_scale_args <- list(expand = ggplot2::expansion(mult = c(0, 0)))
  if (!is.null(x_breaks)) x_scale_args$breaks <- x_breaks
  p <- p + do.call(ggplot2::scale_x_continuous, x_scale_args)

  # Y-axis scale with optional secondary axis
  y_args <- list(
    name = left_label,
    limits = left_range,
    expand = ggplot2::expansion(mult = c(0, 0))
  )
  if (!is.null(left_breaks)) y_args$breaks <- left_breaks

  if (!is.null(area_col)) {
    sec_args <- list(
      trans = ~ . / scale_factor,
      name = right_label
    )
    if (!is.null(right_breaks)) sec_args$breaks <- right_breaks
    y_args$sec.axis <- do.call(ggplot2::sec_axis, sec_args)
  }

  p <- p + do.call(ggplot2::scale_y_continuous, y_args)

  # Theme
  p <- p +
    ggplot2::theme_test() +
    ggplot2::theme(
      axis.title.y.right = ggplot2::element_text(
        angle = 90, color = "black", face = "bold", size = base_size
      ),
      axis.title.y.left = ggplot2::element_text(
        color = "black", face = "bold", size = base_size
      ),
      axis.title.x = ggplot2::element_text(
        color = "black", face = "bold", size = base_size
      ),
      axis.text = ggplot2::element_text(color = "black", size = base_size - 1),
      legend.title = ggplot2::element_blank(),
      legend.position = legend_position,
      legend.text = ggplot2::element_text(color = "black", face = "bold")
    )

  p
}
