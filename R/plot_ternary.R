#' Create a ternary plot for compositional data
#'
#' @description
#' Create a ternary (triangular) plot to visualize three-part compositional data.
#' Supports both discrete and continuous fill variables, automatically detected
#' based on the data type of the fill column. Useful for visualizing cell state
#' proportions, taxonomic composition, or any data where three components sum to
#' a constant.
#'
#' @param data A data frame containing the plotting data.
#' @param x Character string specifying the column name for the x-axis (top vertex).
#' @param y Character string specifying the column name for the y-axis (left vertex).
#' @param z Character string specifying the column name for the z-axis (right vertex).
#' @param fill Character string specifying the column name for fill color.
#'   Automatically uses discrete or continuous scale based on column type.
#'   Default is NULL (no fill mapping).
#' @param point_size Numeric value for point size. Default is 3.
#' @param point_shape Integer specifying the point shape. Default is 21
#'   (filled circle with border).
#' @param point_stroke Numeric value for border width of points. Default is 0.8.
#' @param border_color Character string for point border color. Default is "black".
#' @param title Character string for plot title. Default is NULL.
#' @param x_label Character string for x-axis label. Default is the value of \code{x}.
#' @param y_label Character string for y-axis label. Default is the value of \code{y}.
#' @param z_label Character string for z-axis label. Default is the value of \code{z}.
#' @param show_axis_labels Logical. Whether to show axis labels on the ternary plot.
#'   Default is FALSE.
#' @param fill_colors Named character vector of colors for discrete fill mode.
#'   Names should match the levels of the fill column. If NULL, a default palette
#'   is used. Ignored for continuous fill.
#' @param fill_palette Character string specifying the RColorBrewer palette for
#'   continuous fill mode. Default is "Spectral". Ignored for discrete fill.
#' @param fill_direction Integer (1 or -1) for the direction of the continuous
#'   color palette. Default is -1. Ignored for discrete fill.
#' @param axis_limits Numeric vector of length 2 for ternary axis limits.
#'   Default is NULL (use data range). Commonly set to \code{c(0, 1)}.
#' @param legend_position Numeric vector of length 2 specifying the position
#'   of the legend inside the plot area (0-1 coordinates). Default is
#'   \code{c(0.7, 0.8)}. Set to NULL to use default position.
#'
#' @return A ggtern/ggplot object that can be further customized or printed.
#'
#' @details
#' This function requires the \pkg{ggtern} package. Install it with
#' \code{install.packages("ggtern")}.
#'
#' The fill type (discrete vs continuous) is automatically determined by the
#' data type of the fill column: numeric columns use a continuous gradient scale,
#' while character or factor columns use a discrete categorical scale.
#'
#' For discrete fill, provide \code{fill_colors} as a named vector where names
#' match the unique values of the fill column. For continuous fill, the
#' \code{fill_palette} parameter accepts any RColorBrewer palette name.
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' \dontrun{
#' library(bioRtools)
#'
#' # Discrete fill by group
#' df <- data.frame(
#'   th2 = runif(50, 0.1, 0.6),
#'   th1 = runif(50, 0.1, 0.5),
#'   th17 = runif(50, 0.05, 0.3),
#'   group = sample(c("Young", "Older"), 50, replace = TRUE)
#' )
#'
#' p1 <- plot_ternary(
#'   data = df,
#'   x = "th2", y = "th1", z = "th17",
#'   fill = "group",
#'   fill_colors = c(Young = "#35978f", Older = "#bf812d"),
#'   title = "T cell states"
#' )
#' print(p1)
#'
#' # Continuous fill
#' df$age <- sample(25:80, 50, replace = TRUE)
#' p2 <- plot_ternary(
#'   data = df,
#'   x = "th2", y = "th1", z = "th17",
#'   fill = "age",
#'   fill_palette = "Spectral",
#'   axis_limits = c(0, 1)
#' )
#' print(p2)
#' }
plot_ternary <- function(data,
                         x,
                         y,
                         z,
                         fill = NULL,
                         point_size = 3,
                         point_shape = 21,
                         point_stroke = 0.8,
                         border_color = "black",
                         title = NULL,
                         x_label = x,
                         y_label = y,
                         z_label = z,
                         show_axis_labels = FALSE,
                         fill_colors = NULL,
                         fill_palette = "Spectral",
                         fill_direction = -1,
                         axis_limits = NULL,
                         legend_position = c(0.7, 0.8)) {

  # Check required package
  if (!requireNamespace("ggtern", quietly = TRUE)) {
    stop(
      "Package 'ggtern' is required for ternary plots. ",
      "Install it with: install.packages('ggtern')",
      call. = FALSE
    )
  }

  # Input validation
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame", call. = FALSE)
  }

  required_cols <- c(x, y, z)
  if (!is.null(fill)) required_cols <- c(required_cols, fill)
  missing_cols <- required_cols[!required_cols %in% colnames(data)]
  if (length(missing_cols) > 0) {
    stop(
      "Required columns missing from data: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.numeric(point_size) || length(point_size) != 1 || point_size <= 0) {
    stop("'point_size' must be a single positive numeric value", call. = FALSE)
  }

  # Build aesthetic mapping
  if (!is.null(fill)) {
    aes_mapping <- ggplot2::aes(
      x = !!rlang::sym(x),
      y = !!rlang::sym(y),
      z = !!rlang::sym(z),
      fill = !!rlang::sym(fill)
    )
  } else {
    aes_mapping <- ggplot2::aes(
      x = !!rlang::sym(x),
      y = !!rlang::sym(y),
      z = !!rlang::sym(z)
    )
  }

  # Create base ternary plot
  p <- ggtern::ggtern(data = data, mapping = aes_mapping) +
    ggplot2::geom_point(
      size = point_size,
      shape = point_shape,
      stroke = point_stroke,
      color = border_color
    ) +
    ggtern::theme_rgbw() +
    ggplot2::labs(
      x = x_label,
      y = y_label,
      z = z_label,
      title = title
    )

  # Add fill scale based on data type
  if (!is.null(fill)) {
    is_continuous <- is.numeric(data[[fill]])

    if (is_continuous) {
      p <- p +
        ggplot2::scale_fill_distiller(
          palette = fill_palette,
          direction = fill_direction
        ) +
        ggplot2::guides(
          fill = ggplot2::guide_colourbar()
        )

      if (!is.null(axis_limits)) {
        p <- p +
          ggtern::scale_T_continuous(limits = axis_limits) +
          ggtern::scale_L_continuous(limits = axis_limits) +
          ggtern::scale_R_continuous(limits = axis_limits)
      }
    } else {
      if (!is.null(fill_colors)) {
        p <- p + ggplot2::scale_fill_manual(values = fill_colors)
      } else {
        p <- p + ggplot2::scale_fill_brewer(palette = "Set2")
      }
      p <- p +
        ggplot2::guides(
          fill = ggplot2::guide_legend()
        )
    }
  }

  # Apply theme
  p <- p + ggplot2::theme(
    legend.background = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.key = ggplot2::element_blank()
  )

  if (!is.null(legend_position)) {
    p <- p + ggplot2::theme(
      legend.position.inside = legend_position
    )
  }

  if (!show_axis_labels) {
    p <- p + ggplot2::theme(
      tern.axis.title.T = ggplot2::element_blank(),
      tern.axis.title.L = ggplot2::element_blank(),
      tern.axis.title.R = ggplot2::element_blank()
    )
  }

  if (!is.null(title)) {
    p <- p + ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 12)
    )
  }

  p
}
