#' Prism-inspired ggplot2 Theme with Visible Legend Title
#'
#' A modified version of \code{ggprism::theme_prism()} that preserves the
#' Prism-style aesthetics but keeps the legend title visible instead of
#' hiding it.
#'
#' @param palette Character string specifying the Prism palette.
#'   See \code{names(ggprism::ggprism_data$themes)} for valid options.
#'   Default is "black_and_white".
#' @param base_size Numeric value specifying the base font size in points.
#'   Default is 14.
#' @param base_family Character string specifying the base font family.
#'   Default is "sans".
#' @param base_fontface Character string specifying the base font face.
#'   Default is "bold".
#' @param base_line_size Numeric value specifying the base line width.
#'   Default is \code{base_size / 14}.
#' @param base_rect_size Numeric value specifying the base rectangle line width.
#'   Default is \code{base_size / 14}.
#' @param axis_text_angle Numeric value specifying the axis text rotation angle.
#'   Must be one of 0, 45, 90, or 270. Default is 0.
#' @param border Logical indicating whether to draw a panel border.
#'   Default is FALSE.
#'
#' @return A ggplot2 theme object.
#' @export
#'
#' @examples
#' library(ggplot2)
#'
#' ggplot(iris, aes(x = Species, y = Sepal.Length, fill = Species)) +
#'   geom_boxplot() +
#'   bioRtools::theme_prism()
theme_prism <- function(palette = "black_and_white",
                        base_size = 14,
                        base_family = "sans",
                        base_fontface = "bold",
                        base_line_size = base_size / 14,
                        base_rect_size = base_size / 14,
                        axis_text_angle = 0,
                        border = FALSE) {
  angle <- axis_text_angle[1]

  if (!angle %in% c(0, 45, 90, 270)) {
    stop(sprintf("'axis_text_angle' must be one of [%s]",
      paste(c(0, 45, 90, 270), collapse = ", ")),
      ".\nFor other angles, use the guide_axis() function in ggplot2 instead",
      call. = FALSE
    )
  }

  if (!palette %in% names(ggprism::ggprism_data$themes)) {
    stop("The palette ", paste(palette),
      " does not exist.\n         See names(ggprism_data$themes) for valid palette names")
  }

  colours <- tibble::deframe(ggprism::ggprism_data$themes[[palette]])

  if (!rlang::is_bool(border)) {
    stop("border must be either: TRUE or FALSE")
  }

  if (border) {
    panel.border <- ggplot2::element_rect(fill = NA)
    axis.line <- ggplot2::element_blank()
  } else {
    panel.border <- ggplot2::element_blank()
    axis.line <- ggplot2::element_line()
  }

  t <- ggplot2::theme(
    line = ggplot2::element_line(
      colour = colours["axisColor"],
      linewidth = base_line_size, linetype = 1, lineend = "square"
    ),
    rect = ggplot2::element_rect(
      fill = "white", colour = colours["axisColor"],
      linewidth = base_rect_size, linetype = 1
    ),
    text = ggplot2::element_text(
      family = base_family, face = base_fontface,
      colour = colours["graphTitleColor"],
      size = base_size, lineheight = 0.9, hjust = 0.5,
      vjust = 0.5, angle = 0, margin = ggplot2::margin(), debug = FALSE
    ),
    prism.ticks.length = ggplot2::unit(base_size / 5, "pt"),
    axis.line = axis.line,
    axis.line.x = NULL,
    axis.line.y = NULL,
    axis.text = ggplot2::element_text(
      size = ggplot2::rel(0.95),
      colour = colours["axisLabelColor"]
    ),
    axis.text.x = ggplot2::element_text(
      margin = ggplot2::margin(t = 0.8 * base_size / 4),
      angle = axis_text_angle,
      hjust = ifelse(axis_text_angle %in% c(45, 90, 270), 1, 0.5),
      vjust = ifelse(axis_text_angle %in% c(0, 90, 270), 0.5, 1)
    ),
    axis.text.x.top = ggplot2::element_text(
      margin = ggplot2::margin(b = 0.8 * base_size / 4), vjust = 0
    ),
    axis.text.y = ggplot2::element_text(
      margin = ggplot2::margin(r = 0.5 * base_size / 4), hjust = 1
    ),
    axis.text.y.right = ggplot2::element_text(
      margin = ggplot2::margin(l = 0.5 * base_size / 4), hjust = 0
    ),
    axis.ticks = ggplot2::element_line(),
    axis.ticks.length = ggplot2::unit(base_size / 2.5, "pt"),
    axis.ticks.length.x = NULL,
    axis.ticks.length.x.top = NULL,
    axis.ticks.length.x.bottom = NULL,
    axis.ticks.length.y = NULL,
    axis.ticks.length.y.left = NULL,
    axis.ticks.length.y.right = NULL,
    axis.title = ggplot2::element_text(colour = colours["axisTitleColor"]),
    axis.title.x = ggplot2::element_text(
      margin = ggplot2::margin(t = base_size * 0.6), vjust = 1
    ),
    axis.title.x.top = ggplot2::element_text(
      margin = ggplot2::margin(b = base_size * 0.6), vjust = 0
    ),
    axis.title.y = ggplot2::element_text(
      angle = 90, margin = ggplot2::margin(r = base_size * 0.6), vjust = 1
    ),
    axis.title.y.right = ggplot2::element_text(
      angle = -90, margin = ggplot2::margin(l = base_size * 0.6), vjust = 0
    ),
    legend.background = ggplot2::element_blank(),
    legend.spacing = ggplot2::unit(base_size, "pt"),
    legend.spacing.x = NULL,
    legend.spacing.y = NULL,
    legend.margin = ggplot2::margin(
      base_size / 2, base_size / 2, base_size / 2, base_size / 2
    ),
    legend.key = ggplot2::element_blank(),
    legend.key.size = ggplot2::unit(1.2, "lines"),
    legend.key.height = NULL,
    legend.key.width = ggplot2::unit(base_size * 1.8, "pt"),
    legend.text = ggplot2::element_text(size = ggplot2::rel(0.8), face = "plain"),
    legend.text.align = NULL,
    legend.title = ggplot2::element_text(
      size = ggplot2::rel(0.9), face = "plain",
      colour = colours["axisTitleColor"]
    ),
    legend.title.align = NULL,
    legend.position = "right",
    legend.direction = NULL,
    legend.justification = "center",
    legend.box = NULL,
    legend.box.margin = ggplot2::margin(0, 0, 0, 0, "cm"),
    legend.box.background = ggplot2::element_blank(),
    legend.box.spacing = ggplot2::unit(base_size, "pt"),
    panel.background = ggplot2::element_rect(
      fill = ifelse(palette == "office", colours["plottingAreaColor"], NA),
      colour = NA
    ),
    panel.border = panel.border,
    panel.grid = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.spacing = ggplot2::unit(base_size / 2, "pt"),
    panel.spacing.x = NULL,
    panel.spacing.y = NULL,
    panel.ontop = FALSE,
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(
      colour = colours["axisTitleColor"],
      size = ggplot2::rel(0.8),
      margin = ggplot2::margin(
        base_size / 2.5, base_size / 2.5,
        base_size / 2.5, base_size / 2.5
      )
    ),
    strip.text.x = ggplot2::element_text(
      margin = ggplot2::margin(b = base_size / 3)
    ),
    strip.text.y = ggplot2::element_text(
      angle = -90, margin = ggplot2::margin(l = base_size / 3)
    ),
    strip.text.y.left = ggplot2::element_text(angle = 90),
    strip.placement = "inside",
    strip.placement.x = NULL,
    strip.placement.y = NULL,
    strip.switch.pad.grid = ggplot2::unit(base_size / 4, "pt"),
    strip.switch.pad.wrap = ggplot2::unit(base_size / 4, "pt"),
    plot.background = ggplot2::element_rect(
      fill = colours["pageBackgroundColor"], colour = NA
    ),
    plot.title = ggplot2::element_text(
      size = ggplot2::rel(1.2), hjust = 0.5, vjust = 1,
      margin = ggplot2::margin(b = base_size)
    ),
    plot.title.position = "panel",
    plot.subtitle = ggplot2::element_text(
      hjust = 0.5, vjust = 1,
      margin = ggplot2::margin(b = base_size / 2)
    ),
    plot.caption = ggplot2::element_text(
      size = ggplot2::rel(0.8), hjust = 1, vjust = 1,
      margin = ggplot2::margin(t = base_size / 2)
    ),
    plot.caption.position = "panel",
    plot.tag = ggplot2::element_text(
      size = ggplot2::rel(1.2), hjust = 0.5, vjust = 0.5
    ),
    plot.tag.position = "topleft",
    plot.margin = ggplot2::margin(
      base_size / 2, base_size / 2, base_size / 2, base_size / 2
    ),
    complete = TRUE
  )

  parent <- ggplot2::theme(!!!ggprism::ggprism_data$themes[["all_null"]])

  if (!"legend.text.align" %in% rlang::fn_fmls_names(ggplot2::theme)) {
    t$legend.text.align <- parent$legend.text.align <- NULL
    t$legend.title.align <- parent$legend.title.align <- NULL
  }

  parent %+replace% t
}
