# #' Professional ggplot2 theme for biological data visualization with enhanced parameters
# #'
# #' @description
# #' A clean and professional ggplot2 theme specifically designed for biological 
# #' data visualization. This enhanced version provides explicit parameters for 
# #' common theme modifications, enabling autocomplete functionality while 
# #' maintaining all original features.
# #'
# #' @param base_size Numeric value specifying the base font size in points. 
# #'   All other text elements are scaled relative to this size. Default is 10.
# #' @param base_family Character string specifying the base font family. 
# #'   Default is "" (system default).
# #' @param strip_text_italic Logical indicating whether strip text should be 
# #'   displayed in italic font when using facets. Default is FALSE.
# #' @param show_minor_ticks_x Logical indicating whether to display minor ticks 
# #'   on the X-axis. Default is FALSE.
# #' @param show_minor_ticks_y Logical indicating whether to display minor ticks 
# #'   on the Y-axis. Default is TRUE.
# #' @param show_legend Logical indicating whether to display the legend. 
# #'   Default is TRUE.
# #' @param panel_background_color Character string specifying the panel 
# #'   background color. Default is "white".
# #' @param panel_border_color Character string specifying the panel border 
# #'   color. Default is "black".
# #' @param panel_border_size Numeric value specifying the panel border line 
# #'   width. Default is 0.8.
# #' @param strip_background_color Character string specifying the strip 
# #'   background color for facets. Default is "#f0f0f0".
# #' @param grid_lines Logical indicating whether to show major grid lines. 
# #'   Default is FALSE for a clean look.
# #' @param axis_line Logical indicating whether to show axis lines. Default 
# #'   is FALSE (uses panel border instead).
# #'
# #' @param legend.position Character string or numeric vector specifying legend 
# #'   position. Options: "none", "left", "right", "bottom", "top", or numeric 
# #'   vector c(x, y). Default is NULL (uses show_legend parameter).
# #' @param plot.title.color Character string specifying plot title color. 
# #'   Default is NULL (uses theme default).
# #' @param plot.title.size Numeric value specifying plot title size as multiple 
# #'   of base_size. Default is NULL (uses 1.3).
# #' @param plot.title.hjust Numeric value (0-1) specifying plot title horizontal 
# #'   justification. Default is NULL (uses 0.5).
# #' @param plot.title.face Character string specifying plot title font face. 
# #'   Options: "plain", "bold", "italic", "bold.italic". Default is NULL (uses "bold").
# #' @param plot.subtitle.size Numeric value specifying subtitle size as multiple 
# #'   of base_size. Default is NULL (uses 1.1).
# #' @param plot.caption.size Numeric value specifying caption size as multiple 
# #'   of base_size. Default is NULL (uses 0.8).
# #'
# #' @param axis.text.x.angle Numeric value specifying X-axis text rotation angle 
# #'   in degrees. Default is NULL (no rotation).
# #' @param axis.text.x.hjust Numeric value (0-1) specifying X-axis text horizontal 
# #'   justification. Default is NULL (auto-adjusted based on angle).
# #' @param axis.text.x.vjust Numeric value (0-1) specifying X-axis text vertical 
# #'   justification. Default is NULL (uses 0.5).
# #' @param axis.text.x.size Numeric value specifying X-axis text size as multiple 
# #'   of base_size. Default is NULL (uses 1.0).
# #' @param axis.text.y.angle Numeric value specifying Y-axis text rotation angle 
# #'   in degrees. Default is NULL (no rotation).
# #' @param axis.text.y.size Numeric value specifying Y-axis text size as multiple 
# #'   of base_size. Default is NULL (uses 1.0).
# #'
# #' @param axis.title.x.face Character string specifying X-axis title font face. 
# #'   Options: "plain", "bold", "italic", "bold.italic". Default is NULL.
# #' @param axis.title.x.size Numeric value specifying X-axis title size as 
# #'   multiple of base_size. Default is NULL (uses 1.1).
# #' @param axis.title.x.color Character string specifying X-axis title color. 
# #'   Default is NULL.
# #' @param axis.title.y.face Character string specifying Y-axis title font face. 
# #'   Options: "plain", "bold", "italic", "bold.italic". Default is NULL.
# #' @param axis.title.y.size Numeric value specifying Y-axis title size as 
# #'   multiple of base_size. Default is NULL (uses 1.1).
# #' @param axis.title.y.color Character string specifying Y-axis title color. 
# #'   Default is NULL.
# #'
# #' @param panel.grid.major.x Grid lines for X-axis. Use element_line() or 
# #'   element_blank(). Default is NULL.
# #' @param panel.grid.major.y Grid lines for Y-axis. Use element_line() or 
# #'   element_blank(). Default is NULL.
# #' @param panel.grid.minor.x Minor grid lines for X-axis. Use element_line() or 
# #'   element_blank(). Default is NULL.
# #' @param panel.grid.minor.y Minor grid lines for Y-axis. Use element_line() or 
# #'   element_blank(). Default is NULL.
# #'
# #' @param strip.text.face Character string specifying strip text font face. 
# #'   Options: "plain", "bold", "italic", "bold.italic". Default is NULL.
# #' @param strip.text.size Numeric value specifying strip text size as multiple 
# #'   of base_size. Default is NULL (uses 1.0).
# #' @param strip.text.color Character string specifying strip text color. 
# #'   Default is NULL (uses "black").
# #' @param strip.background.fill Character string specifying strip background 
# #'   fill color. Default is NULL (uses strip_background_color parameter).
# #' @param strip.background.color Character string specifying strip background 
# #'   border color. Default is NULL.
# #'
# #' @param legend.text.size Numeric value specifying legend text size as multiple 
# #'   of base_size. Default is NULL (uses 1.0).
# #' @param legend.title.size Numeric value specifying legend title size as 
# #'   multiple of base_size. Default is NULL (uses 1.0).
# #' @param legend.title.face Character string specifying legend title font face. 
# #'   Options: "plain", "bold", "italic", "bold.italic". Default is NULL.
# #'
# #' @param ... Additional theme arguments passed to \code{ggplot2::theme()}. 
# #'   These won't have autocomplete but provide full flexibility.
# #'
# #' @return A ggplot2 theme object that can be added to any ggplot.
# #'
# #' @export
# #'
# #' @examples
# #' library(dplyr)
# #' library(ggplot2)
# #' library(bioRtools)
# #'
# #' # Basic usage with enhanced parameters (autocomplete supported)
# #' iris %>%
# #'   ggplot(aes(Species, Sepal.Length, fill = Species)) +
# #'   geom_boxplot() +
# #'   theme_bio(
# #'     base_size = 12,
# #'     legend.position = "none",
# #'     plot.title.color = "navy",
# #'     axis.text.x.angle = 45,
# #'     axis.title.x.face = "bold"
# #'   )
# #'
# theme_bio <- function(base_size = 10,
#                       base_family = "",
#                       strip_text_italic = FALSE,
#                       show_minor_ticks_x = FALSE,
#                       show_minor_ticks_y = TRUE,
#                       show_legend = TRUE,
#                       panel_background_color = "white",
#                       panel_border_color = "black", 
#                       panel_border_size = 0.8,
#                       strip_background_color = "#f0f0f0",
#                       grid_lines = FALSE,
#                       axis_line = FALSE,
                      
#                       # Legend parameters
#                       legend.position = NULL,
#                       legend.text.size = NULL,
#                       legend.title.size = NULL,
#                       legend.title.face = NULL,
                      
#                       # Plot title and labels
#                       plot.title.color = NULL,
#                       plot.title.size = NULL,
#                       plot.title.hjust = NULL,
#                       plot.title.face = NULL,
#                       plot.subtitle.size = NULL,
#                       plot.caption.size = NULL,
                      
#                       # Axis text parameters
#                       axis.text.x.angle = NULL,
#                       axis.text.x.hjust = NULL,
#                       axis.text.x.vjust = NULL,
#                       axis.text.x.size = NULL,
#                       axis.text.y.angle = NULL,
#                       axis.text.y.size = NULL,
                      
#                       # Axis title parameters
#                       axis.title.x.face = NULL,
#                       axis.title.x.size = NULL,
#                       axis.title.x.color = NULL,
#                       axis.title.y.face = NULL,
#                       axis.title.y.size = NULL,
#                       axis.title.y.color = NULL,
                      
#                       # Panel grid parameters
#                       panel.grid.major.x = NULL,
#                       panel.grid.major.y = NULL,
#                       panel.grid.minor.x = NULL,
#                       panel.grid.minor.y = NULL,
                      
#                       # Strip parameters
#                       strip.text.face = NULL,
#                       strip.text.size = NULL,
#                       strip.text.color = NULL,
#                       strip.background.fill = NULL,
#                       strip.background.color = NULL,
                      
#                       # Additional parameters
#                       ...) {
  
#   # Input validation
#   if (!is.numeric(base_size) || base_size <= 0) {
#     stop("'base_size' must be a positive number")
#   }
  
#   if (!is.character(base_family)) {
#     stop("'base_family' must be a character string")
#   }
  
#   if (!is.character(panel_background_color)) {
#     stop("'panel_background_color' must be a character string")
#   }
  
#   if (!is.character(panel_border_color)) {
#     stop("'panel_border_color' must be a character string")
#   }
  
#   if (!is.numeric(panel_border_size) || panel_border_size < 0) {
#     stop("'panel_border_size' must be a non-negative number")
#   }
  
#   # Build base theme
#   bio_theme <- ggthemes::theme_foundation(
#     base_size = base_size,
#     base_family = base_family
#   ) +
#     ggplot2::theme(
#       # Line and text aesthetics
#       line = ggplot2::element_line(colour = "black"),
#       rect = ggplot2::element_rect(colour = NA, linetype = 1),
#       text = ggplot2::element_text(colour = "black"),
      
#       # Axis elements
#       axis.line.x = if (axis_line) ggplot2::element_line() else ggplot2::element_blank(),
#       axis.line.y = if (axis_line) ggplot2::element_line() else ggplot2::element_blank(),
#       axis.text = ggplot2::element_text(size = ggplot2::rel(1)),
#       axis.text.x = ggplot2::element_text(
#         vjust = 0, 
#         margin = ggplot2::margin(t = base_size * 0.1, unit = "pt")
#       ),
#       axis.text.x.top = ggplot2::element_text(
#         vjust = 0, 
#         margin = ggplot2::margin(b = base_size, unit = "pt")
#       ),
#       axis.text.y = ggplot2::element_text(
#         hjust = 0.5, 
#         margin = ggplot2::margin(r = base_size * 0.1, unit = "pt")
#       ),
#       axis.ticks = ggplot2::element_line(),
#       axis.title = ggplot2::element_text(size = ggplot2::rel(1.1)),
#       axis.title.x = ggplot2::element_text(
#         margin = ggplot2::margin(t = base_size * 1, unit = "pt")
#       ),
#       axis.title.y = ggplot2::element_text(
#         angle = 90, 
#         margin = ggplot2::margin(r = base_size * 1, unit = "pt")
#       ),
#       axis.ticks.length = ggplot2::unit(base_size * 0.5, "points"),
      
#       # Legend elements
#       legend.background = ggplot2::element_rect(fill = "transparent", linetype = 0),
#       legend.spacing = ggplot2::unit(base_size * 1.5, "points"),
#       legend.key = ggplot2::element_rect(fill = "transparent", linetype = 0),
#       legend.key.size = ggplot2::unit(1.2, "lines"),
#       legend.key.height = NULL,
#       legend.key.width = NULL,
#       legend.text = ggplot2::element_text(size = ggplot2::rel(1.0)),
#       legend.text.align = NULL,
#       legend.title = ggplot2::element_text(size = ggplot2::rel(1.0), hjust = 0),
#       legend.title.align = NULL,
#       legend.direction = NULL,
#       legend.justification = "center",
      
#       # Panel elements
#       panel.background = ggplot2::element_rect(
#         fill = panel_background_color, 
#         linetype = 0
#       ),
#       panel.border = ggplot2::element_rect(
#         color = panel_border_color, 
#         linewidth = panel_border_size, 
#         linetype = "solid",
#         fill = NA
#       ),
#       panel.spacing = ggplot2::unit(0.5, "lines"),
#       panel.grid.major = if (grid_lines) ggplot2::element_line(color = "grey90", linewidth = 0.5) else ggplot2::element_blank(),
#       panel.grid.minor = ggplot2::element_blank(),
      
#       # Strip elements (for facets)
#       strip.background = ggplot2::element_rect(
#         colour = panel_border_color, 
#         fill = strip_background_color, 
#         linewidth = ggplot2::rel(1.0)
#       ),
#       strip.text.y = ggplot2::element_text(
#         color = "black", 
#         size = ggplot2::rel(1.0), 
#         angle = 0
#       ),
      
#       # Plot elements
#       plot.title = ggplot2::element_text(
#         size = ggplot2::rel(1.3), 
#         hjust = 0.5, 
#         vjust = 1.5,
#         face = "bold"
#       ),
#       plot.subtitle = ggplot2::element_text(
#         size = ggplot2::rel(1.1),
#         hjust = 0.5
#       ),
#       plot.caption = ggplot2::element_text(
#         size = ggplot2::rel(0.8),
#         hjust = 1
#       ),
#       plot.margin = ggplot2::unit(c(6, 5, 6, 5) * 2, "points"),
#       complete = FALSE
#     )
  
#   # Add strip text styling based on italic preference
#   if (strip_text_italic) {
#     bio_theme <- bio_theme + ggplot2::theme(
#       strip.text.x = ggplot2::element_text(
#         color = "black",
#         size = ggplot2::rel(1.0),
#         angle = 0,
#         face = "italic"
#       )
#     )
#   } else {
#     bio_theme <- bio_theme + ggplot2::theme(
#       strip.text.x = ggplot2::element_text(
#         color = "black",
#         size = ggplot2::rel(1.0),
#         angle = 0,
#         face = "plain"
#       )
#     )
#   }
  
#   # Add minor ticks if requested
#   if (show_minor_ticks_x) {
#     bio_theme <- bio_theme + ggplot2::theme(
#       axis.minor.ticks.length.x = ggplot2::unit(base_size * 0.25, "points")
#     )
#   }
  
#   if (show_minor_ticks_y) {
#     bio_theme <- bio_theme + ggplot2::theme(
#       axis.minor.ticks.length.y = ggplot2::unit(base_size * 0.25, "points")
#     )
#   }
  
#   # Handle legend display (backward compatibility)
#   if (!show_legend && is.null(legend.position)) {
#     legend.position <- "none"
#   }
  
#   # Build list of additional theme modifications
#   additional_theme_args <- list()
  
#   # Legend parameters
#   if (!is.null(legend.position)) {
#     additional_theme_args$legend.position <- legend.position
#   }
  
#   if (!is.null(legend.text.size)) {
#     additional_theme_args$legend.text <- ggplot2::element_text(size = ggplot2::rel(legend.text.size))
#   }
  
#   if (!is.null(legend.title.size) || !is.null(legend.title.face)) {
#     additional_theme_args$legend.title <- ggplot2::element_text(
#       size = if (!is.null(legend.title.size)) ggplot2::rel(legend.title.size) else ggplot2::rel(1.0),
#       face = if (!is.null(legend.title.face)) legend.title.face else "plain",
#       hjust = 0
#     )
#   }
  
#   # Plot title parameters
#   title_args <- list()
#   if (!is.null(plot.title.color)) title_args$colour <- plot.title.color
#   if (!is.null(plot.title.size)) title_args$size <- ggplot2::rel(plot.title.size)
#   if (!is.null(plot.title.hjust)) title_args$hjust <- plot.title.hjust
#   if (!is.null(plot.title.face)) title_args$face <- plot.title.face
  
#   if (length(title_args) > 0) {
#     # Keep defaults for unspecified parameters
#     if (is.null(title_args$size)) title_args$size <- ggplot2::rel(1.3)
#     if (is.null(title_args$hjust)) title_args$hjust <- 0.5
#     if (is.null(title_args$face)) title_args$face <- "bold"
#     title_args$vjust <- 1.5
    
#     additional_theme_args$plot.title <- do.call(ggplot2::element_text, title_args)
#   }
  
#   if (!is.null(plot.subtitle.size)) {
#     additional_theme_args$plot.subtitle <- ggplot2::element_text(
#       size = ggplot2::rel(plot.subtitle.size),
#       hjust = 0.5
#     )
#   }
  
#   if (!is.null(plot.caption.size)) {
#     additional_theme_args$plot.caption <- ggplot2::element_text(
#       size = ggplot2::rel(plot.caption.size),
#       hjust = 1
#     )
#   }
  
#   # Axis text parameters
#   if (!is.null(axis.text.x.angle) || !is.null(axis.text.x.hjust) || 
#       !is.null(axis.text.x.vjust) || !is.null(axis.text.x.size)) {
    
#     x_text_args <- list()
#     if (!is.null(axis.text.x.angle)) x_text_args$angle <- axis.text.x.angle
#     if (!is.null(axis.text.x.size)) x_text_args$size <- ggplot2::rel(axis.text.x.size)
    
#     # Auto-adjust hjust based on angle if not specified
#     if (!is.null(axis.text.x.hjust)) {
#       x_text_args$hjust <- axis.text.x.hjust
#     } else if (!is.null(axis.text.x.angle) && axis.text.x.angle != 0) {
#       x_text_args$hjust <- 1
#     }
    
#     if (!is.null(axis.text.x.vjust)) {
#       x_text_args$vjust <- axis.text.x.vjust
#     } else {
#       x_text_args$vjust <- 0
#     }
    
#     x_text_args$margin <- ggplot2::margin(t = base_size * 0.1, unit = "pt")
#     additional_theme_args$axis.text.x <- do.call(ggplot2::element_text, x_text_args)
#   }
  
#   if (!is.null(axis.text.y.angle) || !is.null(axis.text.y.size)) {
#     y_text_args <- list()
#     if (!is.null(axis.text.y.angle)) y_text_args$angle <- axis.text.y.angle
#     if (!is.null(axis.text.y.size)) y_text_args$size <- ggplot2::rel(axis.text.y.size)
    
#     y_text_args$hjust <- 0.5
#     y_text_args$margin <- ggplot2::margin(r = base_size * 0.1, unit = "pt")
#     additional_theme_args$axis.text.y <- do.call(ggplot2::element_text, y_text_args)
#   }
  
#   # Axis title parameters
#   if (!is.null(axis.title.x.face) || !is.null(axis.title.x.size) || !is.null(axis.title.x.color)) {
#     x_title_args <- list()
#     if (!is.null(axis.title.x.face)) x_title_args$face <- axis.title.x.face
#     if (!is.null(axis.title.x.size)) x_title_args$size <- ggplot2::rel(axis.title.x.size)
#     if (!is.null(axis.title.x.color)) x_title_args$colour <- axis.title.x.color
    
#     x_title_args$margin <- ggplot2::margin(t = base_size * 1, unit = "pt")
#     additional_theme_args$axis.title.x <- do.call(ggplot2::element_text, x_title_args)
#   }
  
#   if (!is.null(axis.title.y.face) || !is.null(axis.title.y.size) || !is.null(axis.title.y.color)) {
#     y_title_args <- list()
#     if (!is.null(axis.title.y.face)) y_title_args$face <- axis.title.y.face
#     if (!is.null(axis.title.y.size)) y_title_args$size <- ggplot2::rel(axis.title.y.size)
#     if (!is.null(axis.title.y.color)) y_title_args$colour <- axis.title.y.color
    
#     y_title_args$angle <- 90
#     y_title_args$margin <- ggplot2::margin(r = base_size * 1, unit = "pt")
#     additional_theme_args$axis.title.y <- do.call(ggplot2::element_text, y_title_args)
#   }
  
#   # Panel grid parameters
#   if (!is.null(panel.grid.major.x)) {
#     additional_theme_args$panel.grid.major.x <- panel.grid.major.x
#   }
  
#   if (!is.null(panel.grid.major.y)) {
#     additional_theme_args$panel.grid.major.y <- panel.grid.major.y
#   }
  
#   if (!is.null(panel.grid.minor.x)) {
#     additional_theme_args$panel.grid.minor.x <- panel.grid.minor.x
#   }
  
#   if (!is.null(panel.grid.minor.y)) {
#     additional_theme_args$panel.grid.minor.y <- panel.grid.minor.y
#   }
  
#   # Strip parameters
#   if (!is.null(strip.text.face) || !is.null(strip.text.size) || !is.null(strip.text.color)) {
#     strip_text_args <- list()
#     if (!is.null(strip.text.color)) strip_text_args$color <- strip.text.color
#     if (!is.null(strip.text.size)) strip_text_args$size <- ggplot2::rel(strip.text.size)
#     if (!is.null(strip.text.face)) strip_text_args$face <- strip.text.face
    
#     strip_text_args$angle <- 0
#     additional_theme_args$strip.text.x <- do.call(ggplot2::element_text, strip_text_args)
#   }
  
#   if (!is.null(strip.background.fill) || !is.null(strip.background.color)) {
#     strip_bg_args <- list()
#     if (!is.null(strip.background.fill)) strip_bg_args$fill <- strip.background.fill
#     if (!is.null(strip.background.color)) strip_bg_args$colour <- strip.background.color
    
#     strip_bg_args$linewidth <- ggplot2::rel(1.0)
#     additional_theme_args$strip.background <- do.call(ggplot2::element_rect, strip_bg_args)
#   }
  
#   # Apply additional theme arguments
#   if (length(additional_theme_args) > 0) {
#     bio_theme <- bio_theme + do.call(ggplot2::theme, additional_theme_args)
#   }
  
#   # Handle ... parameters (for full flexibility)
#   dots_args <- list(...)
#   if (length(dots_args) > 0) {
#     bio_theme <- bio_theme + do.call(ggplot2::theme, dots_args)
#   }
  
#   return(bio_theme)
# }

#' Professional ggplot2 theme for biological data visualization with enhanced parameters
#'
#' @description
#' A clean and professional ggplot2 theme specifically designed for biological 
#' data visualization. This enhanced version provides explicit parameters for 
#' common theme modifications, enabling autocomplete functionality while 
#' maintaining all original features. Now includes scale control options.
#'
#' @param base_size Numeric value specifying the base font size in points. 
#'   All other text elements are scaled relative to this size. Default is 10.
#' @param base_family Character string specifying the base font family. 
#'   Default is "" (system default).
#' @param strip_text_italic Logical indicating whether strip text should be 
#'   displayed in italic font when using facets. Default is FALSE.
#' @param show_minor_ticks_x Logical indicating whether to display minor ticks 
#'   on the X-axis. Default is FALSE.
#' @param show_minor_ticks_y Logical indicating whether to display minor ticks 
#'   on the Y-axis. Default is TRUE.
#' @param show_legend Logical indicating whether to display the legend. 
#'   Default is TRUE.
#' @param panel_background_color Character string specifying the panel 
#'   background color. Default is "white".
#' @param panel_border_color Character string specifying the panel border 
#'   color. Default is "black".
#' @param panel_border_size Numeric value specifying the panel border line 
#'   width. Default is 0.8.
#' @param strip_background_color Character string specifying the strip 
#'   background color for facets. Default is "#f0f0f0".
#' @param grid_lines Logical indicating whether to show major grid lines. 
#'   Default is FALSE for a clean look.
#' @param axis_line Logical indicating whether to show axis lines. Default 
#'   is FALSE (uses panel border instead).
#'
#' @param y_expand Y-axis expansion setting. Can be a numeric vector of length 2 
#'   (e.g., c(0, 0)) or NULL to use ggplot default. Default is c(0, 0) for no expansion.
#' @param x_expand X-axis expansion setting. Can be a numeric vector of length 2 
#'   or NULL to use ggplot default. Default is NULL.
#' @param return_list Logical indicating whether to return a list containing both 
#'   theme and scale objects, or just the theme. Default is FALSE (theme only).
#'
#' @param legend.position Character string or numeric vector specifying legend 
#'   position. Options: "none", "left", "right", "bottom", "top", or numeric 
#'   vector c(x, y). Default is NULL (uses show_legend parameter).
#' @param plot.title.color Character string specifying plot title color. 
#'   Default is NULL (uses theme default).
#' @param plot.title.size Numeric value specifying plot title size as multiple 
#'   of base_size. Default is NULL (uses 1.3).
#' @param plot.title.hjust Numeric value (0-1) specifying plot title horizontal 
#'   justification. Default is NULL (uses 0.5).
#' @param plot.title.face Character string specifying plot title font face. 
#'   Options: "plain", "bold", "italic", "bold.italic". Default is NULL (uses "bold").
#' @param plot.subtitle.size Numeric value specifying subtitle size as multiple 
#'   of base_size. Default is NULL (uses 1.1).
#' @param plot.caption.size Numeric value specifying caption size as multiple 
#'   of base_size. Default is NULL (uses 0.8).
#'
#' @param axis.text.x.angle Numeric value specifying X-axis text rotation angle 
#'   in degrees. Default is NULL (no rotation).
#' @param axis.text.x.hjust Numeric value (0-1) specifying X-axis text horizontal 
#'   justification. Default is NULL (auto-adjusted based on angle).
#' @param axis.text.x.vjust Numeric value (0-1) specifying X-axis text vertical 
#'   justification. Default is NULL (uses 0.5).
#' @param axis.text.x.size Numeric value specifying X-axis text size as multiple 
#'   of base_size. Default is NULL (uses 1.0).
#' @param axis.text.y.angle Numeric value specifying Y-axis text rotation angle 
#'   in degrees. Default is NULL (no rotation).
#' @param axis.text.y.size Numeric value specifying Y-axis text size as multiple 
#'   of base_size. Default is NULL (uses 1.0).
#'
#' @param axis.title.x.face Character string specifying X-axis title font face. 
#'   Options: "plain", "bold", "italic", "bold.italic". Default is NULL.
#' @param axis.title.x.size Numeric value specifying X-axis title size as 
#'   multiple of base_size. Default is NULL (uses 1.1).
#' @param axis.title.x.color Character string specifying X-axis title color. 
#'   Default is NULL.
#' @param axis.title.y.face Character string specifying Y-axis title font face. 
#'   Options: "plain", "bold", "italic", "bold.italic". Default is NULL.
#' @param axis.title.y.size Numeric value specifying Y-axis title size as 
#'   multiple of base_size. Default is NULL (uses 1.1).
#' @param axis.title.y.color Character string specifying Y-axis title color. 
#'   Default is NULL.
#'
#' @param panel.grid.major.x Grid lines for X-axis. Use element_line() or 
#'   element_blank(). Default is NULL.
#' @param panel.grid.major.y Grid lines for Y-axis. Use element_line() or 
#'   element_blank(). Default is NULL.
#' @param panel.grid.minor.x Minor grid lines for X-axis. Use element_line() or 
#'   element_blank(). Default is NULL.
#' @param panel.grid.minor.y Minor grid lines for Y-axis. Use element_line() or 
#'   element_blank(). Default is NULL.
#'
#' @param strip.text.face Character string specifying strip text font face. 
#'   Options: "plain", "bold", "italic", "bold.italic". Default is NULL.
#' @param strip.text.size Numeric value specifying strip text size as multiple 
#'   of base_size. Default is NULL (uses 1.0).
#' @param strip.text.color Character string specifying strip text color. 
#'   Default is NULL (uses "black").
#' @param strip.background.fill Character string specifying strip background 
#'   fill color. Default is NULL (uses strip_background_color parameter).
#' @param strip.background.color Character string specifying strip background 
#'   border color. Default is NULL.
#'
#' @param legend.text.size Numeric value specifying legend text size as multiple 
#'   of base_size. Default is NULL (uses 1.0).
#' @param legend.title.size Numeric value specifying legend title size as 
#'   multiple of base_size. Default is NULL (uses 1.0).
#' @param legend.title.face Character string specifying legend title font face. 
#'   Options: "plain", "bold", "italic", "bold.italic". Default is NULL.
#'
#' @param ... Additional theme arguments passed to \code{ggplot2::theme()}. 
#'   These won't have autocomplete but provide full flexibility.
#'
#' @return If return_list is FALSE (default): A ggplot2 theme object.
#'   If return_list is TRUE: A list containing 'theme' and scale objects 
#'   that can be added to any ggplot.
#'
#' @export
#'
#' @examples
#' library(dplyr)
#' library(ggplot2)
#' library(bioRtools)
#'
#' # Basic usage (theme only)
#' iris %>%
#'   ggplot(aes(Species, Sepal.Length, fill = Species)) +
#'   geom_boxplot() +
#'   theme_bio(base_size = 12, legend.position = "none") +
#'   scale_y_continuous(expand = c(0, 0))
#'
#' # Enhanced usage with automatic scale application
#' iris %>%
#'   ggplot(aes(Species, Sepal.Length, fill = Species)) +
#'   geom_boxplot() +
#'   theme_bio(
#'     base_size = 12,
#'     legend.position = "none",
#'     y_expand = c(0, 0),
#'     return_list = TRUE
#'   )
#'
theme_bio <- function(base_size = 10,
                      base_family = "",
                      strip_text_italic = FALSE,
                      show_minor_ticks_x = FALSE,
                      show_minor_ticks_y = TRUE,
                      show_legend = TRUE,
                      panel_background_color = "white",
                      panel_border_color = "black", 
                      panel_border_size = 0.8,
                      strip_background_color = "#f0f0f0",
                      grid_lines = FALSE,
                      axis_line = FALSE,
                      
                      # NEW: Scale control parameters
                      y_expand = c(0, 0),
                      x_expand = NULL,
                      return_list = FALSE,
                      
                      # Legend parameters
                      legend.position = NULL,
                      legend.text.size = NULL,
                      legend.title.size = NULL,
                      legend.title.face = NULL,
                      
                      # Plot title and labels
                      plot.title.color = NULL,
                      plot.title.size = NULL,
                      plot.title.hjust = NULL,
                      plot.title.face = NULL,
                      plot.subtitle.size = NULL,
                      plot.caption.size = NULL,
                      
                      # Axis text parameters
                      axis.text.x.angle = NULL,
                      axis.text.x.hjust = NULL,
                      axis.text.x.vjust = NULL,
                      axis.text.x.size = NULL,
                      axis.text.y.angle = NULL,
                      axis.text.y.size = NULL,
                      
                      # Axis title parameters
                      axis.title.x.face = NULL,
                      axis.title.x.size = NULL,
                      axis.title.x.color = NULL,
                      axis.title.y.face = NULL,
                      axis.title.y.size = NULL,
                      axis.title.y.color = NULL,
                      
                      # Panel grid parameters
                      panel.grid.major.x = NULL,
                      panel.grid.major.y = NULL,
                      panel.grid.minor.x = NULL,
                      panel.grid.minor.y = NULL,
                      
                      # Strip parameters
                      strip.text.face = NULL,
                      strip.text.size = NULL,
                      strip.text.color = NULL,
                      strip.background.fill = NULL,
                      strip.background.color = NULL,
                      
                      # Additional parameters
                      ...) {
  
  # Input validation
  if (!is.numeric(base_size) || base_size <= 0) {
    stop("'base_size' must be a positive number")
  }
  
  if (!is.character(base_family)) {
    stop("'base_family' must be a character string")
  }
  
  if (!is.character(panel_background_color)) {
    stop("'panel_background_color' must be a character string")
  }
  
  if (!is.character(panel_border_color)) {
    stop("'panel_border_color' must be a character string")
  }
  
  if (!is.numeric(panel_border_size) || panel_border_size < 0) {
    stop("'panel_border_size' must be a non-negative number")
  }
  
  # Validate new parameters
  if (!is.null(y_expand) && (!is.numeric(y_expand) || length(y_expand) != 2)) {
    stop("'y_expand' must be a numeric vector of length 2 or NULL")
  }
  
  if (!is.null(x_expand) && (!is.numeric(x_expand) || length(x_expand) != 2)) {
    stop("'x_expand' must be a numeric vector of length 2 or NULL")
  }
  
  if (!is.logical(return_list)) {
    stop("'return_list' must be logical (TRUE or FALSE)")
  }
  
  # Build base theme
  bio_theme <- ggthemes::theme_foundation(
    base_size = base_size,
    base_family = base_family
  ) +
    ggplot2::theme(
      # Line and text aesthetics
      line = ggplot2::element_line(colour = "black"),
      rect = ggplot2::element_rect(colour = NA, linetype = 1),
      text = ggplot2::element_text(colour = "black"),
      
      # Axis elements
      axis.line.x = if (axis_line) ggplot2::element_line() else ggplot2::element_blank(),
      axis.line.y = if (axis_line) ggplot2::element_line() else ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = ggplot2::rel(1)),
      axis.text.x = ggplot2::element_text(
        vjust = 0, 
        margin = ggplot2::margin(t = base_size * 0.1, unit = "pt")
      ),
      axis.text.x.top = ggplot2::element_text(
        vjust = 0, 
        margin = ggplot2::margin(b = base_size, unit = "pt")
      ),
      axis.text.y = ggplot2::element_text(
        hjust = 0.5, 
        margin = ggplot2::margin(r = base_size * 0.1, unit = "pt")
      ),
      axis.ticks = ggplot2::element_line(),
      axis.title = ggplot2::element_text(size = ggplot2::rel(1.1)),
      axis.title.x = ggplot2::element_text(
        margin = ggplot2::margin(t = base_size * 1, unit = "pt")
      ),
      axis.title.y = ggplot2::element_text(
        angle = 90, 
        margin = ggplot2::margin(r = base_size * 1, unit = "pt")
      ),
      axis.ticks.length = ggplot2::unit(base_size * 0.5, "points"),
      
      # Legend elements
      legend.background = ggplot2::element_rect(fill = "transparent", linetype = 0),
      legend.spacing = ggplot2::unit(base_size * 1.5, "points"),
      legend.key = ggplot2::element_rect(fill = "transparent", linetype = 0),
      legend.key.size = ggplot2::unit(1.2, "lines"),
      legend.key.height = NULL,
      legend.key.width = NULL,
      legend.text = ggplot2::element_text(size = ggplot2::rel(1.0)),
      legend.text.align = NULL,
      legend.title = ggplot2::element_text(size = ggplot2::rel(1.0), hjust = 0),
      legend.title.align = NULL,
      legend.direction = NULL,
      legend.justification = "center",
      
      # Panel elements
      panel.background = ggplot2::element_rect(
        fill = panel_background_color, 
        linetype = 0
      ),
      panel.border = ggplot2::element_rect(
        color = panel_border_color, 
        linewidth = panel_border_size, 
        linetype = "solid",
        fill = NA
      ),
      panel.spacing = ggplot2::unit(0.5, "lines"),
      panel.grid.major = if (grid_lines) ggplot2::element_line(color = "grey90", linewidth = 0.5) else ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      
      # Strip elements (for facets)
      strip.background = ggplot2::element_rect(
        colour = panel_border_color, 
        fill = strip_background_color, 
        linewidth = ggplot2::rel(1.0)
      ),
      strip.text.y = ggplot2::element_text(
        color = "black", 
        size = ggplot2::rel(1.0), 
        angle = 0
      ),
      
      # Plot elements
      plot.title = ggplot2::element_text(
        size = ggplot2::rel(1.3), 
        hjust = 0.5, 
        vjust = 1.5,
        face = "bold"
      ),
      plot.subtitle = ggplot2::element_text(
        size = ggplot2::rel(1.1),
        hjust = 0.5
      ),
      plot.caption = ggplot2::element_text(
        size = ggplot2::rel(0.8),
        hjust = 1
      ),
      plot.margin = ggplot2::unit(c(6, 5, 6, 5) * 2, "points"),
      complete = FALSE
    )
  
  # Add strip text styling based on italic preference
  if (strip_text_italic) {
    bio_theme <- bio_theme + ggplot2::theme(
      strip.text.x = ggplot2::element_text(
        color = "black",
        size = ggplot2::rel(1.0),
        angle = 0,
        face = "italic"
      )
    )
  } else {
    bio_theme <- bio_theme + ggplot2::theme(
      strip.text.x = ggplot2::element_text(
        color = "black",
        size = ggplot2::rel(1.0),
        angle = 0,
        face = "plain"
      )
    )
  }
  
  # Add minor ticks if requested
  if (show_minor_ticks_x) {
    bio_theme <- bio_theme + ggplot2::theme(
      axis.minor.ticks.length.x = ggplot2::unit(base_size * 0.25, "points")
    )
  }
  
  if (show_minor_ticks_y) {
    bio_theme <- bio_theme + ggplot2::theme(
      axis.minor.ticks.length.y = ggplot2::unit(base_size * 0.25, "points")
    )
  }
  
  # Handle legend display (backward compatibility)
  if (!show_legend && is.null(legend.position)) {
    legend.position <- "none"
  }
  
  # Build list of additional theme modifications
  additional_theme_args <- list()
  
  # Legend parameters
  if (!is.null(legend.position)) {
    additional_theme_args$legend.position <- legend.position
  }
  
  if (!is.null(legend.text.size)) {
    additional_theme_args$legend.text <- ggplot2::element_text(size = ggplot2::rel(legend.text.size))
  }
  
  if (!is.null(legend.title.size) || !is.null(legend.title.face)) {
    additional_theme_args$legend.title <- ggplot2::element_text(
      size = if (!is.null(legend.title.size)) ggplot2::rel(legend.title.size) else ggplot2::rel(1.0),
      face = if (!is.null(legend.title.face)) legend.title.face else "plain",
      hjust = 0
    )
  }
  
  # Plot title parameters
  title_args <- list()
  if (!is.null(plot.title.color)) title_args$colour <- plot.title.color
  if (!is.null(plot.title.size)) title_args$size <- ggplot2::rel(plot.title.size)
  if (!is.null(plot.title.hjust)) title_args$hjust <- plot.title.hjust
  if (!is.null(plot.title.face)) title_args$face <- plot.title.face
  
  if (length(title_args) > 0) {
    # Keep defaults for unspecified parameters
    if (is.null(title_args$size)) title_args$size <- ggplot2::rel(1.3)
    if (is.null(title_args$hjust)) title_args$hjust <- 0.5
    if (is.null(title_args$face)) title_args$face <- "bold"
    title_args$vjust <- 1.5
    
    additional_theme_args$plot.title <- do.call(ggplot2::element_text, title_args)
  }
  
  if (!is.null(plot.subtitle.size)) {
    additional_theme_args$plot.subtitle <- ggplot2::element_text(
      size = ggplot2::rel(plot.subtitle.size),
      hjust = 0.5
    )
  }
  
  if (!is.null(plot.caption.size)) {
    additional_theme_args$plot.caption <- ggplot2::element_text(
      size = ggplot2::rel(plot.caption.size),
      hjust = 1
    )
  }
  
  # Axis text parameters
  if (!is.null(axis.text.x.angle) || !is.null(axis.text.x.hjust) || 
      !is.null(axis.text.x.vjust) || !is.null(axis.text.x.size)) {
    
    x_text_args <- list()
    if (!is.null(axis.text.x.angle)) x_text_args$angle <- axis.text.x.angle
    if (!is.null(axis.text.x.size)) x_text_args$size <- ggplot2::rel(axis.text.x.size)
    
    # Auto-adjust hjust based on angle if not specified
    if (!is.null(axis.text.x.hjust)) {
      x_text_args$hjust <- axis.text.x.hjust
    } else if (!is.null(axis.text.x.angle) && axis.text.x.angle != 0) {
      x_text_args$hjust <- 1
    }
    
    if (!is.null(axis.text.x.vjust)) {
      x_text_args$vjust <- axis.text.x.vjust
    } else {
      x_text_args$vjust <- 0
    }
    
    x_text_args$margin <- ggplot2::margin(t = base_size * 0.1, unit = "pt")
    additional_theme_args$axis.text.x <- do.call(ggplot2::element_text, x_text_args)
  }
  
  if (!is.null(axis.text.y.angle) || !is.null(axis.text.y.size)) {
    y_text_args <- list()
    if (!is.null(axis.text.y.angle)) y_text_args$angle <- axis.text.y.angle
    if (!is.null(axis.text.y.size)) y_text_args$size <- ggplot2::rel(axis.text.y.size)
    
    y_text_args$hjust <- 0.5
    y_text_args$margin <- ggplot2::margin(r = base_size * 0.1, unit = "pt")
    additional_theme_args$axis.text.y <- do.call(ggplot2::element_text, y_text_args)
  }
  
  # Axis title parameters
  if (!is.null(axis.title.x.face) || !is.null(axis.title.x.size) || !is.null(axis.title.x.color)) {
    x_title_args <- list()
    if (!is.null(axis.title.x.face)) x_title_args$face <- axis.title.x.face
    if (!is.null(axis.title.x.size)) x_title_args$size <- ggplot2::rel(axis.title.x.size)
    if (!is.null(axis.title.x.color)) x_title_args$colour <- axis.title.x.color
    
    x_title_args$margin <- ggplot2::margin(t = base_size * 1, unit = "pt")
    additional_theme_args$axis.title.x <- do.call(ggplot2::element_text, x_title_args)
  }
  
  if (!is.null(axis.title.y.face) || !is.null(axis.title.y.size) || !is.null(axis.title.y.color)) {
    y_title_args <- list()
    if (!is.null(axis.title.y.face)) y_title_args$face <- axis.title.y.face
    if (!is.null(axis.title.y.size)) y_title_args$size <- ggplot2::rel(axis.title.y.size)
    if (!is.null(axis.title.y.color)) y_title_args$colour <- axis.title.y.color
    
    y_title_args$angle <- 90
    y_title_args$margin <- ggplot2::margin(r = base_size * 1, unit = "pt")
    additional_theme_args$axis.title.y <- do.call(ggplot2::element_text, y_title_args)
  }
  
  # Panel grid parameters
  if (!is.null(panel.grid.major.x)) {
    additional_theme_args$panel.grid.major.x <- panel.grid.major.x
  }
  
  if (!is.null(panel.grid.major.y)) {
    additional_theme_args$panel.grid.major.y <- panel.grid.major.y
  }
  
  if (!is.null(panel.grid.minor.x)) {
    additional_theme_args$panel.grid.minor.x <- panel.grid.minor.x
  }
  
  if (!is.null(panel.grid.minor.y)) {
    additional_theme_args$panel.grid.minor.y <- panel.grid.minor.y
  }
  
  # Strip parameters
  if (!is.null(strip.text.face) || !is.null(strip.text.size) || !is.null(strip.text.color)) {
    strip_text_args <- list()
    if (!is.null(strip.text.color)) strip_text_args$color <- strip.text.color
    if (!is.null(strip.text.size)) strip_text_args$size <- ggplot2::rel(strip.text.size)
    if (!is.null(strip.text.face)) strip_text_args$face <- strip.text.face
    
    strip_text_args$angle <- 0
    additional_theme_args$strip.text.x <- do.call(ggplot2::element_text, strip_text_args)
  }
  
  if (!is.null(strip.background.fill) || !is.null(strip.background.color)) {
    strip_bg_args <- list()
    if (!is.null(strip.background.fill)) strip_bg_args$fill <- strip.background.fill
    if (!is.null(strip.background.color)) strip_bg_args$colour <- strip.background.color
    
    strip_bg_args$linewidth <- ggplot2::rel(1.0)
    additional_theme_args$strip.background <- do.call(ggplot2::element_rect, strip_bg_args)
  }
  
  # Apply additional theme arguments
  if (length(additional_theme_args) > 0) {
    bio_theme <- bio_theme + do.call(ggplot2::theme, additional_theme_args)
  }
  
  # Handle ... parameters (for full flexibility)
  dots_args <- list(...)
  if (length(dots_args) > 0) {
    bio_theme <- bio_theme + do.call(ggplot2::theme, dots_args)
  }
  
  # NEW: Handle scale parameters and return format
  if (return_list) {
    # Build list of results
    result_list <- list(theme = bio_theme)
    
    # Add Y-axis scale if specified
    if (!is.null(y_expand)) {
      result_list$scale_y <- ggplot2::scale_y_continuous(expand = y_expand)
    }
    
    # Add X-axis scale if specified
    if (!is.null(x_expand)) {
      result_list$scale_x <- ggplot2::scale_x_continuous(expand = x_expand)
    }
    
    # Make the list work with + operator
    class(result_list) <- c("theme_bio_list", "list")
    
    return(result_list)
  } else {
    # Traditional behavior: return theme only
    return(bio_theme)
  }
}

# Custom + operator for theme_bio_list
`+.theme_bio_list` <- function(e1, e2) {
  if (inherits(e2, "theme_bio_list")) {
    # Adding two theme_bio_list objects
    result <- e1
    for (element in e2) {
      result <- result + element
    }
    return(result)
  } else if (inherits(e1, "theme_bio_list")) {
    # Adding theme_bio_list to a ggplot
    result <- e2
    for (element in e1) {
      result <- result + element
    }
    return(result)
  }
  NextMethod()
}

# Convenience function for easy usage
bio_theme_auto <- function(...) {
  theme_bio(..., return_list = TRUE)
}
