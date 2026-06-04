#' Biotools Color Palettes and Scales
#'
#' Standardized color palettes for bioRtools visualization functions.
#' Ensures visual consistency across all plots and provides publication-quality colors.
#'
#' @details
#' Available palettes:
#' - `pal_regulation`: For differential expression/abundance (up, down, not significant)
#' - `pal_groups`: For categorical grouping variables
#' - `pal_significance`: For p-value significance levels
#' - `pal_abundance`: For abundance/count visualization
#' - `pal_heatmap`: For heatmap annotations
#'
#' @keywords internal
#' @noRd

# ============================================================================
# Core Color Definitions
# ============================================================================

#' Regulation direction colors (Up, Down, Not Significant)
#'
#' @param alpha Transparency value (0-1), default: 1 (opaque)
#'
#' @return Named character vector with colors for regulation categories
#' @keywords internal
#' @noRd
pal_regulation <- function(alpha = 1) {
  colors <- c(
    "Up-regulated" = "#ED2939",           # Red
    "Down-regulated" = "#3C5AA6",         # Blue
    "Enriched" = "#ED2939",               # Red (microbiome equivalent)
    "Depleted" = "#3C5AA6",               # Blue
    "Not significant" = "#CCCCCC",        # Gray
    "NS" = "#CCCCCC"                      # Gray (abbreviated)
  )

  if (alpha < 1) {
    # Add transparency
    colors <- scales::alpha(colors, alpha)
  }

  colors
}

#' Significance level colors (p-value thresholds)
#'
#' @return Named character vector with colors for significance levels
#' @keywords internal
#' @noRd
pal_significance <- function() {
  c(
    "****" = "#006400",     # Dark green (p < 0.0001)
    "***" = "#228B22",      # Forest green (p < 0.001)
    "**" = "#FFD700",       # Gold (p < 0.01)
    "*" = "#FF8C00",        # Dark orange (p < 0.05)
    "NS" = "#CCCCCC"        # Gray (not significant)
  )
}

#' Categorical group colors
#'
#' Uses a publication-ready palette suitable for up to 12 groups.
#'
#' @param n_colors Number of colors needed (default: 8)
#' @param palette Palette name: "primary" (default), "pastel", "dark"
#'
#' @return Character vector of colors
#' @keywords internal
#' @noRd
pal_groups <- function(n_colors = 8, palette = "primary") {
  palettes <- list(
    primary = c(
      "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
      "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
      "#bcbd22", "#17becf", "#aec7e8", "#ffbb78"
    ),
    pastel = c(
      "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
      "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
      "#cab2d6", "#6a3d9a", "#ffff99", "#b15928"
    ),
    dark = c(
      "#001f3f", "#0066cc", "#003d99", "#004d99",
      "#003366", "#002b4d", "#1a1a4d", "#4d0000",
      "#660000", "#990000", "#cc0000", "#ff0000"
    )
  )

  selected_palette <- palettes[[palette]] %||% palettes$primary
  selected_palette[1:min(n_colors, length(selected_palette))]
}

#' Abundance visualization colors
#'
#' For visualizing relative abundance or count data.
#'
#' @param type "sequential" (low to high abundance) or "diverging" (low-neutral-high)
#'
#' @return Color palette function or vector
#' @keywords internal
#' @noRd
pal_abundance <- function(type = "sequential") {
  if (type == "sequential") {
    # Yellow to dark blue gradient
    c("#FFF7BC", "#FED976", "#FD8D3C", "#E31A1C", "#800026")
  } else if (type == "diverging") {
    # Blue-white-red diverging scale for centered data
    c("#3C5AA6", "#8FA3D1", "#FFFFFF", "#E5B4B8", "#ED2939")
  }
}

# ============================================================================
# Heatmap Annotation Colors
# ============================================================================

#' Colors for heatmap row/column annotations
#'
#' @param annotation_type Type of annotation: "binary", "groups", or "gradient"
#' @param n_levels Number of levels (for categorical annotations)
#'
#' @return Character vector of colors
#' @keywords internal
#' @noRd
pal_heatmap_annotation <- function(annotation_type = "groups", n_levels = NULL) {
  if (annotation_type == "binary") {
    c("yes" = "#2ca02c", "no" = "#d62728", "1" = "#2ca02c", "0" = "#d62728")
  } else if (annotation_type == "groups" && !is.null(n_levels)) {
    pal_groups(n_levels)
  } else if (annotation_type == "gradient") {
    # Returns a function for continuous gradients
    colorRampPalette(c("#ffffff", "#fee5d9", "#fcae91", "#fb6a4a", "#cb181d", "#99000d"))
  }
}

# ============================================================================
# ggplot2 Scale Wrappers
# ============================================================================

#' Apply regulation color scale to ggplot2
#'
#' Convenience wrapper for regulation color mapping.
#'
#' @param aesthetic Aesthetic to map ("colour", "fill", "color")
#' @param alpha Transparency (0-1)
#' @param guide Guide type ("legend", "none")
#'
#' @return ggplot2 scale object
#' @keywords internal
#' @noRd
scale_regulation <- function(aesthetic = "colour", alpha = 1, guide = "legend") {
  colors <- pal_regulation(alpha)

  if (aesthetic %in% c("colour", "color")) {
    ggplot2::scale_colour_manual(
      values = colors,
      guide = guide,
      name = "Regulation",
      na.value = "#CCCCCC"
    )
  } else if (aesthetic == "fill") {
    ggplot2::scale_fill_manual(
      values = colors,
      guide = guide,
      name = "Regulation",
      na.value = "#CCCCCC"
    )
  }
}

#' Apply significance color scale to ggplot2
#'
#' @param aesthetic Aesthetic to map ("colour", "fill", "color")
#' @param guide Guide type ("legend", "none")
#'
#' @return ggplot2 scale object
#' @keywords internal
#' @noRd
scale_significance <- function(aesthetic = "colour", guide = "legend") {
  colors <- pal_significance()

  if (aesthetic %in% c("colour", "color")) {
    ggplot2::scale_colour_manual(
      values = colors,
      guide = guide,
      name = "Significance",
      na.value = "#CCCCCC"
    )
  } else if (aesthetic == "fill") {
    ggplot2::scale_fill_manual(
      values = colors,
      guide = guide,
      name = "Significance",
      na.value = "#CCCCCC"
    )
  }
}

#' Apply group color scale to ggplot2
#'
#' @param n_groups Number of groups
#' @param palette Palette type ("primary", "pastel", "dark")
#' @param aesthetic Aesthetic to map ("colour", "fill", "color")
#' @param guide Guide type ("legend", "none")
#'
#' @return ggplot2 scale object
#' @keywords internal
#' @noRd
scale_groups <- function(n_groups, palette = "primary", aesthetic = "colour",
                        guide = "legend") {
  colors <- pal_groups(n_groups, palette)

  if (aesthetic %in% c("colour", "color")) {
    ggplot2::scale_colour_manual(
      values = colors,
      guide = guide,
      name = "Group",
      na.value = "#CCCCCC"
    )
  } else if (aesthetic == "fill") {
    ggplot2::scale_fill_manual(
      values = colors,
      guide = guide,
      name = "Group",
      na.value = "#CCCCCC"
    )
  }
}

# ============================================================================
# Color Utilities
# ============================================================================

#' Get hex color code for a specific regulation/significance category
#'
#' @param category Category name (e.g., "Up-regulated", "Down-regulated", "NS")
#' @param alpha Transparency (0-1)
#'
#' @return Hex color code
#' @keywords internal
#' @noRd
get_color <- function(category, alpha = 1) {
  all_colors <- c(
    pal_regulation(alpha),
    pal_significance(),
    pal_groups(12)
  )

  all_colors[category] %||% "#CCCCCC"
}

#' Check color contrast ratio (for accessibility)
#'
#' @param color1 First color (hex or name)
#' @param color2 Second color (hex or name)
#'
#' @return Numeric contrast ratio (should be > 4.5 for good contrast)
#' @keywords internal
#' @noRd
check_color_contrast <- function(color1, color2) {
  # Convert to RGB
  rgb1 <- col2rgb(color1) / 255
  rgb2 <- col2rgb(color2) / 255

  # Calculate relative luminance
  lum <- function(rgb) {
    rgb <- ifelse(rgb <= 0.03928, rgb / 12.92, ((rgb + 0.055) / 1.055)^2.4)
    0.2126 * rgb[1] + 0.7152 * rgb[2] + 0.0722 * rgb[3]
  }

  l1 <- lum(rgb1)
  l2 <- lum(rgb2)

  # Calculate contrast ratio
  (max(l1, l2) + 0.05) / (min(l1, l2) + 0.05)
}

# ============================================================================
# Default Theme Colors
# ============================================================================

#' Default color scheme for bioRtools
#'
#' @return List with default colors for various plot elements
#' @keywords internal
#' @noRd
bioRtools_color_scheme <- function() {
  list(
    primary = "#1f77b4",
    secondary = "#ff7f0e",
    success = "#2ca02c",
    danger = "#d62728",
    warning = "#ff8c00",
    info = "#0099ff",
    light = "#f0f0f0",
    dark = "#333333",
    background = "#ffffff",
    grid = "#e0e0e0"
  )
}
