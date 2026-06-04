#' Plotting Utility Functions
#'
#' Helper functions for common plotting operations across bioRtools visualization functions.
#' These utilities ensure consistent styling and reduce code duplication in plotting code.
#'
#' @keywords internal
#' @noRd

# ============================================================================
# Theme Application
# ============================================================================

#' Apply bioRtools default theme to plot
#'
#' Adds consistent styling with bioRtools theme, including base size and family.
#'
#' @param plot ggplot2 object
#' @param base_size Base font size (default: 12)
#' @param base_family Font family (default: "Arial")
#'
#' @return ggplot2 object with theme applied
#' @keywords internal
#' @noRd
apply_bioRtools_theme <- function(plot, base_size = 12, base_family = "Arial") {
  plot +
    bioRtools::theme_bio(base_size = base_size, base_family = base_family) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = base_size + 2, face = "bold"))
}

#' Apply publication-ready theme to plot
#'
#' Adds ggprism theme for publication quality plots.
#'
#' @param plot ggplot2 object
#' @param palette Color palette name (default: "black_bold")
#'
#' @return ggplot2 object with prism theme applied
#' @keywords internal
#' @noRd
apply_prism_theme <- function(plot, palette = "black_bold") {
  plot + ggprism::theme_prism(palette = palette)
}

# ============================================================================
# Color and Aesthetic Mapping
# ============================================================================

#' Create consistent color scale for continuous data
#'
#' Generates a color gradient scale suitable for continuous variables
#' (e.g., p-values, fold changes, abundance).
#'
#' @param data Data frame or vector of values
#' @param aesthetic Aesthetic name ("colour" or "fill")
#' @param low_color Color for low values (default: "#3C5AA6")
#' @param high_color Color for high values (default: "#ED2939")
#' @param name Scale name for legend
#'
#' @return A ggplot2 scale object
#' @keywords internal
#' @noRd
create_continuous_scale <- function(data, aesthetic = "colour", low_color = "#3C5AA6",
                                    high_color = "#ED2939", name = NULL) {
  scale_fn <- if (aesthetic == "colour") {
    ggplot2::scale_colour_gradient
  } else {
    ggplot2::scale_fill_gradient
  }

  scale_fn(low = low_color, high = high_color, name = name)
}

#' Create consistent color scale for discrete groups
#'
#' Generates discrete color scales for categorical variables with
#' bioRtools color schemes.
#'
#' @param n_groups Number of groups to color
#' @param palette Palette name ("Set1", "Dark2", etc., default: "Set1")
#' @param aesthetic Aesthetic name ("colour" or "fill")
#'
#' @return A ggplot2 scale object
#' @keywords internal
#' @noRd
create_discrete_scale <- function(n_groups, palette = "Set1", aesthetic = "colour") {
  colors <- RColorBrewer::brewer.pal(min(n_groups, 9), palette)

  if (n_groups > length(colors)) {
    # Extend palette if needed
    colors <- colorRampPalette(colors)(n_groups)
  }

  scale_fn <- if (aesthetic == "colour") {
    ggplot2::scale_colour_manual
  } else {
    ggplot2::scale_fill_manual
  }

  scale_fn(values = colors)
}

# ============================================================================
# Common Plot Elements
# ============================================================================

#' Add significance annotation to plot
#'
#' Adds significance stars or brackets to a plot at specified coordinates.
#'
#' @param plot ggplot2 object
#' @param x X position for annotation
#' @param y Y position for annotation
#' @param p_value P-value to display as stars
#' @param bracket Logical, draw bracket annotation (default: TRUE)
#'
#' @return ggplot2 object with annotation
#' @keywords internal
#' @noRd
add_significance_annotation <- function(plot, x, y, p_value, bracket = TRUE) {
  if (p_value < 0.001) {
    sig_text <- "***"
  } else if (p_value < 0.01) {
    sig_text <- "**"
  } else if (p_value < 0.05) {
    sig_text <- "*"
  } else {
    sig_text <- "ns"
  }

  plot + ggplot2::annotate("text", x = x, y = y, label = sig_text,
                          size = 4, fontface = "bold")
}

#' Add grid lines to plot
#'
#' Adds consistent grid styling to plots.
#'
#' @param plot ggplot2 object
#' @param major_only Show only major grid lines (default: TRUE)
#' @param color Grid line color (default: "gray90")
#'
#' @return ggplot2 object with grid lines
#' @keywords internal
#' @noRd
add_grid_lines <- function(plot, major_only = TRUE, color = "gray90") {
  if (major_only) {
    plot + ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = color, linewidth = 0.3)
    )
  } else {
    plot + ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = color, linewidth = 0.3),
      panel.grid.minor = ggplot2::element_line(color = color, linewidth = 0.15)
    )
  }
}

# ============================================================================
# Legend Management
# ============================================================================

#' Position and style legend
#'
#' Applies consistent legend styling across plots.
#'
#' @param plot ggplot2 object
#' @param position Legend position ("right", "bottom", "top", "left")
#' @param ncol Number of columns in legend (default: 1)
#' @param title_size Title font size (default: 10)
#' @param text_size Text font size (default: 9)
#'
#' @return ggplot2 object with styled legend
#' @keywords internal
#' @noRd
style_legend <- function(plot, position = "right", ncol = 1,
                        title_size = 10, text_size = 9) {
  plot + ggplot2::theme(
    legend.position = position,
    legend.title = ggplot2::element_text(size = title_size, face = "bold"),
    legend.text = ggplot2::element_text(size = text_size),
    legend.spacing.y = ggplot2::unit(2, "mm"),
    legend.key = ggplot2::element_blank(),
    legend.background = ggplot2::element_blank(),
    legend.box = "vertical"
  ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(ncol = ncol, byrow = FALSE),
      fill = ggplot2::guide_legend(ncol = ncol, byrow = FALSE),
      shape = ggplot2::guide_legend(ncol = ncol, byrow = FALSE),
      size = ggplot2::guide_legend(ncol = ncol, byrow = FALSE)
    )
}

#' Remove or hide legend
#'
#' Removes all legends from a plot.
#'
#' @param plot ggplot2 object
#'
#' @return ggplot2 object without legend
#' @keywords internal
#' @noRd
remove_legend <- function(plot) {
  plot + ggplot2::theme(legend.position = "none")
}

# ============================================================================
# Axis Formatting
# ============================================================================

#' Format axis labels for scientific notation
#'
#' Applies scientific notation formatting to numeric axis labels.
#'
#' @param plot ggplot2 object
#' @param axis Which axis ("x", "y", or "both")
#'
#' @return ggplot2 object with formatted axes
#' @keywords internal
#' @noRd
format_scientific_axis <- function(plot, axis = "y") {
  formatter <- scales::scientific

  if (axis %in% c("x", "both")) {
    plot <- plot + ggplot2::scale_x_continuous(labels = formatter)
  }

  if (axis %in% c("y", "both")) {
    plot <- plot + ggplot2::scale_y_continuous(labels = formatter)
  }

  plot
}

#' Format axis labels with thousands separator
#'
#' Applies thousands separator formatting to numeric axis labels.
#'
#' @param plot ggplot2 object
#' @param axis Which axis ("x", "y", or "both")
#'
#' @return ggplot2 object with formatted axes
#' @keywords internal
#' @noRd
format_thousands_axis <- function(plot, axis = "y") {
  formatter <- scales::label_number(big.mark = ",")

  if (axis %in% c("x", "both")) {
    plot <- plot + ggplot2::scale_x_continuous(labels = formatter)
  }

  if (axis %in% c("y", "both")) {
    plot <- plot + ggplot2::scale_y_continuous(labels = formatter)
  }

  plot
}

# ============================================================================
# Plot Layout and Combination
# ============================================================================

#' Combine multiple plots with shared title
#'
#' Creates a grid of plots with a common title.
#'
#' @param plot_list List of ggplot2 objects
#' @param title Title for combined plot
#' @param ncol Number of columns (default: 2)
#' @param title_size Title font size (default: 14)
#'
#' @return Combined plot as patchwork object
#' @keywords internal
#' @noRd
combine_plots_with_title <- function(plot_list, title = NULL, ncol = 2,
                                     title_size = 14) {
  combined <- patchwork::wrap_plots(plot_list, ncol = ncol)

  if (!is.null(title)) {
    combined <- combined +
      patchwork::plot_annotation(title = title,
                                theme = ggplot2::theme(
                                  plot.title = ggplot2::element_text(
                                    size = title_size, face = "bold", hjust = 0.5
                                  )
                                ))
  }

  combined
}

# ============================================================================
# Data Preparation for Plotting
# ============================================================================

#' Prepare data for volcano plot
#'
#' Creates a data frame formatted for volcano plot visualization.
#'
#' @param results Results data frame with columns:
#'   log2FoldChange, padj (or pvalue)
#' @param fc_threshold Fold change threshold for coloring
#' @param p_threshold P-value threshold for coloring
#' @param use_log10_p Use -log10(p) for y-axis (default: TRUE)
#'
#' @return Data frame ready for plotting
#' @keywords internal
#' @noRd
prepare_volcano_data <- function(results, fc_threshold = 1, p_threshold = 0.05,
                                use_log10_p = TRUE) {
  df <- results %>%
    dplyr::mutate(
      p_value = if ("padj" %in% colnames(.)) padj else pvalue
    )

  if (use_log10_p) {
    df <- df %>%
      dplyr::mutate(
        neg_log10_p = -log10(p_value + 1e-300),  # Avoid log(0)
        regulation = dplyr::case_when(
          abs(log2FoldChange) >= fc_threshold & p_value < p_threshold ~ "Significant",
          TRUE ~ "Not significant"
        )
      )
  } else {
    df <- df %>%
      dplyr::mutate(
        regulation = dplyr::case_when(
          abs(log2FoldChange) >= fc_threshold & p_value < p_threshold ~ "Significant",
          TRUE ~ "Not significant"
        )
      )
  }

  df
}

#' Prepare data for MA plot (M-A plot)
#'
#' Transforms log2 fold change and abundance data for M-A plotting.
#'
#' @param results Results data frame with log2FoldChange and baseMean columns
#' @param fc_threshold Fold change threshold
#' @param p_threshold P-value threshold
#'
#' @return Data frame with M and A values
#' @keywords internal
#' @noRd
prepare_ma_plot_data <- function(results, fc_threshold = 1, p_threshold = 0.05) {
  results %>%
    dplyr::mutate(
      M = log2FoldChange,
      A = log10(baseMean + 1),
      p_value = if ("padj" %in% colnames(.)) padj else pvalue,
      regulation = dplyr::case_when(
        abs(M) >= fc_threshold & p_value < p_threshold & M > 0 ~ "Up-regulated",
        abs(M) >= fc_threshold & p_value < p_threshold & M < 0 ~ "Down-regulated",
        TRUE ~ "Not significant"
      )
    )
}
