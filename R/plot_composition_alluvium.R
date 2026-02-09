#' Plot Taxonomic Composition as Alluvial Diagram
#'
#' @description
#' Create an alluvial diagram (Sankey-like plot) to visualize taxonomic
#' composition across multiple sample groups. This function is ideal for
#' showing how the relative abundance of different taxa changes across
#' conditions or time points, commonly used in microbiome or virome studies.
#'
#' @param data A data frame in wide format with one taxonomic/identifier column
#'   and multiple sample columns containing relative abundances.
#' @param taxon_column Character string. Name of the column containing taxonomic
#'   identifiers (e.g., "family", "genus", "species"). Default is "family".
#' @param id_column Character string. Alternative name for taxonomic column when
#'   \code{taxon_column} is not found. Default is NULL.
#' @param palette Character string or vector. Color palette for the alluvium.
#'   Can be a single RColorBrewer palette name or a vector of colors.
#'   Default is "Paired".
#' @param alpha Numeric value between 0 and 1 for transparency. Default is 0.7.
#' @param y_labels_format Logical. Whether to format y-axis labels as percentages.
#'   Default is TRUE.
#' @param y_accuracy Numeric value. Precision for y-axis percentage labels.
#'   Default is 1 (i.e., 1 decimal place).
#' @param y_expand Numeric vector of length 2. Y-axis expansion.
#'   Default is c(0, 0) for no expansion.
#' @param legend_position Character string. Position of the legend.
#'   Options: "right", "left", "top", "bottom", "none". Default is "right".
#' @param legend_key_height Numeric value. Height of legend keys. Default is 0.4.
#' @param legend_key_width Numeric value. Width of legend keys. Default is 0.6.
#' @param base_size Numeric value. Base font size for text. Default is 11.
#'
#' @return A ggplot object displaying the alluvial diagram.
#'
#' @details
#' The function performs the following operations:
#' \enumerate{
#'   \item Validates input data format
#'   \item Calculates row sums to determine ordering
#'   \item Converts data from wide to long format
#'   \item Creates alluvial diagram with stratums and flows
#'   \item Applies professional styling
#' }
#'
#' The alluvial diagram shows:
#' \itemize{
#'   \item Stratums: Taxa on the y-axis (ordered by total abundance)
#'   \item Flows: Connections showing relative abundance changes across groups
#'   \item Colors: Distinct colors for each taxon
#' }
#'
#' @note
#' \itemize{
#'   \item Requires \code{ggalluvial} package
#'   \item Y-axis values are automatically scaled to percentages
#'   \item Taxa are ordered by total abundance across all samples
#'   \item Input data should contain relative abundances (not raw counts)
#' }
#'
#' @examples
#' \dontrun{
#' library(bioRtools)
#' library(readxl)
#' library(dplyr)
#'
#' # Read data from Excel file
#' df <- read_excel("data.xlsx", sheet = "Fig.1a", skip = 1)
#'
#' # Basic alluvial plot
#' p <- plot_composition_alluvium(
#'   df,
#'   taxon_column = "Viral Taxonomy (family)"
#' )
#' print(p)
#'
#' # Custom color palette and formatting
#' p_custom <- plot_composition_alluvium(
#'   df,
#'   taxon_column = "Viral Taxonomy (family)",
#'   palette = "Set2",
#'   alpha = 0.8,
#'   y_accuracy = 0.1
#' )
#' print(p_custom)
#'
#' # Custom colors vector
#' p_custom2 <- plot_composition_alluvium(
#'   df,
#'   taxon_column = "Viral Taxonomy (family)",
#'   palette = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
#' )
#' print(p_custom2)
#' }
#'
#' @author Xiang LI \email{lixiang117423@@foxmail.com}
#'
#' @export
plot_composition_alluvium <- function(data,
                                        taxon_column = "family",
                                        id_column = NULL,
                                        palette = "Paired",
                                        alpha = 0.7,
                                        y_labels_format = TRUE,
                                        y_accuracy = 1,
                                        y_expand = c(0, 0),
                                        legend_position = "right",
                                        legend_key_height = 0.4,
                                        legend_key_width = 0.6,
                                        base_size = 11) {

  # Input validation
  validate_alluvium_input(data, taxon_column, id_column, palette, alpha)

  # Prepare data
  df_clean <- clean_alluvium_data(data, taxon_column, id_column)

  # Build the plot
  p <- build_alluvium_plot(
    df_clean,
    palette,
    alpha,
    y_labels_format,
    y_accuracy,
    y_expand,
    legend_position,
    legend_key_height,
    legend_key_width,
    base_size
  )

  p
}


#' Validate Input for Alluvium Plot
#'
#' @param data Input data frame
#' @param taxon_column Taxon column name
#' @param id_column Alternative id column name
#' @param palette Color palette
#' @param alpha Transparency value
#'
#' @return TRUE if validation passes, otherwise stops with error
#' @keywords internal
#' @noRd
validate_alluvium_input <- function(data,
                                     taxon_column,
                                     id_column,
                                     palette,
                                     alpha) {

  # Check data frame
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame", call. = FALSE)
  }

  # Check taxon column
  if (!taxon_column %in% names(data)) {
    # Try alternative column name
    if (!is.null(id_column) && id_column %in% names(data)) {
      taxon_column <- id_column
    } else {
      stop(
        "Taxon column '", taxon_column, "' not found in data. ",
        "Available columns: ", paste(names(data), collapse = ", "),
        call. = FALSE
      )
    }
  }

  # Check for numeric columns (at least 2 needed)
  numeric_cols <- names(data)[sapply(data, is.numeric)]
  if (length(numeric_cols) < 2) {
    stop(
      "Data must contain at least 2 numeric columns with abundance values",
      call. = FALSE
    )
  }

  # Validate alpha
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha < 0 || alpha > 1) {
    stop("'alpha' must be a single numeric value between 0 and 1", call. = FALSE)
  }

  # Validate palette
  if (is.character(palette) && length(palette) == 1) {
    # It's a palette name, ggplot2 will validate it
  } else if (!is.character(palette)) {
    stop("'palette' must be a character vector or palette name", call. = FALSE)
  }

  invisible(TRUE)
}


#' Clean and Prepare Alluvium Data
#'
#' @param data Input data frame
#' @param taxon_column Taxon column name
#' @param id_column Alternative id column name
#'
#' @return Cleaned data frame in long format
#' @keywords internal
#' @noRd
clean_alluvium_data <- function(data, taxon_column, id_column) {

  # Determine actual taxon column name
  if (!taxon_column %in% names(data) && !is.null(id_column)) {
    actual_taxon_col <- id_column
  } else {
    actual_taxon_col <- taxon_column
  }

  # Rename taxon column to standard name
  df_clean <- data %>%
    dplyr::rename(family = !!rlang::sym(actual_taxon_col))

  # Calculate row sums for ordering (sum across all numeric columns except family)
  numeric_cols <- df_clean %>%
    dplyr::select(-family) %>%
    dplyr::select(where(is.numeric)) %>%
    names()

  df_clean <- df_clean %>%
    dplyr::mutate(row_sum = rowSums(dplyr::pick(all_of(numeric_cols)), na.rm = TRUE))

  # Order by row sum
  df_clean <- df_clean %>%
    dplyr::mutate(
      family = factor(
        family,
        levels = family[order(row_sum, decreasing = TRUE)]
      )
    ) %>%
    dplyr::select(-row_sum) %>%
    tidyr::pivot_longer(
      cols = -family,
      names_to = "name",
      values_to = "value"
    )

  df_clean
}


#' Build Alluvium Plot
#'
#' @param df_clean Cleaned data frame in long format
#' @param palette Color palette
#' @param alpha Transparency value
#' @param y_labels_format Whether to format y-axis as percentages
#' @param y_accuracy Precision for y-axis labels
#' @param y_expand Y-axis expansion
#' @param legend_position Legend position
#' @param legend_key_height Legend key height
#' @param legend_key_width Legend key width
#' @param base_size Base font size
#'
#' @return ggplot object
#' @keywords internal
#' @noRd
build_alluvium_plot <- function(df_clean,
                                palette,
                                alpha,
                                y_labels_format,
                                y_accuracy,
                                y_expand,
                                legend_position,
                                legend_key_height,
                                legend_key_width,
                                base_size) {

  # Build y-axis scale
  if (y_labels_format) {
    y_scale <- ggplot2::scale_y_continuous(
      expand = y_expand,
      labels = scales::percent_format(accuracy = y_accuracy)
    )
  } else {
    y_scale <- ggplot2::scale_y_continuous(expand = y_expand)
  }

  # Build the plot
  p <- ggplot2::ggplot(
    df_clean,
    ggplot2::aes(
      x = name,
      y = value,
      alluvium = family,
      stratum = family
    )
  ) +
    ggalluvial::geom_alluvium(
      ggplot2::aes(fill = family, color = family),
      alpha = alpha
    ) +
    ggalluvial::geom_stratum(
      ggplot2::aes(fill = family, color = family)
    ) +
    ggplot2::scale_fill_brewer(palette = palette) +
    ggplot2::scale_color_brewer(palette = palette) +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    y_scale +
    ggplot2::labs(y = "Relative abundance (%)", x = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.background = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank(),
      legend.key.height = ggplot2::unit(legend_key_height, "cm"),
      legend.key.width = ggplot2::unit(legend_key_width, "cm"),
      legend.key.spacing.y = ggplot2::unit(0.1, "cm"),
      legend.position = legend_position,
      axis.text = ggplot2::element_text(color = "black", size = base_size)
    )

  p
}
