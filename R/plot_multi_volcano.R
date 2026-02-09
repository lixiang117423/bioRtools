#' Plot Multi-Group Volcano Plot
#'
#' @description
#' Create a multi-group volcano plot to visualize differential expression
#' analysis results across multiple clusters or groups. The plot displays
#' fold change values on the y-axis with different groups shown as separate
#' columns on the x-axis. Points are color-coded by up/down regulation and
#' significant genes can be labeled.
#'
#' @param data A data frame containing the volcano plot data with columns
#'   specified by \code{fc_column}, \code{cluster_column}, \code{gene_column},
#'   and \code{pval_column}.
#' @param fc_column Character string. Name of the column containing fold change
#'   values (typically log2 fold change). Default is "avg_log2FC".
#' @param cluster_column Character string. Name of the column containing group
#'   or cluster identifiers. Default is "cluster".
#' @param gene_column Character string. Name of the column containing gene
#'   names or identifiers. Default is "gene".
#' @param pval_column Character string. Name of the column containing adjusted
#'   p-values. Default is "p_val_adj".
#' @param fc_threshold Numeric value. Absolute fold change threshold for
#'   labeling significant genes. Default is 1.
#' @param pval_threshold Numeric value between 0 and 1. P-value threshold for
#'   filtering genes to label. Default is 0.05.
#' @param label_n Integer. Number of top up-regulated and down-regulated genes
#'   to label per cluster. Default is 5.
#' @param cluster_colors Character vector. Color palette for cluster
#'   background strips. Default uses scientific color palette.
#' @param up_color Character string. Color for up-regulated genes.
#'   Default is "#0073C2FF" (blue).
#' @param down_color Character string. Color for down-regulated genes.
#'   Default is "#EE0000FF" (red).
#' @param y_limits Numeric vector of length 2. Y-axis limits as c(min, max).
#'   Default is c(-3, 7).
#' @param y_breaks Numeric vector. Y-axis break points. Default is
#'   c(-3, -2, -1, 0, 2, 4, 6).
#' @param legend_position Character string or numeric vector. Position of the
#'   legend. Default is c(0.08, 0.9) for top-left.
#'
#' @return A ggplot object that can be further customized or saved using
#'   \code{ggplot2::ggsave()}.
#'
#' @details
#' The function creates a volcano plot with the following features:
#' \itemize{
#'   \item Each cluster/group is shown as a separate column
#'   \item Background strips indicate cluster boundaries
#'   \item Points are color-coded by regulation direction
#'   \item Significant genes are labeled based on fold change and p-value
#'   \item Cluster names are displayed in the center of each column
#' }
#'
#' Genes are labeled based on the following criteria:
#' \enumerate{
#'   \item Adjusted p-value < \code{pval_threshold}
#'   \item Absolute fold change > \code{fc_threshold}
#'   \item Top \code{label_n} up-regulated and down-regulated genes per cluster
#' }
#'
#' @note
#' \itemize{
#'   \item The function requires \code{ggrepel} and \code{ggprism} packages
#'   \item Cluster order is preserved from the input data
#'   \item Labels may overlap if many genes meet the criteria
#'   \item Adjust \code{label_n} to control the number of labels
#' }
#'
#' @examples
#' \dontrun{
#' library(bioRtools)
#' library(dplyr)
#' library(readr)
#'
#' # Read data
#' df <- read_tsv("data.tsv")
#'
#' # Basic multi-group volcano plot
#' p <- plot_multi_volcano(df)
#' print(p)
#'
#' # Customize with different thresholds
#' p_custom <- plot_multi_volcano(
#'   df,
#'   fc_threshold = 1.5,
#'   pval_threshold = 0.01,
#'   label_n = 3
#' )
#' print(p_custom)
#'
#' # Save the plot
#' ggplot2::ggsave("multi_volcano.png", p, width = 10, height = 6, dpi = 300)
#'
#' # Custom colors and axis limits
#' p_custom2 <- plot_multi_volcano(
#'   df,
#'   cluster_colors = c("#E41A1C", "#377EB8", "#4DAF4A"),
#'   y_limits = c(-5, 10),
#'   up_color = "#FF6B6B",
#'   down_color = "#6B6BFF"
#' )
#' print(p_custom2)
#' }
#'
#' @export
plot_multi_volcano <- function(data,
                                fc_column = "avg_log2FC",
                                cluster_column = "cluster",
                                gene_column = "gene",
                                pval_column = "p_val_adj",
                                fc_threshold = 1,
                                pval_threshold = 0.05,
                                label_n = 5,
                                cluster_colors = c("#3B9AB2", "#78B7C5",
                                                  "#EBCC2A", "#E1AF00",
                                                  "#F21A00", "#C51B7D",
                                                  "#7F3B08", "#B2ABD2",
                                                  "#ABDDA4", "#FC8D62"),
                                up_color = "#0073C2FF",
                                down_color = "#EE0000FF",
                                y_limits = c(-3, 7),
                                y_breaks = c(-3, -2, -1, 0, 2, 4, 6),
                                legend_position = c(0.08, 0.9)) {

  # Input validation
  validate_multi_volcano_input(
    data,
    fc_column,
    cluster_column,
    gene_column,
    pval_column,
    fc_threshold,
    pval_threshold,
    label_n,
    cluster_colors,
    up_color,
    down_color,
    y_limits,
    y_breaks,
    legend_position
  )

  # Prepare data
  df <- prepare_multi_volcano_data(
    data,
    fc_column,
    cluster_column
  )

  # Create background strips
  bg_df <- create_background_strips(df, cluster_column)

  # Create vertical background strips
  bg_vertical <- create_vertical_backgrounds(df, fc_column, cluster_column)

  # Select genes to label
  label_df <- select_genes_to_label(
    df,
    fc_column,
    cluster_column,
    gene_column,
    pval_column,
    fc_threshold,
    pval_threshold,
    label_n
  )

  # Build the plot
  p <- build_multi_volcano_plot(
    df,
    label_df,
    bg_df,
    bg_vertical,
    fc_column,
    cluster_column,
    gene_column,
    cluster_colors,
    up_color,
    down_color,
    y_limits,
    y_breaks,
    legend_position
  )

  p
}


#' Validate Input for Multi-Group Volcano Plot
#'
#' @param data Input data frame
#' @param fc_column Fold change column name
#' @param cluster_column Cluster column name
#' @param gene_column Gene column name
#' @param pval_column P-value column name
#' @param fc_threshold Fold change threshold
#' @param pval_threshold P-value threshold
#' @param label_n Number of labels
#' @param cluster_colors Cluster colors
#' @param up_color Up-regulated color
#' @param down_color Down-regulated color
#' @param y_limits Y-axis limits
#' @param y_breaks Y-axis breaks
#' @param legend_position Legend position
#'
#' @return TRUE if validation passes, otherwise stops with error
#' @keywords internal
#' @noRd
validate_multi_volcano_input <- function(data,
                                          fc_column,
                                          cluster_column,
                                          gene_column,
                                          pval_column,
                                          fc_threshold,
                                          pval_threshold,
                                          label_n,
                                          cluster_colors,
                                          up_color,
                                          down_color,
                                          y_limits,
                                          y_breaks,
                                          legend_position) {

  # Check data frame
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame", call. = FALSE)
  }

  # Check required columns
  required_cols <- c(cluster_column, fc_column, gene_column, pval_column)
  missing_cols <- setdiff(required_cols, names(data))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  # Validate thresholds
  if (!is.numeric(fc_threshold) || length(fc_threshold) != 1 || fc_threshold <= 0) {
    stop("'fc_threshold' must be a single positive numeric value", call. = FALSE)
  }

  if (!is.numeric(pval_threshold) || length(pval_threshold) != 1 ||
      pval_threshold <= 0 || pval_threshold >= 1) {
    stop("'pval_threshold' must be a single numeric value between 0 and 1",
         call. = FALSE)
  }

  if (!is.numeric(label_n) || length(label_n) != 1 || label_n < 1 ||
      label_n != as.integer(label_n)) {
    stop("'label_n' must be a positive integer", call. = FALSE)
  }

  # Validate colors
  if (!is.character(cluster_colors) || length(cluster_colors) < 1) {
    stop("'cluster_colors' must be a character vector", call. = FALSE)
  }

  if (!is.character(up_color) || length(up_color) != 1) {
    stop("'up_color' must be a single character string", call. = FALSE)
  }

  if (!is.character(down_color) || length(down_color) != 1) {
    stop("'down_color' must be a single character string", call. = FALSE)
  }

  # Validate axis parameters
  if (!is.numeric(y_limits) || length(y_limits) != 2) {
    stop("'y_limits' must be a numeric vector of length 2", call. = FALSE)
  }

  if (!is.numeric(y_breaks) || length(y_breaks) < 1) {
    stop("'y_breaks' must be a numeric vector", call. = FALSE)
  }

  # Validate legend position
  if (!is.character(legend_position) && !is.numeric(legend_position)) {
    stop("'legend_position' must be a character string or numeric vector",
         call. = FALSE)
  }

  if (is.numeric(legend_position) && length(legend_position) != 2) {
    stop("'legend_position' must be a numeric vector of length 2 when numeric",
         call. = FALSE)
  }

  invisible(TRUE)
}


#' Prepare Data for Multi-Group Volcano Plot
#'
#' @param data Input data frame
#' @param fc_column Fold change column name
#' @param cluster_column Cluster column name
#'
#' @return Prepared data frame with chr and type columns
#' @keywords internal
#' @noRd
prepare_multi_volcano_data <- function(data, fc_column, cluster_column) {
  df <- data %>%
    dplyr::mutate(
      chr = as.character(!!rlang::sym(cluster_column)),
      type = dplyr::if_else(
        !!rlang::sym(fc_column) > 0,
        "UP_Highly",
        "Down_Highly"
      )
    )

  # Set cluster factor order
  df[[cluster_column]] <- factor(
    df[[cluster_column]],
    levels = unique(df[[cluster_column]])
  )

  df
}


#' Create Background Strips for Clusters
#'
#' @param df Data frame
#' @param cluster_column Cluster column name
#'
#' @return Data frame with background strip coordinates
#' @keywords internal
#' @noRd
create_background_strips <- function(df, cluster_column) {
  tibble::tibble(
    cluster = levels(df[[cluster_column]]),
    xmin = seq_along(levels(df[[cluster_column]])) - 0.48,
    xmax = seq_along(levels(df[[cluster_column]])) + 0.48,
    ymin = -0.5,
    ymax = 0.5
  )
}


#' Create Vertical Background Strips
#'
#' @param df Data frame
#' @param fc_column Fold change column name
#' @param cluster_column Cluster column name
#'
#' @return Data frame with vertical background strip coordinates
#' @keywords internal
#' @noRd
create_vertical_backgrounds <- function(df, fc_column, cluster_column) {
  df %>%
    dplyr::group_by(!!rlang::sym(cluster_column)) %>%
    dplyr::summarise(
      ymin = min(!!rlang::sym(fc_column), na.rm = TRUE) - 0.1,
      ymax = max(!!rlang::sym(fc_column), na.rm = TRUE) + 0.1,
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      xmin = dplyr::row_number() - 0.48,
      xmax = dplyr::row_number() + 0.48
    )
}


#' Select Genes to Label
#'
#' @param df Data frame
#' @param fc_column Fold change column name
#' @param cluster_column Cluster column name
#' @param gene_column Gene column name
#' @param pval_column P-value column name
#' @param fc_threshold Fold change threshold
#' @param pval_threshold P-value threshold
#' @param label_n Number of labels per cluster
#'
#' @return Data frame with genes to label
#' @keywords internal
#' @noRd
select_genes_to_label <- function(df,
                                  fc_column,
                                  cluster_column,
                                  gene_column,
                                  pval_column,
                                  fc_threshold,
                                  pval_threshold,
                                  label_n) {

  df %>%
    dplyr::filter(
      !!rlang::sym(pval_column) < pval_threshold,
      abs(!!rlang::sym(fc_column)) > fc_threshold
    ) %>%
    dplyr::group_by(chr) %>%
    dplyr::group_modify(~ {
      dplyr::bind_rows(
        dplyr::slice_max(.x, !!rlang::sym(fc_column), n = label_n),
        dplyr::slice_min(.x, !!rlang::sym(fc_column), n = label_n)
      )
    }) %>%
    dplyr::ungroup()
}


#' Build Multi-Group Volcano Plot
#'
#' @param df Data frame
#' @param label_df Data frame with labels
#' @param bg_df Background strips data frame
#' @param bg_vertical Vertical background strips data frame
#' @param fc_column Fold change column name
#' @param cluster_column Cluster column name
#' @param gene_column Gene column name
#' @param cluster_colors Color palette
#' @param up_color Up-regulated color
#' @param down_color Down-regulated color
#' @param y_limits Y-axis limits
#' @param y_breaks Y-axis breaks
#' @param legend_position Legend position
#'
#' @return ggplot object
#' @keywords internal
#' @noRd
build_multi_volcano_plot <- function(df,
                                      label_df,
                                      bg_df,
                                      bg_vertical,
                                      fc_column,
                                      cluster_column,
                                      gene_column,
                                      cluster_colors,
                                      up_color,
                                      down_color,
                                      y_limits,
                                      y_breaks,
                                      legend_position) {

  ggplot2::ggplot(df, ggplot2::aes(.data[[cluster_column]],
                                     .data[[fc_column]],
                                     color = .data$type)) +
    # Vertical background strips
    ggplot2::geom_rect(
      data = bg_vertical,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = "grey95"
    ) +
    # Cluster background strips
    ggplot2::geom_rect(
      data = bg_df,
      ggplot2::aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax,
        fill = cluster
      ),
      inherit.aes = FALSE
    ) +
    # Scatter points
    ggplot2::geom_jitter(stroke = 0) +
    # Gene labels
    ggrepel::geom_text_repel(
      data = label_df,
      ggplot2::aes(
        .data[[cluster_column]],
        .data[[fc_column]],
        label = .data[[gene_column]]
      ),
      size = 2.5,
      color = "black",
      box.padding = 0.2
    ) +
    # Cluster labels
    ggplot2::geom_text(
      ggplot2::aes(.data[[cluster_column]], 0, label = chr),
      size = 3,
      color = "white",
      show.legend = FALSE
    ) +
    # Color scales
    ggplot2::scale_fill_manual(
      values = cluster_colors,
      guide = "none"
    ) +
    ggplot2::scale_color_manual(
      values = c("UP_Highly" = up_color, "Down_Highly" = down_color)
    ) +
    # Y-axis scale
    ggplot2::scale_y_continuous(
      limits = y_limits,
      breaks = y_breaks,
      guide = "prism_offset_minor"
    ) +
    # Labels
    ggplot2::labs(x = NULL, y = "average log2FC") +
    # Legend
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(size = 5, shape = 19)
      )
    ) +
    # Theme
    ggprism::theme_prism(base_line_size = 0.3) +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 10),
      axis.title = ggplot2::element_text(size = 11),
      legend.position = legend_position,
      legend.key = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(margin = ggplot2::margin(l = 0))
    )
}
