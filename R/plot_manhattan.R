#' Plot Manhattan plot from GWAS results or genomic data
#'
#' @description
#' A comprehensive function for creating Manhattan plots from GWAS results or other
#' genome-wide data. Supports both data frames and file input, with customizable
#' visualization options. Optionally displays an rMVP-style SNP density panel below
#' the main plot.
#'
#' @param data A data frame containing the plotting data. If \code{input_file} is
#'   provided, this parameter can be omitted.
#' @param input_file Character string specifying the path to the input data file.
#'   If \code{data} is provided, this parameter can be omitted. Supports common
#'   delimited formats (csv, tsv, txt).
#' @param chr_col Character string specifying the chromosome column name.
#'   Default is "chr".
#' @param pos_col Character string specifying the physical position column name.
#'   Default is "pos".
#' @param val_col Character string specifying the value column name (e.g., p-values).
#'   Default is "value".
#' @param transform_log10 Logical value indicating whether to apply -log10
#'   transformation to the values. Set to TRUE for p-values. Default is FALSE.
#' @param title Character string for the plot title. Default is "Manhattan Plot".
#' @param ylab Character string for the Y-axis label. If NULL, will be set
#'   automatically based on \code{transform_log10}.
#' @param colors Character vector of colors for alternating chromosome coloring.
#'   Default is c("#D85A30", "#1D9E75").
#' @param threshold_line Numeric value for drawing a horizontal threshold line.
#'   Default is NULL (no line).
#' @param threshold_color Character string specifying the threshold line color.
#'   Default is "red".
#' @param point_size Numeric value for point size. Default is 1.2.
#' @param point_alpha Numeric value for point transparency (0-1). Default is 0.8.
#' @param show_density Logical. Whether to display an rMVP-style SNP density panel
#'   below the main Manhattan plot. Requires the \code{patchwork} package.
#'   Default is FALSE.
#' @param density_colors Character vector of 3 colors used to interpolate the
#'   density gradient from low to high SNP density.
#'   Default is c("darkgreen", "yellow", "red"), matching rMVP's style.
#' @param bin_size Integer. Window size (in bp) for counting SNPs when computing
#'   density. Smaller values give finer resolution. Default is 1e5 (100 kb).
#' @param density_height Numeric. Relative height of the density panel compared
#'   to the main plot (passed to \code{patchwork::plot_layout}). Default is 0.12.
#' @param show_density_legend Logical. Whether to show a color legend for the
#'   SNP density panel. Default is FALSE.
#'
#' @return A list containing:
#' \describe{
#'   \item{plot.manhattan}{The main ggplot Manhattan plot object.}
#'   \item{plot.density}{The density ggplot object (NULL if show_density = FALSE).}
#'   \item{plot.combined}{The patchwork-combined plot (NULL if show_density = FALSE).}
#'   \item{data.processed}{The processed data frame used for plotting.}
#'   \item{chromosome.centers}{A data frame of chromosome center positions.}
#' }
#'
#' @details
#' When \code{show_density = TRUE}, the function computes per-bin SNP counts across
#' the genome and renders them as a color-coded density strip beneath the main plot,
#' mimicking the style produced by \code{rMVP::MVP.Report}. The combined figure is
#' assembled with \code{patchwork}.
#'
#' @note
#' Required packages: dplyr, ggplot2, readr, gtools, rlang.
#' Additional package for density panel: patchwork.
#'
#' @export
plot_manhattan <- function(data = NULL,
                           input_file = NULL,
                           chr_col = "chr",
                           pos_col = "pos",
                           val_col = "value",
                           transform_log10 = FALSE,
                           title = "Manhattan Plot",
                           ylab = NULL,
                           colors = pal_groups(2, palette = "manhattan"),
                           threshold_line = NULL,
                           threshold_color = "red",
                           point_size = 1.2,
                           point_alpha = 0.8,
                           # --- density panel parameters ---
                           show_density = FALSE,
                           density_colors = c("darkgreen", "yellow", "red"),
                           bin_size = 1e5,
                           density_height = 0.12,
                           show_density_legend = FALSE) {

  # ── Input validation ────────────────────────────────────────────────────────
  if (is.null(data) && is.null(input_file)) {
    stop("Either 'data' or 'input_file' must be provided")
  }
  if (!is.null(data) && !is.null(input_file)) {
    warning("Both 'data' and 'input_file' provided. Using 'data' and ignoring 'input_file'")
  }

  # ── Package checks ───────────────────────────────────────────────────────────
  required_packages <- c("dplyr", "ggplot2", "gtools", "rlang")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop("Required packages missing: ", paste(missing_packages, collapse = ", "))
  }
  if (show_density && !requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required for show_density = TRUE. ",
         "Install it with: install.packages('patchwork')")
  }

  # ── Data loading ─────────────────────────────────────────────────────────────
  if (!is.null(data)) {
    df <- as.data.frame(data)
  } else {
    if (!file.exists(input_file)) stop("Input file not found: ", input_file)
    tryCatch(
      { df <- readr::read_delim(input_file, show_col_types = FALSE) },
      error = function(e) stop("Failed to read input file: ", e$message)
    )
  }

  # ── Column validation ─────────────────────────────────────────────────────────
  required_cols <- c(chr_col, pos_col, val_col)
  missing_cols <- required_cols[!required_cols %in% colnames(df)]
  if (length(missing_cols) > 0) {
    stop("Required columns missing from data: ", paste(missing_cols, collapse = ", "))
  }

  # ── Process data ──────────────────────────────────────────────────────────────
  df_processed <- df %>%
    dplyr::select(
      chr   = !!rlang::sym(chr_col),
      pos   = !!rlang::sym(pos_col),
      value = !!rlang::sym(val_col)
    ) %>%
    dplyr::mutate(
      pos   = as.numeric(pos),
      value = as.numeric(value)
    ) %>%
    dplyr::filter(!is.na(pos), !is.na(value), pos > 0, value > 0)

  if (nrow(df_processed) == 0) stop("No valid data remaining after filtering")

  # ── Log10 transform ───────────────────────────────────────────────────────────
  if (transform_log10) {
    df_processed <- df_processed %>% dplyr::mutate(value = -log10(value))
    if (is.null(ylab)) ylab <- expression(-log[10]("p-value"))
  } else {
    if (is.null(ylab)) ylab <- val_col
  }

  # ── Cumulative chromosome positions ──────────────────────────────────────────
  data_cum <- df_processed %>%
    dplyr::group_by(chr) %>%
    dplyr::summarise(max_pos = max(pos, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(chr = factor(chr, levels = gtools::mixedsort(unique(chr)))) %>%
    dplyr::arrange(chr) %>%
    dplyr::mutate(
      cum_pos_end   = cumsum(max_pos),
      cum_pos_start = dplyr::lag(cum_pos_end, default = 0)
    ) %>%
    dplyr::select(chr, cum_pos_start) %>%
    dplyr::left_join(df_processed, ., by = "chr") %>%
    dplyr::arrange(chr, pos) %>%
    dplyr::mutate(
      pos_cum    = pos + cum_pos_start,
      chr_factor = factor(chr, levels = gtools::mixedsort(unique(chr)))
    )

  chromosome_centers <- data_cum %>%
    dplyr::group_by(chr_factor) %>%
    dplyr::summarise(
      center = (min(pos_cum, na.rm = TRUE) + max(pos_cum, na.rm = TRUE)) / 2,
      .groups = "drop"
    ) %>%
    dplyr::arrange(chr_factor)

  # ── Main Manhattan plot ───────────────────────────────────────────────────────
  p_manhattan <- data_cum %>%
    ggplot2::ggplot(ggplot2::aes(x = pos_cum, y = value)) +
    ggplot2::geom_point(
      ggplot2::aes(color = chr_factor),
      alpha = point_alpha,
      size  = point_size
    ) +
    ggplot2::scale_color_manual(
      values = rep(colors, length.out = nlevels(data_cum$chr_factor))
    ) +
    ggplot2::scale_x_continuous(
      breaks = chromosome_centers$center,
      labels = chromosome_centers$chr_factor,
      expand = ggplot2::expansion(mult = c(0.01, 0.01))
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.1))
    ) +
    ggplot2::labs(x = NULL, y = ylab, title = title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position        = "none",
      panel.grid.major.x     = ggplot2::element_blank(),
      panel.grid.minor.x     = ggplot2::element_blank(),
      panel.grid.minor.y     = ggplot2::element_blank(),
      axis.text.x            = ggplot2::element_blank(),   # labels moved to density strip
      axis.ticks.x           = ggplot2::element_blank(),
      plot.title             = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.margin            = ggplot2::margin(5, 5, 0, 5)
    )

  if (!is.null(threshold_line)) {
    p_manhattan <- p_manhattan +
      ggplot2::geom_hline(
        yintercept = threshold_line,
        color      = threshold_color,
        linetype   = "dashed",
        linewidth  = 0.8
      )
  }

  # ── Density panel (optional) ──────────────────────────────────────────────────
  p_density <- NULL
  p_combined <- NULL

  if (show_density) {

    # 1. Compute per-bin SNP counts using cumulative positions
    density_df <- data_cum %>%
      dplyr::mutate(bin = floor(pos_cum / bin_size) * bin_size) %>%
      dplyr::group_by(chr_factor, bin) %>%
      dplyr::summarise(n_snp = dplyr::n(), .groups = "drop")

    # 2. x-axis limits shared with the main plot
    x_range <- range(data_cum$pos_cum, na.rm = TRUE)

    # 3. Build density plot — use n_snp as fill aesthetic so ggplot manages the
    #    gradient scale and can render a proper colour bar legend
    p_density <- ggplot2::ggplot(density_df,
        ggplot2::aes(xmin = bin, xmax = bin + bin_size,
                     ymin = 0,   ymax = 1,
                     fill = n_snp)) +
      ggplot2::geom_rect(color = NA) +
      ggplot2::scale_fill_gradientn(
        colours = density_colors,
        name    = "SNP\nDensity"
      ) +
      ggplot2::scale_x_continuous(
        limits = x_range,
        breaks = chromosome_centers$center,
        labels = chromosome_centers$chr_factor,
        expand = ggplot2::expansion(mult = c(0.01, 0.01))
      ) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::labs(x = "Chromosome", y = NULL) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.grid        = ggplot2::element_blank(),
        axis.text.y       = ggplot2::element_blank(),
        axis.ticks.y      = ggplot2::element_blank(),
        axis.text.x       = ggplot2::element_text(angle = 45, hjust = 1, size = 9),
        plot.margin       = ggplot2::margin(0, 5, 5, 5),
        legend.position   = if (show_density_legend) "right" else "none",
        legend.key.height = ggplot2::unit(0.8, "cm"),
        legend.key.width  = ggplot2::unit(0.3, "cm"),
        legend.title      = ggplot2::element_text(size = 8),
        legend.text       = ggplot2::element_text(size = 7)
      )

    # 4. Assemble with patchwork
    p_combined <- patchwork::wrap_plots(
      p_manhattan, p_density,
      ncol    = 1,
      heights = c(1, density_height)
    )
  }

  # ── Return ────────────────────────────────────────────────────────────────────
  return(list(
    plot.manhattan     = p_manhattan,
    plot.density       = p_density,
    plot.combined      = p_combined,
    data.processed     = data_cum,
    chromosome.centers = chromosome_centers
  ))
}