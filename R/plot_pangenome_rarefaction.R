#' Pan-genome Rarefaction Curve with Heaps' Law Fit
#'
#' @description
#' Fit Heaps' law power models to pan-genome and core-genome rarefaction data
#' and visualize the fitted curves with error bands and data points.
#'
#' @param data A data frame with rarefaction curve data containing sample sizes,
#'   pan-genome/core-genome means and standard deviations.
#' @param sample_col Character string specifying the column name for sample
#'   size (x-axis). Default is "Sample_Size".
#' @param pan_mean_col Character string specifying the column name for
#'   pan-genome mean values. Default is "Pan_Mean".
#' @param pan_std_col Character string specifying the column name for
#'   pan-genome standard deviation. Default is "Pan_Std".
#' @param core_mean_col Character string specifying the column name for
#'   core-genome mean values. Default is "Core_Mean".
#' @param core_std_col Character string specifying the column name for
#'   core-genome standard deviation. Default is "Core_Std".
#' @param pan_color Color for pan-genome elements.
#'   Default is "#D85A30" (coral).
#' @param core_color Color for core-genome elements.
#'   Default is "#1D9E75" (teal).
#' @param title Character string for plot title.
#'   Default is "Pan-genome and Core-genome Rarefaction Curve".
#' @param subtitle Character string for plot subtitle. Default is NULL.
#' @param x_lab Character string for x-axis label.
#'   Default is "Number of genomes".
#' @param y_lab Character string for y-axis label.
#'   Default is "Number of gene families".
#' @param ribbon_alpha Numeric value between 0 and 1 for error ribbon
#'   transparency (±1 SD). Default is 0.15. Set to 0 to hide.
#' @param point_size Numeric value for data point size. Default is 3.
#' @param errorbar_width Numeric value for error bar width in data units.
#'   Default is 0.2.
#' @param line_size Numeric value for fitted curve line width. Default is 1.
#' @param n_smooth Integer specifying the number of points for smooth fitted
#'   curves. Default is 200.
#' @param filename Character string for file path to save the plot.
#'   Default is NULL (no saving).
#' @param width Numeric value for plot width in inches when saving.
#'   Default is 8.
#' @param height Numeric value for plot height in inches when saving.
#'   Default is 5.
#'
#' @return A list containing:
#' \describe{
#'   \item{plot}{A ggplot object of the rarefaction curve.}
#'   \item{fit_pan}{The nls fit object for pan-genome.}
#'   \item{fit_core}{The nls fit object for core-genome.}
#'   \item{r2_pan}{R-squared for pan-genome fit.}
#'   \item{r2_core}{R-squared for core-genome fit.}
#'   \item{heaps_pan}{Named vector of Heaps' law parameters (a, gamma)
#'     for pan-genome.}
#'   \item{heaps_core}{Named vector of Heaps' law parameters (a, alpha)
#'     for core-genome.}
#'   \item{smooth_data}{Data frame of smooth fitted values.}
#' }
#'
#' @details
#' Fits Heaps' law power models using nonlinear least squares:
#' \itemize{
#'   \item Pan-genome: \code{y = a * x^gamma} (gamma > 0)
#'   \item Core-genome: \code{y = a * x^alpha} (alpha < 0)
#' }
#'
#' The Heaps' exponent gamma indicates pan-genome openness:
#' \itemize{
#'   \item gamma close to 0: Closed pan-genome (near saturation)
#'   \item 0 < gamma < 1: Open pan-genome
#'   \item gamma close to 1: Highly open pan-genome
#' }
#'
#' Starting values for the nonlinear fit are auto-estimated via log-log
#' linear regression. Uses \code{minpack.lm::nlsLM()} for robust convergence.
#'
#' @note
#' Reference: Tettelin et al. (2008) Comparative genomics: the bacterial
#' pan-genome. Current Opinion in Microbiology, 11(5), 472-477.
#'
#' @export
#'
#' @examples
#' set.seed(42)
#' n <- 20
#' df <- data.frame(
#'   Sample_Size = 1:n,
#'   Pan_Mean = 4000 * (1:n)^0.05 + rnorm(n, 0, 50),
#'   Pan_Std = abs(rnorm(n, 100, 30)),
#'   Core_Mean = 4500 * (1:n)^(-0.1) + rnorm(n, 0, 30),
#'   Core_Std = abs(rnorm(n, 80, 20))
#' )
#'
#' result <- plot_pangenome_rarefaction(df)
#' result$plot
#' result$r2_pan
#' result$r2_core
#'
plot_pangenome_rarefaction <- function(
    data,
    sample_col = "Sample_Size",
    pan_mean_col = "Pan_Mean",
    pan_std_col = "Pan_Std",
    core_mean_col = "Core_Mean",
    core_std_col = "Core_Std",
    pan_color = "#D85A30",
    core_color = "#1D9E75",
    title = "Pan-genome and Core-genome Rarefaction Curve",
    subtitle = NULL,
    x_lab = "Number of genomes",
    y_lab = "Number of gene families",
    ribbon_alpha = 0.15,
    point_size = 3,
    errorbar_width = 0.2,
    line_size = 1,
    n_smooth = 200,
    filename = NULL,
    width = 8,
    height = 5) {

  # ── Input validation ────────────────────────────────────────────────────────
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }

  required_cols <- c(sample_col, pan_mean_col, pan_std_col,
                     core_mean_col, core_std_col)
  missing_cols <- required_cols[!required_cols %in% colnames(data)]
  if (length(missing_cols) > 0) {
    stop("Required columns missing from data: ",
         paste(missing_cols, collapse = ", "))
  }

  if (nrow(data) < 3) {
    stop("Need at least 3 data points for nonlinear fitting")
  }

  # ── Extract columns ─────────────────────────────────────────────────────────
  x <- data[[sample_col]]
  pan_y <- data[[pan_mean_col]]
  pan_sd <- data[[pan_std_col]]
  core_y <- data[[core_mean_col]]
  core_sd <- data[[core_std_col]]

  # ── Estimate start values via log-log linear regression ──────────────────────
  est_start <- function(x, y) {
    fit <- stats::lm(log(y) ~ log(x))
    coefs <- stats::coef(fit)
    list(a = unname(exp(coefs[1])), b = unname(coefs[2]))
  }

  start_pan <- est_start(x, pan_y)
  start_core <- est_start(x, core_y)

  # ── Fit Heaps' law: y ~ a * x^b ─────────────────────────────────────────────
  fit_df <- data.frame(x = x, pan_y = pan_y, core_y = core_y)

  fit_pan <- minpack.lm::nlsLM(
    pan_y ~ a * x^b,
    data = fit_df,
    start = start_pan
  )

  fit_core <- minpack.lm::nlsLM(
    core_y ~ a * x^b,
    data = fit_df,
    start = start_core
  )

  # ── R-squared ───────────────────────────────────────────────────────────────
  calc_r2 <- function(actual, fitted_obj) {
    ss_res <- sum(stats::residuals(fitted_obj)^2)
    ss_tot <- sum((actual - mean(actual))^2)
    1 - ss_res / ss_tot
  }

  r2_pan <- calc_r2(pan_y, fit_pan)
  r2_core <- calc_r2(core_y, fit_core)

  # ── Heaps' law parameters ───────────────────────────────────────────────────
  heaps_pan <- stats::coef(fit_pan)
  heaps_core <- stats::coef(fit_core)
  names(heaps_pan) <- c("a", "gamma")
  names(heaps_core) <- c("a", "alpha")

  # ── Smooth predictions ──────────────────────────────────────────────────────
  smooth_x <- seq(min(x), max(x), length.out = n_smooth)
  smooth_df <- data.frame(
    x    = smooth_x,
    pan  = stats::predict(fit_pan, newdata = data.frame(x = smooth_x)),
    core = stats::predict(fit_core, newdata = data.frame(x = smooth_x))
  )

  # ── Build plot data ─────────────────────────────────────────────────────────
  plot_df <- data.frame(
    x        = x,
    pan_y    = pan_y,
    pan_ymin = pan_y - pan_sd,
    pan_ymax = pan_y + pan_sd,
    core_y   = core_y,
    core_ymin = core_y - core_sd,
    core_ymax = core_y + core_sd
  )

  # Label positions: right edge, aligned with last fitted values
  label_x <- max(x)
  pan_label_y <- utils::tail(pan_y, 1)
  core_label_y <- utils::tail(core_y, 1)

  # ── Build ggplot ────────────────────────────────────────────────────────────
  p <- ggplot2::ggplot() +

    # Error ribbons (±1 SD)
    ggplot2::geom_ribbon(
      data = plot_df,
      ggplot2::aes(x = .data$x, ymin = .data$pan_ymin, ymax = .data$pan_ymax),
      fill = pan_color, alpha = ribbon_alpha
    ) +
    ggplot2::geom_ribbon(
      data = plot_df,
      ggplot2::aes(x = .data$x, ymin = .data$core_ymin, ymax = .data$core_ymax),
      fill = core_color, alpha = ribbon_alpha
    ) +

    # Fitted curves
    ggplot2::geom_line(
      data = smooth_df,
      ggplot2::aes(x = .data$x, y = .data$pan),
      color = pan_color, linewidth = line_size
    ) +
    ggplot2::geom_line(
      data = smooth_df,
      ggplot2::aes(x = .data$x, y = .data$core),
      color = core_color, linewidth = line_size
    ) +

    # Data points
    ggplot2::geom_point(
      data = plot_df,
      ggplot2::aes(x = .data$x, y = .data$pan_y),
      color = pan_color, size = point_size, shape = 16
    ) +
    ggplot2::geom_point(
      data = plot_df,
      ggplot2::aes(x = .data$x, y = .data$core_y),
      color = core_color, size = point_size, shape = 16
    ) +

    # Error bars
    ggplot2::geom_errorbar(
      data = plot_df,
      ggplot2::aes(x = .data$x, ymin = .data$pan_ymin, ymax = .data$pan_ymax),
      color = pan_color, width = errorbar_width, linewidth = 0.5
    ) +
    ggplot2::geom_errorbar(
      data = plot_df,
      ggplot2::aes(x = .data$x, ymin = .data$core_ymin, ymax = .data$core_ymax),
      color = core_color, width = errorbar_width, linewidth = 0.5
    ) +

    # Curve labels
    ggplot2::annotate("text",
      x = label_x, y = pan_label_y,
      label = "Pan-genome",
      color = pan_color, fontface = "bold", hjust = 1, size = 3.8
    ) +
    ggplot2::annotate("text",
      x = label_x, y = core_label_y,
      label = "Core genome",
      color = core_color, fontface = "bold", hjust = 1, size = 3.8
    ) +

    # Axes
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.02, 0.05))
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::comma,
      expand = ggplot2::expansion(mult = 0.02)
    ) +
    ggplot2::labs(
      title    = title,
      subtitle = subtitle,
      x        = x_lab,
      y        = y_lab
    ) +

    # Theme
    theme_prism() +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(color = "grey50", size = 10),
      panel.grid.major.y = ggplot2::element_line(
        color = "grey92", linewidth = 0.3
      )
    )

  # ── Save if requested ───────────────────────────────────────────────────────
  if (!is.null(filename)) {
    ggplot2::ggsave(filename, plot = p, width = width, height = height)
  }

  # ── Return ──────────────────────────────────────────────────────────────────
  list(
    plot        = p,
    fit_pan     = fit_pan,
    fit_core    = fit_core,
    r2_pan      = r2_pan,
    r2_core     = r2_core,
    heaps_pan   = heaps_pan,
    heaps_core  = heaps_core,
    smooth_data = smooth_df
  )
}
