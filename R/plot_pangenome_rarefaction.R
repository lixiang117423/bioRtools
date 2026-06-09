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
#' @author Xiang LI <lixiang117423@gmail.com>
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

  # ── Fit Heaps' law or constant model ────────────────────────────────────────
  fit_df <- data.frame(x = x, pan_y = pan_y, core_y = core_y)
  pan_var <- stats::var(pan_y)
  core_var <- stats::var(core_y)

  # Pan-genome
  if (is.na(pan_var) || pan_var < .Machine$double.eps) {
    fit_pan <- NULL
    r2_pan <- NA_real_
    heaps_pan <- c(a = pan_y[1L], gamma = 0)
    pan_smooth_val <- rep(pan_y[1L], n_smooth)
  } else {
    start_pan <- list(a = pan_y[which.min(x)], b = 0.1)
    fit_pan <- stats::nls(pan_y ~ a * x^b, data = fit_df, start = start_pan)
    ss_res <- sum(stats::residuals(fit_pan)^2)
    ss_tot <- sum((pan_y - mean(pan_y))^2)
    r2_pan <- 1 - ss_res / ss_tot
    heaps_pan <- stats::coef(fit_pan)
    names(heaps_pan) <- c("a", "gamma")
    pan_smooth_val <- NULL
  }

  # Core-genome
  if (is.na(core_var) || core_var < .Machine$double.eps) {
    fit_core <- NULL
    r2_core <- NA_real_
    heaps_core <- c(a = core_y[1L], alpha = 0)
    core_smooth_val <- rep(core_y[1L], n_smooth)
  } else {
    start_core <- list(a = core_y[which.min(x)], b = -0.1)
    fit_core <- stats::nls(core_y ~ a * x^b, data = fit_df, start = start_core)
    ss_res <- sum(stats::residuals(fit_core)^2)
    ss_tot <- sum((core_y - mean(core_y))^2)
    r2_core <- 1 - ss_res / ss_tot
    heaps_core <- stats::coef(fit_core)
    names(heaps_core) <- c("a", "alpha")
    core_smooth_val <- NULL
  }

  # ── Smooth predictions ──────────────────────────────────────────────────────
  smooth_x <- seq(min(x), max(x), length.out = n_smooth)
  smooth_df <- data.frame(
    x    = smooth_x,
    pan  = if (!is.null(fit_pan))
             stats::predict(fit_pan, newdata = data.frame(x = smooth_x))
           else pan_smooth_val,
    core = if (!is.null(fit_core))
             stats::predict(fit_core, newdata = data.frame(x = smooth_x))
           else core_smooth_val
  )

  # ── Print fitted parameters ─────────────────────────────────────────────────
  if (!is.na(r2_pan)) {
    message(sprintf("Pan-genome  | a = %.4f, gamma = %.4f, R2 = %.4f",
                    heaps_pan["a"], heaps_pan["gamma"], r2_pan))
  } else {
    message("Pan-genome  | constant data (gamma = 0)")
  }
  if (!is.na(r2_core)) {
    message(sprintf("Core-genome | a = %.4f, alpha = %.4f, R2 = %.4f",
                    heaps_core["a"], heaps_core["alpha"], r2_core))
  } else {
    message("Core-genome | constant data (alpha = 0)")
  }

  # ── Build plot data (long format for legend) ────────────────────────────────
  n_pts <- length(x)
  plot_long <- data.frame(
    x    = rep(x, 2),
    y    = c(pan_y, core_y),
    ymin = c(pan_y - pan_sd, core_y - core_sd),
    ymax = c(pan_y + pan_sd, core_y + core_sd),
    type = rep(c("Pan", "Core"), each = n_pts),
    stringsAsFactors = FALSE
  )
  plot_long$type <- factor(plot_long$type,
                           levels = c("Pan", "Core"))

  smooth_long <- data.frame(
    x    = rep(smooth_x, 2),
    y    = c(smooth_df$pan, smooth_df$core),
    type = rep(c("Pan", "Core"), each = n_smooth),
    stringsAsFactors = FALSE
  )
  smooth_long$type <- factor(smooth_long$type,
                             levels = c("Pan", "Core"))

  color_map <- c("Pan" = pan_color, "Core" = core_color)

  # ── Build ggplot ────────────────────────────────────────────────────────────
  p <- ggplot2::ggplot() +

    # Error ribbons (±1 SD)
    ggplot2::geom_ribbon(
      data = plot_long,
      ggplot2::aes(x = .data$x, ymin = .data$ymin, ymax = .data$ymax,
                   fill = .data$type),
      alpha = ribbon_alpha
    ) +

    # Fitted curves
    ggplot2::geom_line(
      data = smooth_long,
      ggplot2::aes(x = .data$x, y = .data$y, color = .data$type),
      linewidth = line_size
    ) +

    # Data points
    ggplot2::geom_point(
      data = plot_long,
      ggplot2::aes(x = .data$x, y = .data$y, color = .data$type),
      size = point_size, shape = 16
    ) +

    # Error bars
    ggplot2::geom_errorbar(
      data = plot_long,
      ggplot2::aes(x = .data$x, ymin = .data$ymin, ymax = .data$ymax,
                   color = .data$type),
      width = errorbar_width, linewidth = 0.5
    ) +

    # Color / fill scales with legend
    ggplot2::scale_color_manual(values = color_map, name = NULL) +
    ggplot2::scale_fill_manual(values = color_map, guide = "none") +

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
    ggprism::theme_prism() +
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
