#' Plot Linkage Disequilibrium (LD) Decay with Nonlinear Fitting
#'
#' Bins LD data by distance, fits an exponential decay model per population,
#' calculates half-decay distances, and produces a publication-ready LD decay
#' plot with optional half-decay reference lines.
#'
#' @param data Data frame with at least three columns: population identifier,
#'   pairwise distance, and an LD statistic (e.g., r²).
#' @param pop_col Column name for population grouping (default: "Population").
#' @param dist_col Column name for pairwise physical distance in bp
#'   (default: "Dist").
#' @param value_col Column name for LD statistic (default: "Mean_r2").
#' @param bin_size Numeric. Binning interval in bp (default: 1000).
#'   Distances are grouped into bins of this size and averaged.
#' @param max_dist Numeric or NULL. Maximum distance (in kb) to include in
#'   fitting and plotting (default: NULL, no limit).
#' @param model Character. Decay model to fit. Options:
#'   \itemize{
#'     \item "exponential": \eqn{y ~ a + b * exp(-c * x)} (default)
#'     \item "hill_weir": Hill & Weir (1988) expectation (requires \code{n_ind})
#'   }
#' @param n_ind Integer. Number of individuals, required when
#'   \code{model = "hill_weir"}. Ignored for exponential model.
#' @param show_half_decay Logical. Whether to draw half-decay reference lines
#'   on the plot (default: TRUE).
#' @param palette Character or named vector. Either a viridis palette name
#'   (e.g., "turbo", "viridis") or a named vector of colors keyed by
#'   population (default: "turbo").
#' @param verbose Logical. Print progress messages (default: TRUE).
#'
#' @return A named list containing:
#'   \describe{
#'     \item{\code{plot}}{ggplot2 object of the LD decay curve(s).}
#'     \item{\code{binned_data}}{Data frame of binned mean r² per population.}
#'     \item{\code{fitted_data}}{Data frame of model predictions for smooth
#'       curve plotting.}
#'     \item{\code{half_decay}}{Data frame with columns \code{Population} and
#'       \code{half_decay_kb} (distance at which r² drops to half its
#'       initial value).}
#'     \item{\code{models}}{Named list of fitted model objects per population.}
#'   }
#'
#' @details
#' \strong{Exponential decay model:}
#' \deqn{r^2 = a + b \cdot \exp(-c \cdot d)}
#' where \eqn{a} is the asymptotic background LD, \eqn{b} is the amplitude,
#' \eqn{c} is the decay rate, and \eqn{d} is distance in kb.
#'
#' \strong{Half-decay distance} is defined as the distance (kb) at which
#' the fitted r² drops to half of its initial value (at distance = 0).
#'
#' @references
#' Hill, W. G. & Weir, B. S. (1988). Variances and covariances of squared
#' linkage disequilibria in finite populations. \emph{Theoretical Population
#' Biology}, 33(1), 54-78.
#'
#' @seealso \code{\link{ld_decay_threshold}} for threshold-based LD analysis
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' \dontrun{
#' library(bioRtools)
#'
#' df.ld <- bioRtools::read_data("ld_result.tsv") %>%
#'   dplyr::filter(Population != "Unknown")
#'
#' res <- bioRtools::plot_ld_decay(df.ld, bin_size = 1000)
#'
#' # View plot
#' print(res$plot)
#'
#' # Half-decay distances
#' print(res$half_decay)
#'
#' # Custom palette
#' res <- bioRtools::plot_ld_decay(df.ld, palette = c(
#'   PopA = "#E64B35", PopB = "#4DBBD5", PopC = "#00A087"
#' ))
#' }
plot_ld_decay <- function(data, pop_col = "Population", dist_col = "Dist",
                          value_col = "Mean_r2", bin_size = 1000,
                          max_dist = NULL, model = "exponential",
                          n_ind = NULL, show_half_decay = TRUE,
                          palette = "turbo", verbose = TRUE) {

  # --- Input validation ---
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  for (col in c(pop_col, dist_col, value_col)) {
    if (!col %in% names(data)) {
      stop("Column '", col, "' not found in data")
    }
  }

  model <- match.arg(model, c("exponential", "hill_weir"))
  if (model == "hill_weir" && is.null(n_ind)) {
    stop("'n_ind' is required when model = 'hill_weir'")
  }

  # --- Prepare data ---
  df <- data %>%
    dplyr::select(
      Population = !!rlang::sym(pop_col),
      Dist = !!rlang::sym(dist_col),
      Mean_r2 = !!rlang::sym(value_col)
    ) %>%
    dplyr::filter(!is.na(Dist), !is.na(Mean_r2), Mean_r2 > 0, Dist >= 0) %>%
    dplyr::mutate(Dist_kb = Dist / 1000)

  if (!is.null(max_dist)) {
    df <- df %>% dplyr::filter(Dist_kb <= max_dist)
  }

  populations <- unique(df$Population)
  if (length(populations) == 0) {
    stop("No valid data remaining after filtering")
  }

  # --- Bin data ---
  df_binned <- df %>%
    dplyr::mutate(Dist_bin_kb = floor(Dist_kb / (bin_size / 1000)) * (bin_size / 1000)) %>%
    dplyr::group_by(Population, Dist_bin_kb) %>%
    dplyr::summarise(Mean_r2_bin = mean(Mean_r2, na.rm = TRUE), .groups = "drop")

  if (verbose) {
    message("Binned data: ", nrow(df_binned), " bins across ", length(populations), " population(s)")
  }

  # --- Fit models per population ---
  models_list <- list()
  fitted_dfs <- list()
  half_decay_rows <- list()

  for (pop in populations) {
    d <- df_binned %>% dplyr::filter(Population == pop)

    if (nrow(d) < 4) {
      if (verbose) message("  Skipping '", pop, "': too few bins (", nrow(d), ")")
      next
    }

    # Initial values for exponential model
    y_max <- max(d$Mean_r2_bin)
    y_min <- min(d$Mean_r2_bin)
    start_a <- max(0, y_min * 0.5)
    start_b <- y_max - start_a
    start_c <- 0.01

    fit <- tryCatch({
      if (model == "exponential") {
        minpack.lm::nlsLM(
          Mean_r2_bin ~ a + b * exp(-c * Dist_bin_kb),
          data = d,
          start = list(a = start_a, b = start_b, c = start_c),
          control = list(maxiter = 500)
        )
      } else {
        n <- n_ind
        minpack.lm::nlsLM(
          Mean_r2_bin ~ ((10 + rho * Dist_bin_kb) /
            (22 + 13 * rho * Dist_bin_kb + rho^2 * Dist_bin_kb^2)) *
            (1 + ((3 + rho * Dist_bin_kb) /
              (n * (22 + 13 * rho * Dist_bin_kb + rho^2 * Dist_bin_kb^2)))),
          data = d,
          start = list(rho = 0.1),
          control = list(maxiter = 500)
        )
      }
    }, error = function(e) {
      if (verbose) message("  Model fitting failed for '", pop, "': ", e$message)
      NULL
    })

    if (is.null(fit)) next

    models_list[[pop]] <- fit

    # Generate prediction curve
    pred_dist <- seq(0, max(d$Dist_bin_kb), length.out = 500)
    pred_r2 <- stats::predict(fit, newdata = data.frame(Dist_bin_kb = pred_dist))
    fitted_dfs[[pop]] <- data.frame(
      Population = pop,
      Dist_bin_kb = pred_dist,
      Mean_r2_fitted = pred_r2
    )

    # Calculate half-decay distance
    r2_initial <- pred_r2[1]
    target <- r2_initial / 2
    idx <- which(pred_r2 <= target)[1]
    if (!is.na(idx)) {
      half_decay_rows[[pop]] <- data.frame(
        Population = pop,
        half_decay_kb = round(pred_dist[idx], 2)
      )
    } else {
      half_decay_rows[[pop]] <- data.frame(
        Population = pop,
        half_decay_kb = NA_real_
      )
    }

    if (verbose) {
      hd <- half_decay_rows[[pop]]$half_decay_kb
      message("  ", pop, ": fitted (half-decay = ",
              if (is.na(hd)) "NA" else paste0(hd, " kb"), ")")
    }
  }

  if (length(models_list) == 0) {
    stop("All model fits failed. Check data quality and bin sizes.")
  }

  fitted_df <- do.call(rbind, fitted_dfs)
  rownames(fitted_df) <- NULL
  half_decay_df <- do.call(rbind, half_decay_rows)
  rownames(half_decay_df) <- NULL

  # --- Build plot ---
  p <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = df_binned,
      ggplot2::aes(x = Dist_bin_kb, y = Mean_r2_bin, color = Population),
      linewidth = 0.8, alpha = 0.5
    ) +
    ggplot2::geom_line(
      data = fitted_df,
      ggplot2::aes(x = Dist_bin_kb, y = Mean_r2_fitted, color = Population),
      linewidth = 1
    ) +
    ggplot2::scale_x_continuous(expand = c(0.01, 0)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::labs(
      x = "Distance (kb)",
      y = expression(LD ~ decay ~ (r^2)),
      color = NULL
    ) +
    ggplot2::theme(
      legend.position = "top",
      legend.key.width = ggplot2::unit(1.2, "cm"),
      panel.grid.major.y = ggplot2::element_line(color = "grey90", linewidth = 0.3),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(color = "black")
    )

  # Color palette
  fitted_pops <- names(models_list)
  if (is.character(palette) && length(palette) == 1 && palette %in% c("turbo", "viridis", "plasma", "inferno", "magma", "cividis")) {
    p <- p + ggplot2::scale_color_viridis_d(option = palette)
  } else if (is.named(palette)) {
    p <- p + ggplot2::scale_color_manual(values = palette)
  }

  # Half-decay reference lines
  if (show_half_decay && nrow(half_decay_df) > 0) {
    hd_lines <- half_decay_df[!is.na(half_decay_df$half_decay_kb), ]
    if (nrow(hd_lines) > 0) {
      pop_colors <- ggplot2::ggplot_build(p)$data[[1]]
      # Use dashed lines with matching colors
      for (i in seq_len(nrow(hd_lines))) {
        p <- p + ggplot2::geom_vline(
          xintercept = hd_lines$half_decay_kb[i],
          linetype = "dashed", alpha = 0.4, linewidth = 0.5,
          color = "grey40"
        )
      }
    }
  }

  # --- Return ---
  list(
    plot = p,
    binned_data = df_binned,
    fitted_data = fitted_df,
    half_decay = half_decay_df,
    models = models_list
  )
}
