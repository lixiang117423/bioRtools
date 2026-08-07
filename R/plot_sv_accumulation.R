#' Pan- and Core-SV Accumulation Curve
#'
#' @description
#' Compute pan- and core-SV accumulation curves from a structural variant
#' presence/absence matrix via random permutation of sample order, and
#' visualize the curves with ±1 SD ribbons.
#'
#' @param data A data.frame, data.table, or a file path (character string).
#'   If a file path, it is read with \code{data.table::fread()} using
#'   tab-separated, header, and check.names = FALSE.
#' @param meta_cols Character vector of column names treated as SV metadata
#'   (not sample columns). All remaining columns are treated as sample
#'   presence/absence columns. Default identifies SV_ID, SVTYPE, SVLEN,
#'   CHROM, START, END.
#' @param n_perm Integer. Number of random permutations of sample order.
#'   Larger values produce smoother curves. Default is 100.
#' @param core_threshold Numeric between 0 and 1. An SV is considered "core"
#'   when present in at least this proportion of the first k samples.
#'   Default is 0.95 (i.e., ≥95% of samples). Set to 1 for strict "all
#'   samples" core definition.
#' @param seed Integer. Random seed for reproducible permutations.
#'   Default is 1234.
#' @param title Character. Plot title.
#'   Default is "Pan- and Core-SV accumulation curve".
#' @param x_lab Character. X-axis label. Default is "Number of genomes".
#' @param y_lab Character. Y-axis label. Default is "Number of SVs".
#' @param pan_color Character. Color for pan-SV elements.
#'   Default is "#D85A30".
#' @param core_color Character. Color for core-SV elements.
#'   Default is "#1D9E75".
#' @param ribbon_alpha Numeric. Transparency for ±1 SD ribbons, between 0
#'   and 1. Default is 0.2. Set to 0 to hide ribbons.
#' @param line_size Numeric. Line width for accumulation curves.
#'   Default is 1.
#' @param point_size Numeric. Point size for data markers.
#'   Default is 1.
#' @param point_alpha Numeric. Point transparency. Default is 0.6.
#'
#' @return A list containing:
#' \describe{
#'   \item{plot}{A ggplot object of the pan/core accumulation curves.}
#'   \item{data}{Data frame with columns n_genome, n_SV, sd_SV, type.}
#'   \item{pan_mat}{Integer matrix of pan-SV counts across permutations
#'     (n_perm rows by n_sample columns).}
#'   \item{core_mat}{Integer matrix of core-SV counts across permutations
#'     (n_perm rows by n_sample columns).}
#'   \item{n_sample}{Number of samples.}
#'   \item{n_sv}{Total number of SVs.}
#' }
#'
#' @details
#' For each permutation, sample order is randomly shuffled and cumulative
#' pan-SV (present in at least 1 sample) and core-SV (present in at least
#' \code{core_threshold} proportion of the first k samples) counts are
#' computed across 1 to all genomes. The mean and standard deviation across
#' all permutations are then plotted as curves with ±1 SD ribbons.
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' \dontrun{
#' # From file
#' result <- plot_sv_accumulation(
#'   "~/data/sv_results.tsv",
#'   n_perm = 200
#' )
#' result$plot
#' head(result$data)
#' }
plot_sv_accumulation <- function(data,
                                 meta_cols = c("SV_ID", "SVTYPE", "SVLEN",
                                               "CHROM", "START", "END"),
                                 n_perm = 100,
                                 core_threshold = 0.95,
                                 seed = 1234,
                                 title = "Pan- and Core-SV accumulation curve",
                                 x_lab = "Number of genomes",
                                 y_lab = "Number of SVs",
                                 pan_color = "#D85A30",
                                 core_color = "#1D9E75",
                                 ribbon_alpha = 0.2,
                                 line_size = 1,
                                 point_size = 1,
                                 point_alpha = 0.6) {

  # ── Read file if path provided ──────────────────────────────────────────────
  if (is.character(data) && length(data) == 1L) {
    if (!file.exists(data)) {
      stop("File not found: ", data)
    }
    data <- data.table::fread(data, sep = "\t", header = TRUE,
                              check.names = FALSE)
  }

  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame, data.table, or a file path (character)")
  }

  # ── Identify sample columns ─────────────────────────────────────────────────
  sample_cols <- setdiff(colnames(data), meta_cols)
  if (length(sample_cols) == 0L) {
    stop("No sample columns found after excluding meta_cols")
  }

  n_sample <- length(sample_cols)
  if (n_sample < 2L) {
    stop("Need at least 2 sample columns, found ", n_sample)
  }

  if (!is.numeric(n_perm) || n_perm < 1L) {
    stop("'n_perm' must be a positive integer, got ", n_perm)
  }

  if (!is.numeric(core_threshold) || core_threshold <= 0 || core_threshold > 1) {
    stop("'core_threshold' must be between 0 and 1 (exclusive of 0)")
  }

  if (!is.numeric(ribbon_alpha) || ribbon_alpha < 0 || ribbon_alpha > 1) {
    stop("'ribbon_alpha' must be between 0 and 1")
  }

  # ── Convert to integer presence/absence matrix ──────────────────────────────
  mat <- as.matrix(data[, sample_cols, drop = FALSE])
  mode(mat) <- "integer"

  # ── Random permutation: Pan/Core accumulation ───────────────────────────────
  set.seed(seed)

  pan_mat <- matrix(NA_integer_, nrow = n_perm, ncol = n_sample)
  core_mat <- matrix(NA_integer_, nrow = n_perm, ncol = n_sample)

  for (i in seq_len(n_perm)) {
    perm_order <- sample(sample_cols, n_sample)
    cum_mat <- mat[, perm_order, drop = FALSE]
    cum_presence <- t(apply(cum_mat, 1, cumsum))

    for (k in seq_len(n_sample)) {
      pan_mat[i, k] <- sum(cum_presence[, k] > 0)
      core_mat[i, k] <- sum(cum_presence[, k] >= ceiling(core_threshold * k))
    }
  }

  # ── Build plot data frame ───────────────────────────────────────────────────
  pan_df <- data.frame(
    n_genome = rep(seq_len(n_sample), 2),
    n_SV     = c(colMeans(pan_mat), colMeans(core_mat)),
    sd_SV    = c(apply(pan_mat, 2, stats::sd), apply(core_mat, 2, stats::sd)),
    type     = rep(c("Pan", "Core"), each = n_sample),
    stringsAsFactors = FALSE
  )

  # ── Build ggplot ────────────────────────────────────────────────────────────
  color_map <- c("Pan" = pan_color, "Core" = core_color)

  p <- ggplot2::ggplot(
    pan_df,
    ggplot2::aes(
      x = .data$n_genome, y = .data$n_SV,
      color = .data$type, fill = .data$type
    )
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$n_SV - .data$sd_SV,
                   ymax = .data$n_SV + .data$sd_SV),
      alpha = ribbon_alpha, color = NA
    ) +
    ggplot2::geom_line(linewidth = line_size) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    ggplot2::scale_color_manual(values = color_map) +
    ggplot2::scale_fill_manual(values = color_map) +
    ggplot2::labs(
      title = title,
      x     = x_lab,
      y     = y_lab,
      color = NULL,
      fill  = NULL
    ) +
    ggprism::theme_prism() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "top"
    )

  # ── Return ──────────────────────────────────────────────────────────────────
  list(
    plot     = p,
    data     = pan_df,
    pan_mat  = pan_mat,
    core_mat = core_mat,
    n_sample = n_sample,
    n_sv     = nrow(mat)
  )
}
