#' Mantel Test Correlogram
#'
#' @description
#' Create a Mantel test correlogram combining environmental correlation
#' heatmap with Mantel test results. For each sample group (defined by
#' \code{group_col}), characteristic features are selected and tested
#' against environmental variables via \pkg{linkET}.
#'
#' @param data_spec Feature abundance table. First column is feature ID
#'   (e.g., "asv"), remaining columns are samples (features × samples).
#'   Auto-transposed internally.
#' @param data_env Environmental data table. Must contain a \code{sample}
#'   column and environmental variable columns (either in long format with
#'   \code{sample}, \code{indicator}, \code{value}, or in wide format).
#'   Auto-detected.
#' @param data_sample Sample metadata table. Must contain a \code{sample}
#'   column matching samples in \code{data_spec} and \code{data_env},
#'   and a group column.
#' @param group_col Column name in \code{data_sample} defining sample groups.
#'   Default is \code{"group"}. Falls back to auto-detection if not found.
#' @param spec_select_method How to select characteristic features per group:
#'   \code{"prevalence"} (features present above threshold, default) or
#'   \code{"abundance"} (top N by mean abundance).
#' @param prevalence_threshold Prevalence threshold for feature selection
#'   when \code{spec_select_method = "prevalence"}. Default is 0.3.
#' @param top_n Number of top features to select per group when
#'   \code{spec_select_method = "abundance"}. Default is 50.
#' @param spec_dist Distance method for feature data. Default is
#'   \code{"bray"}.
#' @param env_dist Distance method for environmental data. Default is
#'   \code{"euclidean"}.
#' @param cor_method Correlation method for environmental variables.
#'   Default is \code{"pearson"}.
#' @param type Plot type passed to \code{linkET::qcorrplot()}. Default is
#'   \code{"lower"}.
#' @param diag Logical. Show diagonal in correlation plot. Default is
#'   \code{FALSE}.
#' @param r_breaks Numeric vector. Break points for Mantel's r categories.
#'   Default is \code{c(-Inf, 0.2, 0.4, Inf)}.
#' @param r_labels Character vector. Labels for r categories.
#'   Default is \code{c("< 0.2", "0.2 - 0.4", ">= 0.4")}.
#' @param p_breaks Numeric vector. Break points for Mantel's p categories.
#'   Default is \code{c(-Inf, 0.01, 0.05, Inf)}.
#' @param p_labels Character vector. Labels for p categories.
#'   Default is \code{c("< 0.01", "0.01 - 0.05", ">= 0.05")}.
#' @param fill_colors Vector of colors for the correlation fill scale.
#'   Default uses \code{RColorBrewer::brewer.pal(11, "RdBu")}.
#' @param size_values Numeric vector. Size values for Mantel's r categories.
#'   Default is \code{c(0.5, 1, 2)}.
#' @param p_colors Vector of colors for p-value categories. Default uses
#'   \code{linkET::color_pal(3)}.
#' @param verbose Logical. Print progress messages. Default is \code{TRUE}.
#'
#' @return A list containing:
#' \describe{
#'   \item{plot}{A ggplot object with the Mantel correlogram.}
#'   \item{mantel}{The Mantel test result data frame.}
#'   \item{spec_select}{The feature selection list used.}
#' }
#'
#' @details
#' The function accepts flexible input formats:
#' \itemize{
#'   \item \code{data_env}: either long format (\code{sample}, \code{indicator},
#'         \code{value}) or wide format (\code{sample} + indicators as columns).
#'         Auto-detected.
#'   \item \code{data_spec}: features × samples, first column is feature ID.
#'         Auto-transposed to samples × features.
#'   \item \code{data_sample}: must have \code{sample} and group columns.
#'         \code{group_col} is auto-detected if \code{NULL}.
#' }
#'
#' The workflow:
#' \enumerate{
#'   \item Match samples across all three input tables
#'   \item For each sample group, select characteristic features
#'         (prevalence-based or abundance-based)
#'   \item Create \code{spec_select} list mapping group names to feature
#'         column indices
#'   \item Run \code{linkET::mantel_test()} to compute Mantel correlations
#'         between each feature group and each environmental variable
#'   \item Run \code{linkET::correlate()} on environmental data
#'   \item Combine with \code{linkET::qcorrplot()} and
#'         \code{linkET::geom_couple()}
#' }
#'
#' Requires \pkg{linkET} and \pkg{RColorBrewer}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' df.phy <- readxl::read_excel("理化性质.xlsx") %>%
#'   dplyr::select(sample, indicator, value) %>%
#'   tidyr::pivot_wider(names_from = indicator, values_from = value)
#'
#' df.asv <- readxl::read_excel("ASV丰度表.xlsx") %>%
#'   dplyr::select(asv, df.phy$sample)
#'
#' df.sample <- readxl::read_excel("样品信息表.xlsx") %>%
#'   dplyr::filter(sample %in% df.phy$sample)
#'
#' result <- plot_mantel(df.asv, df.phy, df.sample)
#' result$plot
#' }
#'
plot_mantel <- function(data_spec,
                        data_env,
                        data_sample,
                        group_col = "group",
                        spec_select_method = c("prevalence", "abundance"),
                        prevalence_threshold = 0.3,
                        top_n = 50,
                        spec_dist = "bray",
                        env_dist = "euclidean",
                        cor_method = "pearson",
                        type = "lower",
                        diag = FALSE,
                        r_breaks = c(-Inf, 0.2, 0.4, Inf),
                        r_labels = c("< 0.2", "0.2 - 0.4", ">= 0.4"),
                        p_breaks = c(-Inf, 0.01, 0.05, Inf),
                        p_labels = c("< 0.01", "0.01 - 0.05", ">= 0.05"),
                        fill_colors = NULL,
                        size_values = c(0.5, 1, 2),
                        p_colors = NULL,
                        verbose = TRUE) {

  if (!requireNamespace("linkET", quietly = TRUE)) {
    stop("Package 'linkET' is required. Install with: install.packages('linkET')")
  }
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("Package 'RColorBrewer' is required. Install with: install.packages('RColorBrewer')")
  }

  spec_select_method <- match.arg(spec_select_method)

  # ── Input validation ──────────────────────────────────────────────────────
  if (!is.data.frame(data_spec))   stop("'data_spec' must be a data frame")
  if (!is.data.frame(data_env))    stop("'data_env' must be a data frame")
  if (!is.data.frame(data_sample)) stop("'data_sample' must be a data frame")

  if (!"sample" %in% colnames(data_sample)) {
    stop("'data_sample' must contain a 'sample' column")
  }

  # ── Auto-detect group column ──────────────────────────────────────────────
  if (!group_col %in% colnames(data_sample)) {
    if (verbose) message("'", group_col, "' not found, auto-detecting group column...")
    cand <- setdiff(colnames(data_sample), "sample")
    if (length(cand) == 0) stop("No group column found in 'data_sample'")
    n_uniq <- vapply(data_sample[cand], function(x) length(unique(x)), integer(1))
    cand <- cand[n_uniq > 1]
    if (length(cand) == 0) stop("No group column with >1 unique values found in 'data_sample'")
    n_uniq <- vapply(data_sample[cand], function(x) length(unique(x)), integer(1))
    group_col <- cand[which.min(n_uniq)]
    if (verbose) message("Auto-detected group column: '", group_col, "'")
  }

  if (!group_col %in% colnames(data_sample)) {
    stop("'", group_col, "' column not found in 'data_sample'")
  }

  # ── Parse data_env: auto-detect long vs wide format ────────────────────────
  data_env <- as.data.frame(data_env)
  has_indicator_value <- all(c("indicator", "value") %in% colnames(data_env)) &&
    "sample" %in% colnames(data_env)

  if (has_indicator_value) {
    if (verbose) message("Detected data_env in long format, pivoting...")
    data_env <- data_env %>%
      dplyr::select(sample, indicator, value) %>%
      tidyr::pivot_wider(names_from = indicator, values_from = value)
  }

  if (!"sample" %in% colnames(data_env)) {
    stop("'data_env' must contain a 'sample' column")
  }

  # ── Parse data_spec: features × samples → samples × features ──────────────
  data_spec <- as.data.frame(data_spec)
  feature_col <- colnames(data_spec)[1]
  feature_names <- as.character(data_spec[[feature_col]])

  if (verbose) message("Transposing data_spec: features × samples → samples × features")

  data_spec <- data_spec %>%
    dplyr::select(-dplyr::all_of(feature_col)) %>%
    t() %>%
    as.data.frame()
  colnames(data_spec) <- feature_names
  data_spec$sample <- rownames(data_spec)

  # ── Harmonize samples ─────────────────────────────────────────────────────
  common_samples <- Reduce(intersect, list(
    data_spec$sample,
    data_env$sample,
    as.character(data_sample$sample)
  ))

  if (length(common_samples) == 0) {
    stop("No common samples found across input tables")
  }

  if (verbose) message("Common samples: ", length(common_samples))

  data_spec <- data_spec[match(common_samples, data_spec$sample), , drop = FALSE]
  data_env  <- data_env[match(common_samples, data_env$sample), , drop = FALSE]
  data_sample <- data_sample[match(common_samples, data_sample$sample), , drop = FALSE]

  # Set sample as rowname for downstream use
  rownames(data_spec) <- data_spec$sample
  data_spec$sample <- NULL
  rownames(data_env) <- data_env$sample
  data_env$sample <- NULL

  # Keep only numeric columns in env
  env_numeric <- vapply(data_env, is.numeric, logical(1))
  if (!all(env_numeric)) {
    data_env <- data_env[, env_numeric, drop = FALSE]
    if (verbose) message("Dropped non-numeric columns from 'data_env'")
  }

  # ── Build spec_select from sample groups ───────────────────────────────────
  groups_vec <- as.character(data_sample[[group_col]])
  unique_groups <- unique(groups_vec)

  if (verbose) {
    message("Sample groups: ", paste(unique_groups, collapse = ", "))
  }

  spec_select <- list()

  for (grp in unique_groups) {
    idx <- which(groups_vec == grp)
    grp_spec <- data_spec[idx, , drop = FALSE]

    if (spec_select_method == "prevalence") {
      prevalence <- colSums(grp_spec > 0) / nrow(grp_spec)
      selected <- which(prevalence >= prevalence_threshold)
    } else {
      mean_abund <- colMeans(grp_spec)
      selected <- order(mean_abund, decreasing = TRUE)[seq_len(min(top_n, ncol(grp_spec)))]
    }

    if (length(selected) == 0) {
      warning("No features selected for group '", grp, "'. Skipping.")
      next
    }

    spec_select[[grp]] <- selected

    if (verbose) {
      message("  Group '", grp, "': ", length(selected), " features selected")
    }
  }

  if (length(spec_select) == 0) {
    stop("No features selected for any group")
  }

  # ── Run Mantel test ────────────────────────────────────────────────────────
  if (verbose) message("Running Mantel test...")

  mantel <- linkET::mantel_test(
    data_spec,
    data_env,
    spec_select = spec_select,
    spec_dist = spec_dist,
    env_dist = env_dist
  )

  mantel <- mantel %>%
    dplyr::mutate(
      rd = cut(r, breaks = r_breaks, labels = r_labels),
      pd = cut(p, breaks = p_breaks, labels = p_labels)
    )

  if (verbose) {
    message("Mantel test: ", nrow(mantel), " comparisons")
  }

  # ── Correlate environmental variables ──────────────────────────────────────
  if (verbose) message("Computing environmental correlations...")

  cor_env <- linkET::correlate(data_env, method = cor_method)

  # ── Build plot ────────────────────────────────────────────────────────────
  if (is.null(fill_colors)) {
    fill_colors <- RColorBrewer::brewer.pal(11, "RdBu")
  }
  if (is.null(p_colors)) {
    p_colors <- linkET::color_pal(length(p_labels))
  }

  p <- linkET::qcorrplot(cor_env, type = type, diag = diag) +
    linkET::geom_square() +
    linkET::geom_couple(
      ggplot2::aes(colour = pd, size = rd),
      data = mantel,
      curvature = linkET::nice_curvature()
    ) +
    ggplot2::scale_fill_gradientn(colours = fill_colors) +
    ggplot2::scale_size_manual(values = size_values) +
    ggplot2::scale_colour_manual(values = p_colors) +
    ggplot2::guides(
      size = ggplot2::guide_legend(
        title = "Mantel's r",
        override.aes = list(colour = "grey35"),
        order = 2
      ),
      colour = ggplot2::guide_legend(
        title = "Mantel's p",
        override.aes = list(size = 3),
        order = 1
      ),
      fill = ggplot2::guide_colorbar(
        title = "Pearson's r",
        order = 3
      )
    )

  list(
    plot        = p,
    mantel      = mantel,
    spec_select = spec_select
  )
}
