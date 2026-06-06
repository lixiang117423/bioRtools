#' Mantel Test Correlogram
#'
#' @description
#' Create a Mantel test correlogram combining environmental correlation
#' heatmap with Mantel test results. For each sample group (defined by
#' \code{group_col}), characteristic features are selected and tested
#' against environmental variables via \pkg{linkET}.
#'
#' @param data.spec Feature abundance table. Samples as rows, features
#'   (ASVs/OTUs) as columns.
#' @param data.env Environmental data table. Samples as rows, environmental
#'   variables as columns. Must contain the same samples as \code{data.spec}.
#' @param data.sample Sample metadata table. Must contain a sample ID column
#'   (first column) matching rownames of \code{data.spec} and \code{data.env},
#'   and a group column.
#' @param group_col Column name in \code{data.sample} defining sample groups.
#'   Default is \code{"group"}.
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
#'   tidyr::pivot_wider(names_from = indicator, values_from = value) %>%
#'   tibble::column_to_rownames(var = "sample")
#'
#' df.asv <- readxl::read_excel("ASV丰度表.xlsx") %>%
#'   dplyr::select(asv, rownames(df.phy)) %>%
#'   tibble::column_to_rownames(var = "asv") %>%
#'   t() %>% as.data.frame()
#'
#' df.sample <- readxl::read_excel("样品信息表.xlsx") %>%
#'   dplyr::filter(sample %in% rownames(df.asv)) %>%
#'   dplyr::select(sample, alt) %>%
#'   dplyr::rename(group = alt)
#'
#' result <- plot_mantel(df.asv, df.phy, df.sample, group_col = "group")
#' result$plot
#' }
#'
plot_mantel <- function(data.spec,
                        data.env,
                        data.sample,
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
  if (!is.data.frame(data.spec)) stop("'data.spec' must be a data frame")
  if (!is.data.frame(data.env))  stop("'data.env' must be a data frame")
  if (!is.data.frame(data.sample)) stop("'data.sample' must be a data frame")
  if (!group_col %in% colnames(data.sample)) {
    stop("'", group_col, "' column not found in 'data.sample'")
  }

  # ── Harmonize sample IDs ──────────────────────────────────────────────────
  data.spec  <- as.data.frame(data.spec)
  data.env   <- as.data.frame(data.env)
  sample_id_col <- colnames(data.sample)[1]

  common_samples <- Reduce(intersect, list(
    rownames(data.spec),
    rownames(data.env),
    as.character(data.sample[[sample_id_col]])
  ))

  if (length(common_samples) == 0) {
    stop("No common samples found across 'data.spec', 'data.env' and 'data.sample'")
  }

  if (verbose) {
    message("Common samples: ", length(common_samples))
  }

  data.spec  <- data.spec[common_samples, , drop = FALSE]
  data.env   <- data.env[common_samples, , drop = FALSE]
  data.sample <- data.sample[match(common_samples, as.character(data.sample[[sample_id_col]])), , drop = FALSE]

  # Keep only numeric columns in env
  env_numeric <- vapply(data.env, is.numeric, logical(1))
  if (!all(env_numeric)) {
    data.env <- data.env[, env_numeric, drop = FALSE]
    if (verbose) message("Dropped non-numeric columns from 'data.env'")
  }

  # ── Build spec_select from sample groups ───────────────────────────────────
  groups_vec <- as.character(data.sample[[group_col]])
  unique_groups <- unique(groups_vec)

  if (verbose) {
    message("Sample groups: ", paste(unique_groups, collapse = ", "))
  }

  spec_select <- list()

  for (grp in unique_groups) {
    idx <- which(groups_vec == grp)
    grp_spec <- data.spec[idx, , drop = FALSE]

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
    data.spec,
    data.env,
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

  cor_env <- linkET::correlate(data.env, method = cor_method)

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
