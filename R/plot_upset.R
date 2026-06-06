#' UpSet Plot of Shared and Unique Taxa
#'
#' @description
#' Create an UpSet plot showing how many taxa (OTUs/ASVs) are
#' shared or unique across groups. Uses abundance table, sample metadata,
#' and optional taxonomy table. Uses \pkg{ggupset} to produce
#' publication-quality combination matrix plots.
#'
#' @param data ASV/OTU abundance table. Can be in either orientation:
#'   features × samples (first column is feature ID) or
#'   samples × features (rownames are sample IDs). Auto-detected.
#' @param sample Sample metadata table. Must contain sample names matching
#'   those in data.
#' @param taxo Optional taxonomy table. If provided and \code{which} is
#'   specified, aggregation is done at that taxonomic level instead of
#'   OTU level.
#' @param group_col Column name in the sample table defining groups.
#'   Default is "treatment".
#' @param which Taxonomic level to analyze (e.g., "phylum", "genus").
#'   Must be a column name in \code{taxo}. Default is NULL (OTU level).
#' @param n_intersections Integer, maximum number of intersections to show.
#'   Default is 20.
#' @param order_by How to order intersections: "freq" (by count, default)
#'   or "degree" (by number of groups).
#' @param fill Bar fill color. Default is "#2874A6".
#' @param show_counts Logical. Show count labels on top of bars.
#'   Default is TRUE.
#' @param verbose Logical. Print progress messages. Default is TRUE.
#'
#' @return A list containing:
#' \describe{
#'   \item{plot}{A ggplot object with the UpSet plot.}
#'   \item{data.upset}{A data frame with the feature-group membership data
#'     used for plotting.}
#'   \item{data.pav}{A presence/absence matrix (feature × group).}
#' }
#'
#' @details
#' The workflow:
#' \enumerate{
#'   \item Merge abundance, sample metadata, and (optionally) taxonomy
#'   \item Convert to presence/absence
#'   \item For each feature, record which groups it appears in
#'   \item Create UpSet plot via \code{ggupset::scale_x_upset()}
#' }
#'
#' Requires \pkg{ggupset}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- plot_tax_upset(df.asv, df.sample, df.tax,
#'   which = "phylum", group_col = "treatment"
#' )
#' result$plot
#' }
#'
plot_tax_upset <- function(data,
                       sample,
                       taxo = NULL,
                       group_col = "treatment",
                       which = NULL,
                       n_intersections = 20,
                       order_by = "freq",
                       fill = "#2874A6",
                       show_counts = TRUE,
                       verbose = TRUE) {

  if (!requireNamespace("ggupset", quietly = TRUE)) {
    stop("Package 'ggupset' is required. Install with: install.packages('ggupset')")
  }

  # ── Input validation ──────────────────────────────────────────────────────
  if (!is.data.frame(data)) stop("'data' must be a data frame")
  if (!is.data.frame(sample)) stop("'sample' must be a data frame")
  if (!group_col %in% colnames(sample)) {
    stop("'", group_col, "' column not found in sample table")
  }
  if (!is.null(which) && is.null(taxo)) {
    stop("'taxo' must be provided when 'which' is specified")
  }
  if (!is.null(which) && !which %in% colnames(taxo)) {
    stop("'", which, "' column not found in taxo table")
  }

  # ── Auto-detect and transpose data orientation ────────────────────────────
  data <- as.data.frame(data)
  sample_names <- as.character(sample[[1]])
  data_colnames <- colnames(data)
  data_rownames <- rownames(data)

  cols_match <- if (ncol(data) > 1)
    sum(data_colnames[-1] %in% sample_names) else 0
  rows_match <- if (!is.null(data_rownames))
    sum(data_rownames %in% sample_names) else 0

  if (cols_match > rows_match && cols_match >= 2) {
    if (verbose) message("Detected features x samples format, transposing...")
    feature_names <- as.character(data[[1]])
    feature_col <- data_colnames[1]
    data_matrix <- as.matrix(data[, -1, drop = FALSE])
    storage.mode(data_matrix) <- "numeric"
    data_matrix <- t(data_matrix)
    colnames(data_matrix) <- feature_names
    df_long <- as.data.frame(data_matrix) %>%
      tibble::rownames_to_column(var = "sample") %>%
      tidyr::pivot_longer(-sample, names_to = "feature", values_to = "abundance")
    df_long$.feature_id <- df_long$feature
  } else {
    df_long <- data %>%
      tibble::rownames_to_column(var = "sample") %>%
      tidyr::pivot_longer(-sample, names_to = "feature", values_to = "abundance")
    df_long$.feature_id <- df_long$feature
  }

  # ── Join sample metadata ──────────────────────────────────────────────────
  sample_id_col <- colnames(sample)[1]
  sample_select <- sample %>%
    dplyr::select(dplyr::all_of(c(sample_id_col, group_col)))
  colnames(sample_select)[1] <- "sample"
  df_long <- df_long %>%
    dplyr::left_join(sample_select, by = "sample")

  # ── Join taxonomy if provided ─────────────────────────────────────────────
  if (!is.null(taxo) && !is.null(which)) {
    taxo_df <- as.data.frame(taxo)
    taxo_col <- colnames(taxo_df)[1]
    taxo_df <- taxo_df %>%
      dplyr::select(dplyr::all_of(c(taxo_col, which))) %>%
      dplyr::rename(.feature_id = !!rlang::sym(taxo_col))

    df_long <- df_long %>%
      dplyr::left_join(taxo_df, by = ".feature_id")

    # Handle NA taxonomy
    df_long[[which]] <- ifelse(is.na(df_long[[which]]) |
                                 df_long[[which]] == "", "Unknown",
                               df_long[[which]])
    group_var <- which
  } else {
    group_var <- "feature"
  }

  # ── Convert to presence/absence and aggregate by group ────────────────────
  df_pa <- df_long %>%
    dplyr::filter(abundance > 0) %>%
    dplyr::select(group = !!rlang::sym(group_col), feature = !!rlang::sym(group_var)) %>%
    dplyr::distinct()

  # ── Create list column: which groups each feature belongs to ──────────────
  df_upset <- df_pa %>%
    dplyr::group_by(feature) %>%
    dplyr::summarize(groups = list(sort(unique(group))), .groups = "drop")

  if (verbose) {
    message("Features: ", nrow(df_upset),
            " | Groups: ", length(unique(df_pa$group)))
  }

  # ── Build UpSet plot ──────────────────────────────────────────────────────
  p <- df_upset %>%
    ggplot2::ggplot(ggplot2::aes(x = groups)) +
    ggplot2::geom_bar(fill = fill, width = 0.6) +
    ggupset::scale_x_upset(
      n_intersections = n_intersections,
      order_by = order_by
    ) +
    ggplot2::labs(
      y = "Number of features"
    )

  if (show_counts) {
    p <- p +
      ggplot2::geom_text(
        stat = "count",
        ggplot2::aes(label = ggplot2::after_stat(count)),
        vjust = -0.3, size = 3.5
      )
  }

  p <- p +
    ggupset::theme_combmatrix() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "grey92", linewidth = 0.3)
    )

  # ── PAV table: feature x group presence/absence matrix ────────────────────
  all_groups <- sort(unique(df_pa$group))
  df_pav <- df_pa %>%
    dplyr::mutate(present = 1L) %>%
    tidyr::pivot_wider(
      names_from = group, values_from = present,
      values_fill = 0L
    ) %>%
    dplyr::select(feature, dplyr::all_of(all_groups))

  list(
    plot = p,
    data.upset = df_upset,
    data.pav = df_pav
  )
}


#' UpSet Plot from a Group-Value Data Frame
#'
#' @description
#' Create an UpSet plot from a data frame containing group and value
#' columns (e.g., comparison and feature_id). Shows which values are
#' shared or unique across groups. Uses \pkg{ggupset} to produce
#' publication-quality combination matrix plots.
#'
#' @param data A data frame with at least two columns: one for group
#'   labels and one for value/feature identifiers.
#' @param group Column name in \code{data} containing the group labels.
#'   Default is \code{"comparison"}.
#' @param value Column name in \code{data} containing the value/feature
#'   identifiers. Default is \code{"feature_id"}.
#' @param n_intersections Integer, maximum number of intersections to show.
#'   Default is 20.
#' @param order_by How to order intersections: \code{"freq"} (by count,
#'   default) or \code{"degree"} (by number of groups).
#' @param fill Bar fill color. Default is \code{"#2874A6"}.
#' @param show_counts Logical. Show count labels on top of bars.
#'   Default is \code{TRUE}.
#'
#' @return A list containing:
#' \describe{
#'   \item{plot}{A ggplot object with the UpSet plot.}
#'   \item{data.upset}{A data frame with the value-group membership data
#'     used for plotting.}
#'   \item{data.pav}{A presence/absence matrix (value × group).}
#' }
#'
#' @details
#' The workflow:
#' \enumerate{
#'   \item Select group and value columns from the input data frame
#'   \item Deduplicate to get unique group-value pairs
#'   \item For each value, record which groups it appears in
#'   \item Create UpSet plot via \code{ggupset::scale_x_upset()}
#'   \item Build a presence/absence matrix
#' }
#'
#' Requires \pkg{ggupset}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' df.deg %>%
#'   dplyr::select(comparison, feature_id) %>%
#'   bioRtools::plot_upset()
#' }
#'
plot_upset <- function(data,
                       group = "comparison",
                       value = "feature_id",
                       n_intersections = 20,
                       order_by = "freq",
                       fill = "#2874A6",
                       show_counts = TRUE) {

  if (!requireNamespace("ggupset", quietly = TRUE)) {
    stop("Package 'ggupset' is required. Install with: install.packages('ggupset')")
  }

  if (!is.data.frame(data)) stop("'data' must be a data frame")
  if (!group %in% colnames(data)) {
    stop("'", group, "' column not found in 'data'")
  }
  if (!value %in% colnames(data)) {
    stop("'", value, "' column not found in 'data'")
  }

  df_pa <- data %>%
    dplyr::select(group = !!rlang::sym(group),
                  feature = !!rlang::sym(value)) %>%
    dplyr::distinct()

  df_upset <- df_pa %>%
    dplyr::group_by(feature) %>%
    dplyr::summarize(groups = list(sort(unique(group))), .groups = "drop")

  all_groups <- sort(unique(df_pa$group))
  df_pav <- df_pa %>%
    dplyr::mutate(present = 1L) %>%
    tidyr::pivot_wider(
      names_from = group, values_from = present,
      values_fill = 0L
    ) %>%
    dplyr::select(feature, dplyr::all_of(all_groups))

  p <- df_upset %>%
    ggplot2::ggplot(ggplot2::aes(x = groups)) +
    ggplot2::geom_bar(fill = fill, width = 0.6) +
    ggupset::scale_x_upset(
      n_intersections = n_intersections,
      order_by = order_by
    ) +
    ggplot2::labs(
      y = "Number of features"
    )

  if (show_counts) {
    p <- p +
      ggplot2::geom_text(
        stat = "count",
        ggplot2::aes(label = ggplot2::after_stat(count)),
        vjust = -0.3, size = 3.5
      )
  }

  p <- p +
    ggupset::theme_combmatrix() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "grey92", linewidth = 0.3)
    )

  list(
    plot = p,
    data.upset = df_upset,
    data.pav = df_pav
  )
}
