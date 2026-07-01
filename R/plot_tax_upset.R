#' UpSet Plot of Shared and Unique Taxa
#'
#' @description
#' Create an UpSet plot showing how many features (OTUs/ASVs/taxa) are
#' shared or unique across groups from an abundance table. Supports two
#' plotting engines: \pkg{ggupset} (default) and \pkg{ggVennDiagram}.
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
#'   Only used when \code{engine = "ggupset"}.
#' @param fill Bar fill color. Default is "#2874A6".
#'   Only used when \code{engine = "ggupset"}.
#' @param show_counts Logical. Show count labels on top of bars.
#'   Default is TRUE. Only used when \code{engine = "ggupset"}.
#' @param engine Plotting engine: "ggupset" (default) or "ggVennDiagram".
#' @param relative_height Height ratio between upper and lower panels (default: 2).
#'   Only used when \code{engine = "ggVennDiagram"}.
#' @param relative_width Width ratio for the upset layout (default: 0.3).
#'   Only used when \code{engine = "ggVennDiagram"}.
#' @param verbose Logical. Print progress messages. Default is TRUE.
#'
#' @return A list containing:
#' \describe{
#'   \item{plot}{A ggplot object with the UpSet plot.}
#'   \item{data_upset}{A data frame with the feature-group membership data
#'     used for plotting.}
#'   \item{data_pav}{A presence/absence matrix (feature × group).}
#' }
#'
#' @details
#' The workflow:
#' \enumerate{
#'   \item Merge abundance, sample metadata, and (optionally) taxonomy
#'   \item Convert to presence/absence
#'   \item For each feature, record which groups it appears in
#'   \item Create UpSet plot via selected engine
#' }
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' \dontrun{
#' result <- plot_tax_upset(df.asv, df.sample, df.tax,
#'   which = "phylum", group_col = "treatment"
#' )
#' result$plot
#'
#' # Using ggVennDiagram engine
#' result <- plot_tax_upset(df.asv, df.sample,
#'   group_col = "treatment", engine = "ggVennDiagram"
#' )
#' }
plot_tax_upset <- function(data,
                            sample,
                            taxo = NULL,
                            group_col = "treatment",
                            which = NULL,
                            n_intersections = 20,
                            order_by = "freq",
                            fill = "#2874A6",
                            show_counts = TRUE,
                            engine = "ggupset",
                            relative_height = 2,
                            relative_width = 0.3,
                            verbose = TRUE) {

  if (!engine %in% c("ggupset", "ggVennDiagram")) {
    stop("'engine' must be 'ggupset' or 'ggVennDiagram'")
  }

  if (engine == "ggupset" && !requireNamespace("ggupset", quietly = TRUE)) {
    stop("Package 'ggupset' is required. Install with: install.packages('ggupset')")
  }
  if (engine == "ggVennDiagram" && !requireNamespace("ggVennDiagram", quietly = TRUE)) {
    stop("Package 'ggVennDiagram' is required. Install with: install.packages('ggVennDiagram')")
  }

  # -- Input validation -------------------------------------------------------
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

  # -- Auto-detect and transpose data orientation ------------------------------
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

  # -- Join sample metadata ---------------------------------------------------
  sample_id_col <- colnames(sample)[1]
  sample_select <- sample %>%
    dplyr::select(dplyr::all_of(c(sample_id_col, group_col)))
  colnames(sample_select)[1] <- "sample"
  df_long <- df_long %>%
    dplyr::left_join(sample_select, by = "sample")

  # -- Join taxonomy if provided ----------------------------------------------
  if (!is.null(taxo) && !is.null(which)) {
    taxo_df <- as.data.frame(taxo)
    taxo_col <- colnames(taxo_df)[1]
    taxo_df <- taxo_df %>%
      dplyr::select(dplyr::all_of(c(taxo_col, which))) %>%
      dplyr::rename(.feature_id = !!rlang::sym(taxo_col))

    df_long <- df_long %>%
      dplyr::left_join(taxo_df, by = ".feature_id")

    df_long[[which]] <- ifelse(is.na(df_long[[which]]) |
                                 df_long[[which]] == "", "Unknown",
                               df_long[[which]])
    group_var <- which
  } else {
    group_var <- "feature"
  }

  # -- Convert to presence/absence and aggregate by group ---------------------
  df_pa <- df_long %>%
    dplyr::filter(abundance > 0) %>%
    dplyr::select(group = !!rlang::sym(group_col), feature = !!rlang::sym(group_var)) %>%
    dplyr::distinct()

  # -- Create list column: which groups each feature belongs to ---------------
  df_upset <- df_pa %>%
    dplyr::group_by(feature) %>%
    dplyr::summarize(groups = list(sort(unique(group))), .groups = "drop")

  if (verbose) {
    message("Features: ", nrow(df_upset),
            " | Groups: ", length(unique(df_pa$group)))
  }

  # -- Build UpSet plot -------------------------------------------------------
  all_groups <- sort(unique(df_pa$group))

  if (engine == "ggupset") {
    p <- df_upset %>%
      ggplot2::ggplot(ggplot2::aes(x = groups)) +
      ggplot2::geom_bar(fill = fill, width = 0.6) +
      ggupset::scale_x_upset(
        n_intersections = n_intersections,
        order_by = order_by
      ) +
      ggplot2::labs(y = "Number of features")

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
  } else {
    groups_list <- split(df_pa$feature, df_pa$group)
    groups_list <- groups_list[all_groups]

    venn <- ggVennDiagram::Venn(groups_list)
    p <- ggVennDiagram::plot_upset(
      venn,
      nintersections = n_intersections,
      relative_height = relative_height,
      relative_width = relative_width
    )
  }

  # -- PAV table --------------------------------------------------------------
  df_pav <- df_pa %>%
    dplyr::mutate(present = 1L) %>%
    tidyr::pivot_wider(
      names_from = group, values_from = present,
      values_fill = 0L
    ) %>%
    dplyr::select(feature, dplyr::all_of(all_groups))

  list(
    plot = p,
    data_upset = df_upset,
    data_pav = df_pav
  )
}
