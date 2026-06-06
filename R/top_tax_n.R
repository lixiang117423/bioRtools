#' Analyze and visualize top abundant taxa at different taxonomic levels
#'
#' @description
#' Identify the most abundant microorganisms at specified taxonomic levels and
#' perform statistical comparisons between groups. This function creates both
#' relative abundance visualizations and statistical analysis results for
#' microbiome community composition studies.
#'
#' @param data Feature abundance table (OTU/ASV table) where rows represent
#'   samples and columns represent features. Row names should be sample IDs
#'   matching those in the sample table.
#' @param sample Sample metadata table where the first column contains sample
#'   names and must include a group column for statistical comparisons.
#' @param taxo Taxonomic classification table where rows represent features
#'   (matching column names in data) and columns represent taxonomic levels
#'   (e.g., phylum, class, order, family, genus, species).
#' @param which Character string specifying which taxonomic level to analyze.
#'   Must be a column name in the taxo table. Default is "phylum".
#' @param n_top Integer specifying the number of top taxa to display separately.
#'   Remaining taxa will be grouped as "Other". Default is 9 (showing top 9 + Other = 10 total).
#' @param by Character string specifying how to calculate abundance rankings.
#'   Options: "sum" (total abundance), "mean" (average abundance), "prevalence"
#'   (number of samples present). Default is "sum".
#' @param method Character string specifying the statistical test method.
#'   Options: "wilcox" (Wilcoxon rank-sum test), "t.test" (t-test), "kruskal"
#'   (Kruskal-Wallis test). Default is "wilcox".
#' @param ref Character string specifying the reference group for pairwise
#'   comparisons. Must be a level in the group column. Default is "CK".
#' @param p_threshold Numeric value for p-value significance threshold.
#'   Default is 0.05.
#' @param abundance_transform Character string specifying abundance transformation.
#'   Options: "none", "log10", "sqrt", "arcsin_sqrt". Default is "none".
#' @param show_statistics Logical indicating whether to display significance
#'   symbols on the plot. Default is TRUE.
#' @param color_palette Character string specifying color palette. Options:
#'   "d3", "viridis", "brewer", "custom". Default is "d3".
#' @param plot_type Character string specifying plot type. Options: "stacked_bar",
#'   "alluvial". Default is "stacked_bar". "alluvial" creates a flow/alluvial
#'   diagram using \pkg{ggalluvial}.
#' @param group_col Character string specifying the column name in the sample
#'   table that defines groups. Default is "group".
#' @param plot_by Character string specifying the aggregation level for plotting.
#'   Options: "sample" (plot each sample separately), "group" (aggregate to
#'   group level and plot by group). Default is "sample".
#' @param show_labels Logical indicating whether to show percentage labels
#'   inside bars. Default is FALSE.
#' @param verbose Logical indicating whether to print progress information.
#'   Default is TRUE.
#'
#' @return A list containing five components:
#' \describe{
#'   \item{result.statistics}{A data frame containing statistical test results
#'     for each taxon including p-values and significance levels. NULL when
#'     \code{plot_by = "group"}.}
#'   \item{plot.abundance}{A ggplot object showing relative abundance
#'     visualization with optional significance annotations.}
#'   \item{data.processed}{A data frame containing the processed abundance
#'     data used for plotting and analysis.}
#'   \item{taxa.summary}{A summary table showing the top taxa selected and
#'     their average abundances across groups.}
#'   \item{analysis.parameters}{A data frame recording all analysis parameters
#'     for reproducibility.}
#' }
#'
#' @details
#' The analysis workflow includes:
#' \enumerate{
#'   \item Merging abundance, sample, and taxonomic data
#'   \item Aggregating abundance by taxonomic level
#'   \item Selecting top N most abundant taxa
#'   \item Grouping remaining taxa as "Other"
#'   \item Performing statistical tests between groups (when \code{plot_by = "sample"})
#'   \item Creating relative abundance visualization
#' }
#'
#' When \code{plot_by = "group"}, the aggregation follows a two-step approach:
#' first averaging across samples within each group for each feature, then
#' summing across features within each taxon-group combination. This produces
#' a group-level composition bar chart.
#'
#' Statistical tests available (only when \code{plot_by = "sample"}):
#' \itemize{
#'   \item Wilcoxon rank-sum test: Non-parametric, robust for microbiome data
#'   \item t-test: Parametric test for normally distributed data
#'   \item Kruskal-Wallis test: Non-parametric for multiple groups
#' }
#'
#' @note
#' \itemize{
#'   \item Sample names must match between data row names and sample table
#'   \item Feature names must match between data column names and taxo row names
#'   \item Missing taxonomic information will be labeled as "Unclassified"
#'   \item Very low abundance taxa may be filtered out automatically
#' }
#'
#' @export
#'
#' @examples
#' library(dplyr)
#' library(bioRtools)
#'
#' # Load example microbiome data
#' data(df.top10.otu)
#' data(df.top10.sample)
#' data(df.top10.class)
#'
#' # Basic top taxa analysis at phylum level
#' top_result <- top_tax_n(
#'   data = df.top10.otu,
#'   sample = df.top10.sample,
#'   taxo = df.top10.class,
#'   which = "phylum"
#' )
#'
#' # View the plot
#' top_result$plot.abundance
#'
#' # Check statistical results
#' top_result$result.statistics
#'
#' # View taxa summary
#' top_result$taxa.summary
#'
#' # Group-level composition with percentage labels
#' group_result <- top_tax_n(
#'   data = df.top10.otu,
#'   sample = df.top10.sample,
#'   taxo = df.top10.class,
#'   which = "phylum",
#'   group_col = "group",
#'   plot_by = "group",
#'   show_labels = TRUE
#' )
#' group_result$plot.abundance
#'
top_tax_n <- function(data,
                     sample,
                     taxo,
                     which = "phylum",
                     n_top = 9,
                     by = "sum",
                     method = "wilcox",
                     ref = "CK",
                     p_threshold = 0.05,
                     abundance_transform = "none",
                     show_statistics = TRUE,
                     color_palette = "d3",
                     plot_type = "alluvial",
                     group_col = "treatment",
                     plot_by = "group",
                     show_labels = FALSE,
                     verbose = TRUE) {

  # ── Input validation ──────────────────────────────────────────────────────
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("'data' must be a data frame or matrix")
  }
  if (!is.data.frame(sample)) {
    stop("'sample' must be a data frame")
  }
  if (!is.data.frame(taxo)) {
    stop("'taxo' must be a data frame")
  }
  if (!group_col %in% colnames(sample)) {
    stop("'", group_col, "' column not found in sample table")
  }
  if (!which %in% colnames(taxo)) {
    stop("Taxonomic level '", which, "' not found in 'taxo' table")
  }

  valid_methods <- c("wilcox", "t.test", "kruskal")
  if (!method %in% valid_methods) {
    stop("'method' must be one of: ", paste(valid_methods, collapse = ", "))
  }
  valid_by <- c("sum", "mean", "prevalence")
  if (!by %in% valid_by) {
    stop("'by' must be one of: ", paste(valid_by, collapse = ", "))
  }
  valid_transforms <- c("none", "log10", "sqrt", "arcsin_sqrt")
  if (!abundance_transform %in% valid_transforms) {
    stop("'abundance_transform' must be one of: ",
         paste(valid_transforms, collapse = ", "))
  }
  valid_plot_by <- c("sample", "group")
  if (!plot_by %in% valid_plot_by) {
    stop("'plot_by' must be one of: ", paste(valid_plot_by, collapse = ", "))
  }
  if (!is.numeric(n_top) || n_top < 1) {
    stop("'n_top' must be a positive integer")
  }

  # ── Normalize group column name ──────────────────────────────────────────
  # Internally always use "group" regardless of the actual column name
  if (group_col != "group") {
    if ("group" %in% colnames(sample)) {
      sample <- sample %>% dplyr::select(-group)
    }
    sample <- sample %>% dplyr::rename(group = !!rlang::sym(group_col))
  }

  # ── Auto-detect and transpose data orientation ───────────────────────────
  sample_names <- sample[[1]]
  data_colnames <- colnames(data)
  data_rownames <- rownames(data)

  # Check if column names (excluding first col) match sample names
  cols_match <- if (ncol(data) > 1)
    sum(data_colnames[-1] %in% sample_names) else 0
  rows_match <- if (!is.null(data_rownames))
    sum(data_rownames %in% sample_names) else 0

  if (cols_match > rows_match && cols_match >= 2) {
    # Data is features × samples: first column is feature IDs, rest are samples
    if (verbose) message("Detected features x samples format, transposing...")
    feature_names <- as.character(data[[1]])
    data_matrix <- as.matrix(data[, -1, drop = FALSE])
    storage.mode(data_matrix) <- "numeric"
    data <- t(data_matrix)
    colnames(data) <- feature_names
    data <- as.data.frame(data)
  } else {
    data <- as.data.frame(data)
  }

  if (!ref %in% sample$group) {
    stop("Reference group '", ref, "' not found in sample groups")
  }

  if (verbose) {
    message("Starting top taxa analysis...")
    message("Taxonomic level: ", which)
    message("Number of top taxa: ", n_top)
    message("Plot by: ", plot_by)
    if (plot_by == "sample") {
      message("Statistical method: ", method)
      message("Reference group: ", ref)
    }
  }

  # ── Normalize taxo table: ensure feature IDs are in a column ──────────────
  taxo <- as.data.frame(taxo)
  feature_names <- colnames(data)
  # If rownames are feature IDs, convert to column
  if (all(feature_names %in% rownames(taxo))) {
    taxo <- taxo %>% tibble::rownames_to_column(var = ".feature_id")
  } else {
    # Use first column as feature ID
    taxo_join_col <- colnames(taxo)[1]
    taxo <- taxo %>% dplyr::rename(.feature_id = !!rlang::sym(taxo_join_col))
  }

  # ── Step 1: Merge data and convert to long format ─────────────────────────
  data_long <- data %>%
    tibble::rownames_to_column(var = "sample") %>%
    tidyr::pivot_longer(
      cols = -sample,
      names_to = "feature",
      values_to = "abundance"
    ) %>%
    dplyr::left_join(sample, by = "sample") %>%
    dplyr::left_join(taxo, by = c("feature" = ".feature_id")) %>%
    dplyr::select(sample, feature, abundance, group, !!rlang::sym(which))

  data_long <- data_long %>%
    dplyr::mutate(
      !!which := ifelse(is.na(!!rlang::sym(which)) | !!rlang::sym(which) == "",
        "Unclassified",
        !!rlang::sym(which))
    )

  # Apply abundance transformation
  if (abundance_transform != "none") {
    data_long <- data_long %>%
      dplyr::mutate(
        abundance = switch(abundance_transform,
          "log10" = log10(abundance + 1),
          "sqrt" = sqrt(abundance),
          "arcsin_sqrt" = asin(sqrt(abundance / max(abundance, na.rm = TRUE))),
          abundance
        )
      )
  }

  # ── Step 2: Calculate ranking metric and select top taxa ──────────────────
  if (by == "sum") {
    ranking_data <- data_long %>%
      dplyr::group_by(!!rlang::sym(which)) %>%
      dplyr::summarise(metric = sum(abundance, na.rm = TRUE), .groups = "drop")
  } else if (by == "mean") {
    ranking_data <- data_long %>%
      dplyr::group_by(!!rlang::sym(which)) %>%
      dplyr::summarise(metric = mean(abundance, na.rm = TRUE), .groups = "drop")
  } else if (by == "prevalence") {
    ranking_data <- data_long %>%
      dplyr::group_by(!!rlang::sym(which)) %>%
      dplyr::summarise(metric = sum(abundance > 0), .groups = "drop")
  }

  top_taxa_vec <- ranking_data %>%
    dplyr::arrange(desc(metric)) %>%
    dplyr::slice_head(n = n_top) %>%
    dplyr::pull(!!rlang::sym(which))

  if (verbose) {
    message("Top ", n_top, " taxa selected: ",
            paste(top_taxa_vec, collapse = ", "))
  }

  # ── Step 3: Create final dataset ──────────────────────────────────────────
  if (plot_by == "group") {
    # Group-level aggregation: mean per sample, then sum by taxon+group
    data_processed <- data_long %>%
      dplyr::mutate(
        taxon = ifelse(!!rlang::sym(which) %in% top_taxa_vec,
                       !!rlang::sym(which), "Other"),
        taxon = factor(taxon, levels = c(top_taxa_vec, "Other"))
      ) %>%
      dplyr::group_by(taxon, feature, group) %>%
      dplyr::summarise(mean.sample = mean(abundance, na.rm = TRUE),
                       .groups = "drop") %>%
      dplyr::group_by(taxon, group) %>%
      dplyr::summarise(abundance = sum(mean.sample, na.rm = TRUE),
                       .groups = "drop") %>%
      dplyr::group_by(group) %>%
      dplyr::mutate(
        relative_abundance = abundance / sum(abundance),
        percentage = relative_abundance * 100
      ) %>%
      dplyr::ungroup()

  } else {
    # Sample-level (original behavior)
    data_processed <- data_long %>%
      dplyr::mutate(
        taxon = ifelse(!!rlang::sym(which) %in% top_taxa_vec,
                       !!rlang::sym(which), "Other"),
        taxon = factor(taxon, levels = c(top_taxa_vec, "Other"))
      ) %>%
      dplyr::group_by(sample, taxon, group) %>%
      dplyr::summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
      dplyr::group_by(sample) %>%
      dplyr::mutate(
        relative_abundance = abundance / sum(abundance),
        percentage = relative_abundance * 100
      ) %>%
      dplyr::ungroup()
  }

  # ── Step 4: Statistical analysis ──────────────────────────────────────────
  if (plot_by == "sample") {
    if (verbose) message("Performing statistical tests...")

    stat_results <- data_processed %>%
      dplyr::group_by(taxon) %>%
      dplyr::group_modify(~ {
        tryCatch(
          {
            if (method == "wilcox") {
              test_result <- rstatix::wilcox_test(.x,
                relative_abundance ~ group, ref.group = ref)
            } else if (method == "t.test") {
              test_result <- rstatix::t_test(.x,
                relative_abundance ~ group, ref.group = ref)
            } else if (method == "kruskal") {
              test_result <- rstatix::kruskal_test(.x,
                relative_abundance ~ group)
            }
            return(test_result)
          },
          error = function(e) {
            data.frame(
              group1 = ref, group2 = "Error", p = NA,
              method = paste("Failed:", method)
            )
          })
      }) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        p_adj = p.adjust(p, method = "fdr"),
        significance = dplyr::case_when(
          p_adj < 0.001 ~ "***",
          p_adj < 0.01 ~ "**",
          p_adj < 0.05 ~ "*",
          p_adj < 0.1 ~ ".",
          TRUE ~ ""
        ),
        significant = p_adj < p_threshold
      )
  } else {
    stat_results <- NULL
    if (verbose) message("Skipping statistics (plot_by = 'group')")
  }

  # ── Step 5: Taxa summary ──────────────────────────────────────────────────
  if (plot_by == "sample") {
    taxa_summary <- data_processed %>%
      dplyr::group_by(taxon, group) %>%
      dplyr::summarise(
        mean_abundance = mean(relative_abundance),
        sd_abundance = sd(relative_abundance),
        n_samples = dplyr::n(),
        .groups = "drop"
      ) %>%
      tidyr::pivot_wider(
        names_from = group,
        values_from = c(mean_abundance, sd_abundance, n_samples),
        names_sep = "_"
      )
  } else {
    taxa_summary <- data_processed %>%
      dplyr::select(taxon, group, relative_abundance) %>%
      tidyr::pivot_wider(
        names_from = group,
        values_from = relative_abundance,
        names_prefix = "prop_"
      )
  }

  # ── Step 6: Visualization ─────────────────────────────────────────────────
  if (verbose) message("Creating visualization...")

  if (plot_by == "group") {

    if (plot_type == "alluvial") {
      # Alluvial / flow diagram
      if (!requireNamespace("ggalluvial", quietly = TRUE)) {
        stop("Package 'ggalluvial' is required for alluvial plots. ",
             "Install with: install.packages('ggalluvial')")
      }
      abundance_plot <- data_processed %>%
        ggplot2::ggplot(ggplot2::aes(
          x = group, y = relative_abundance,
          alluvium = taxon, stratum = taxon
        )) +
        ggalluvial::geom_alluvium(
          ggplot2::aes(fill = taxon), color = NA, width = 0.4
        ) +
        ggalluvial::geom_stratum(
          ggplot2::aes(fill = taxon), color = NA, width = 0.4
        ) +
        ggplot2::scale_y_continuous(
          labels = scales::percent, limits = c(0, 1), expand = c(0, 0)
        ) +
        ggplot2::scale_x_discrete(expand = c(0, 0)) +
        ggplot2::labs(
          x = stringr::str_to_title(group_col),
          y = "Relative abundance (%)",
          fill = stringr::str_to_title(which),
          color = stringr::str_to_title(which)
        )

    } else {
      # Group-level stacked bar
      abundance_plot <- data_processed %>%
        ggplot2::ggplot(ggplot2::aes(
          x = group, y = relative_abundance, fill = taxon
        )) +
        ggplot2::geom_col(width = 0.8) +
        ggplot2::scale_x_discrete(expand = c(0, 0)) +
        ggplot2::scale_y_continuous(
          labels = scales::percent, limits = c(0, 1), expand = c(0, 0)
        ) +
        ggplot2::labs(
          x = stringr::str_to_title(group_col),
          y = "Relative abundance (%)",
          fill = stringr::str_to_title(which)
        )

      if (show_labels) {
        abundance_plot <- abundance_plot +
          ggplot2::geom_text(
            ggplot2::aes(
              label = scales::percent(relative_abundance, accuracy = 1)
            ),
            position = ggplot2::position_stack(vjust = 0.5),
            size = 3
          )
      }
    }

  } else {
    # Sample-level (original behavior)
    plot_data <- data_processed %>%
      dplyr::left_join(
        stat_results %>% dplyr::select(taxon, significance),
        by = "taxon"
      )

    abundance_plot <- plot_data %>%
      ggplot2::ggplot(
        ggplot2::aes(x = sample, y = relative_abundance, fill = taxon)
      )

    if (plot_type == "stacked_bar") {
      abundance_plot <- abundance_plot +
        ggplot2::geom_bar(stat = "identity", position = "fill", width = 0.8)

      if (show_statistics && !is.null(stat_results)) {
        abundance_plot <- abundance_plot +
          ggplot2::geom_text(
            ggplot2::aes(label = significance),
            position = ggplot2::position_fill(vjust = 0.5),
            size = 3, color = "black"
          )
      }

      if (show_labels) {
        abundance_plot <- abundance_plot +
          ggplot2::geom_text(
            ggplot2::aes(
              label = scales::percent(relative_abundance, accuracy = 1)
            ),
            position = ggplot2::position_fill(vjust = 0.5),
            size = 2.5
          )
      }

      abundance_plot <- abundance_plot +
        ggplot2::scale_y_continuous(
          expand = c(0, 0),
          labels = scales::percent_format()
        ) +
        ggplot2::labs(
          x = "Sample",
          y = "Relative Abundance",
          fill = stringr::str_to_title(which),
          title = paste("Top", n_top, stringr::str_to_title(which), "Composition")
        )
    }
  }

  # Apply color palette
  if (color_palette == "d3") {
    abundance_plot <- abundance_plot + ggsci::scale_fill_d3()
  } else if (color_palette == "viridis") {
    abundance_plot <- abundance_plot +
      ggplot2::scale_fill_viridis_d()
  } else if (color_palette == "brewer") {
    abundance_plot <- abundance_plot +
      ggplot2::scale_fill_brewer(type = "qual")
  }

  # Apply theme
  if (plot_by == "group") {
    abundance_plot <- abundance_plot +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 20, hjust = 1, vjust = 1),
        legend.position = "right"
      )
  } else {
    abundance_plot <- abundance_plot +
      ggprism::theme_prism() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        legend.position = "right"
      )
  }

  # ── Analysis parameters ──────────────────────────────────────────────────
  analysis_params <- data.frame(
    parameter = c("taxonomic_level", "n_top", "ranking_method",
                  "statistical_method", "reference_group", "p_threshold",
                  "abundance_transform", "color_palette", "group_col",
                  "plot_by", "show_labels"),
    value = c(which, n_top, by, method, ref, p_threshold,
              abundance_transform, color_palette, group_col,
              plot_by, show_labels),
    stringsAsFactors = FALSE
  )

  if (verbose) {
    message("Analysis completed successfully!")
    if (plot_by == "sample" && !is.null(stat_results)) {
      message("Significant taxa (p < ", p_threshold, "): ",
              sum(stat_results$significant, na.rm = TRUE))
    }
  }

  list(
    result.statistics  = stat_results,
    plot.abundance     = abundance_plot,
    data.processed     = data_processed,
    taxa.summary       = taxa_summary,
    analysis.parameters = analysis_params
  )
}
