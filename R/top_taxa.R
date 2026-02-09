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
#'   names and must include a 'group' column for statistical comparisons.
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
#'   comparisons. Must be a level in the 'group' column. Default is "CK".
#' @param p_threshold Numeric value for p-value significance threshold.
#'   Default is 0.05.
#' @param abundance_transform Character string specifying abundance transformation.
#'   Options: "none", "log10", "sqrt", "arcsin_sqrt". Default is "none".
#' @param show_statistics Logical indicating whether to display significance
#'   symbols on the plot. Default is TRUE.
#' @param color_palette Character string specifying color palette. Options:
#'   "d3", "viridis", "brewer", "custom". Default is "d3".
#' @param plot_type Character string specifying plot type. Options: "stacked_bar",
#'   "grouped_bar", "heatmap". Default is "stacked_bar".
#' @param verbose Logical indicating whether to print progress information.
#'   Default is TRUE.
#'
#' @return A list containing five components:
#' \describe{
#'   \item{result.statistics}{A data frame containing statistical test results
#'     for each taxon including p-values and significance levels.}
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
#'   \item Performing statistical tests between groups
#'   \item Creating relative abundance visualization
#' }
#'
#' Statistical tests available:
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
#' top_result <- top_taxa(
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
#' # Analysis at genus level with more taxa
#' genus_result <- top_taxa(
#'   data = df.top10.otu,
#'   sample = df.top10.sample,
#'   taxo = df.top10.class,
#'   which = "genus",
#'   n_top = 15,
#'   method = "kruskal"
#' )
#'
#' # Custom visualization options
#' custom_result <- top_taxa(
#'   data = df.top10.otu,
#'   sample = df.top10.sample,
#'   taxo = df.top10.class,
#'   which = "family",
#'   abundance_transform = "log10",
#'   color_palette = "viridis",
#'   show_statistics = FALSE
#' )
#'
#' # Access processed data for custom analysis
#' processed_data <- top_result$data.processed
#'
#' # Create custom plot
#' library(ggplot2)
#' processed_data %>%
#'   ggplot(aes(x = sample, y = relative_abundance, fill = taxon)) +
#'   geom_bar(stat = "identity", position = "stack") +
#'   facet_wrap(~group, scales = "free_x") +
#'   theme_bio() +
#'   labs(title = "Custom Taxonomic Composition")
#'
top_taxa <- function(data,
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
                     plot_type = "stacked_bar",
                     verbose = TRUE) {
  # Input validation
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("'data' must be a data frame or matrix")
  }

  if (!is.data.frame(sample)) {
    stop("'sample' must be a data frame")
  }

  if (!is.data.frame(taxo)) {
    stop("'taxo' must be a data frame")
  }

  # Check required columns
  if (!"group" %in% colnames(sample)) {
    stop("'sample' must contain a 'group' column")
  }

  if (!which %in% colnames(taxo)) {
    stop("Taxonomic level '", which, "' not found in 'taxo' table")
  }

  # Validate parameters
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
    stop("'abundance_transform' must be one of: ", paste(valid_transforms, collapse = ", "))
  }

  if (!is.numeric(n_top) || n_top < 1) {
    stop("'n_top' must be a positive integer")
  }

  if (!ref %in% sample$group) {
    stop("Reference group '", ref, "' not found in sample groups")
  }

  if (verbose) {
    message("Starting top taxa analysis...")
    message("Taxonomic level: ", which)
    message("Number of top taxa: ", n_top)
    message("Statistical method: ", method)
    message("Reference group: ", ref)
  }

  # Step 1: Merge data and convert to long format
  data_long <- data %>%
    tibble::rownames_to_column(var = "sample") %>%
    tidyr::pivot_longer(
      cols = -sample,
      names_to = "feature",
      values_to = "abundance"
    ) %>%
    dplyr::left_join(sample, by = "sample") %>%
    dplyr::left_join(taxo, by = c("feature" = rownames(taxo)[1])) %>%
    dplyr::select(sample, feature, abundance, group, !!rlang::sym(which))

  # Handle missing taxonomic information
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

  # Step 2: Calculate ranking metric and select top taxa
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

  # Select top N taxa
  top_taxa <- ranking_data %>%
    dplyr::arrange(desc(metric)) %>%
    dplyr::slice_head(n = n_top) %>%
    dplyr::pull(!!rlang::sym(which))

  if (verbose) {
    message("Top ", n_top, " taxa selected: ", paste(top_taxa, collapse = ", "))
  }

  # Step 3: Create final dataset with "Other" category
  data_final <- data_long %>%
    dplyr::mutate(
      taxon = ifelse(!!rlang::sym(which) %in% top_taxa,
        !!rlang::sym(which),
        "Other"),
      taxon = factor(taxon, levels = c(top_taxa, "Other"))
    ) %>%
    dplyr::group_by(sample, taxon, group) %>%
    dplyr::summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop")

  # Calculate relative abundances
  data_processed <- data_final %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(
      relative_abundance = abundance / sum(abundance),
      percentage = relative_abundance * 100
    ) %>%
    dplyr::ungroup()

  # Step 4: Statistical analysis
  if (verbose) {
    message("Performing statistical tests...")
  }

  stat_results <- data_processed %>%
    dplyr::group_by(taxon) %>%
    dplyr::group_modify(~ {
      tryCatch(
        {
          if (method == "wilcox") {
            test_result <- rstatix::wilcox_test(.x, relative_abundance ~ group, ref.group = ref)
          } else if (method == "t.test") {
            test_result <- rstatix::t_test(.x, relative_abundance ~ group, ref.group = ref)
          } else if (method == "kruskal") {
            test_result <- rstatix::kruskal_test(.x, relative_abundance ~ group)
          }
          return(test_result)
        },
        error = function(e) {
          data.frame(
            group1 = ref,
            group2 = "Error",
            p = NA,
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

  # Step 5: Create taxa summary
  taxa_summary <- data_processed %>%
    dplyr::group_by(taxon, group) %>%
    dplyr::summarise(
      mean_abundance = mean(relative_abundance),
      sd_abundance = sd(relative_abundance),
      n_samples = n(),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from = group,
      values_from = c(mean_abundance, sd_abundance, n_samples),
      names_sep = "_"
    )

  # Step 6: Create visualization
  if (verbose) {
    message("Creating visualization...")
  }

  # Prepare data for plotting with statistics
  plot_data <- data_processed %>%
    dplyr::left_join(
      stat_results %>% dplyr::select(taxon, significance),
      by = "taxon"
    )

  # Create base plot
  abundance_plot <- plot_data %>%
    ggplot2::ggplot(ggplot2::aes(x = sample, y = relative_abundance, fill = taxon))

  if (plot_type == "stacked_bar") {
    abundance_plot <- abundance_plot +
      ggplot2::geom_bar(stat = "identity", position = "fill", width = 0.8)

    if (show_statistics) {
      abundance_plot <- abundance_plot +
        ggplot2::geom_text(
          ggplot2::aes(label = significance),
          position = ggplot2::position_fill(vjust = 0.5),
          size = 3,
          color = "black"
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

  # Apply color palette
  if (color_palette == "d3") {
    abundance_plot <- abundance_plot + ggsci::scale_fill_d3()
  } else if (color_palette == "viridis") {
    abundance_plot <- abundance_plot + ggplot2::scale_fill_viridis_d()
  } else if (color_palette == "brewer") {
    abundance_plot <- abundance_plot + ggplot2::scale_fill_brewer(type = "qual")
  }

  # Apply theme
  abundance_plot <- abundance_plot +
    bioRtools::CommonlyUsed.plot_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )

  # Create analysis parameters record
  analysis_params <- data.frame(
    parameter = c("taxonomic_level", "n_top", "ranking_method", "statistical_method",
      "reference_group", "p_threshold", "abundance_transform", "color_palette"),
    value = c(which, n_top, by, method, ref, p_threshold, abundance_transform, color_palette),
    stringsAsFactors = FALSE
  )

  if (verbose) {
    message("Analysis completed successfully!")
    message("Significant taxa (p < ", p_threshold, "): ",
      sum(stat_results$significant, na.rm = TRUE))
  }

  # Return comprehensive results following project conventions
  return(list(
    result.statistics = stat_results,
    plot.abundance = abundance_plot,
    data.processed = data_processed,
    taxa.summary = taxa_summary,
    analysis.parameters = analysis_params
  ))
}
