#' Random Forest Taxonomic Classification
#'
#' @description
#' Perform random forest classification at multiple taxonomic levels using
#' \pkg{mlr3}. Aggregates ASV/OTU abundance at each taxonomic level, trains
#' a tuned random forest model, and returns accuracy, predictions, and
#' variable importance.
#'
#' @param data ASV/OTU abundance table (features × samples, first column is
#'   feature ID).
#' @param taxo Taxonomic classification table (first column is feature ID
#'   matching data).
#' @param sample Sample metadata table. Must contain a column matching sample
#'   names in data and a grouping column.
#' @param tax_cols Named character vector mapping English level names to actual
#'   column names in the taxo table. Default maps standard English names
#'   (phylum, class, order, family, genus, species).
#' @param feature_col Column name for feature IDs in data and taxo.
#'   Default is the first column name.
#' @param sample_col Column name for sample IDs in the sample table.
#'   Default is "sample".
#' @param group_col Column name for the grouping/response variable in the
#'   sample table.
#' @param tax_levels Character vector specifying which taxonomic levels to
#'   analyze. Use "OTU" for feature-level. Default analyzes all levels.
#' @param ratio Train/test split ratio. Default is 0.7.
#' @param seed Random seed for data splitting. Default is 1288.
#' @param seed_tune Random seed for tuning. Default is 8812.
#' @param num_trees_range Integer vector of length 2: range for number of
#'   trees search. Default is c(1, 25) (internally multiplied by 20).
#' @param min_node_size_range Integer vector of length 2: range for minimum
#'   node size. Default is c(3, 30).
#' @param term_evals Number of tuning evaluations. Default is 10.
#' @param cv_folds Number of cross-validation folds for tuning. Default is 10.
#' @param top_n Number of top important features to display in boxplot.
#'   Default is 9.
#' @param fill_palette Fill palette for the top features boxplot.
#'   Default uses \code{scale_fill_research()}.
#' @param verbose Logical. Print progress messages. Default is TRUE.
#'
#' @return A list containing:
#' \describe{
#'   \item{accuracy}{Data frame of prediction accuracy at each taxonomic level.}
#'   \item{predictions}{Named list of mlr3 prediction objects per level.}
#'   \item{importance}{Data frame of variable importance across all levels.}
#'   \item{plot_accuracy}{ggplot bar chart of accuracy by taxonomic level.}
#'   \item{plot_top_features}{ggplot boxplot of top N features at OTU level.}
#' }
#'
#' @details
#' Requires packages: \pkg{mlr3}, \pkg{mlr3learners}, \pkg{mlr3tuning},
#' \pkg{paradox}, \pkg{ranger}.
#'
#' The workflow for each taxonomic level:
#' \enumerate{
#'   \item Aggregate abundance by sample and taxon
#'   \item Create wide-format matrix (samples × taxa)
#'   \item Train random forest with hyperparameter tuning via random search
#'   \item Evaluate on held-out test set
#'   \item Extract variable importance
#' }
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' \dontrun{
#' result <- rf_taxa_classification(
#'   data = df_asv,
#'   taxo = df_tax,
#'   sample = df_sample,
#'   group_col = "treatment"
#' )
#' result$plot_accuracy
#' result$plot_top_features
#' result$accuracy
#' }
#'
rf_taxa_classification <- function(data,
                                    taxo,
                                    sample,
                                    tax_cols = c(
                                      phylum = "phylum", class = "class",
                                      order = "order", family = "family",
                                      genus = "genus", species = "species"
                                    ),
                                    feature_col = NULL,
                                    sample_col = "sample",
                                    group_col = NULL,
                                    tax_levels = c("OTU", names(tax_cols)),
                                    ratio = 0.7,
                                    seed = 1288,
                                    seed_tune = 8812,
                                    num_trees_range = c(1, 25),
                                    min_node_size_range = c(3, 30),
                                    term_evals = 10,
                                    cv_folds = 10,
                                    top_n = 9,
                                    fill_palette = "research",
                                    verbose = TRUE) {

  # ── Check required packages ──────────────────────────────────────────────
  pkgs <- c("mlr3", "mlr3learners", "mlr3tuning", "paradox")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Required packages not installed: ", paste(missing, collapse = ", "),
         "\nInstall with: install.packages(c('",
         paste(missing, collapse = "', '"), "'))")
  }

  # ── Input validation ─────────────────────────────────────────────────────
  if (!is.data.frame(data)) stop("'data' must be a data frame")
  if (!is.data.frame(taxo)) stop("'taxo' must be a data frame")
  if (!is.data.frame(sample)) stop("'sample' must be a data frame")
  if (is.null(group_col)) stop("'group_col' must be specified")

  # Feature column defaults to first column
  if (is.null(feature_col)) feature_col <- colnames(data)[1]

  # ── Prepare data ─────────────────────────────────────────────────────────
  # Pivot ASV table to long format
  id_cols <- feature_col
  value_cols <- setdiff(colnames(data), feature_col)

  df_long <- data %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(value_cols),
      names_to = sample_col,
      values_to = "value"
    )

  # Join taxonomy (select only relevant columns)
  taxo_select <- taxo %>%
    dplyr::select(dplyr::all_of(c(feature_col, tax_cols)))
  colnames(taxo_select)[1] <- ".feature_id"
  df_long <- df_long %>%
    dplyr::left_join(taxo_select, by = setNames(".feature_id", feature_col))

  # Join sample metadata (only group column)
  sample_select <- sample %>%
    dplyr::select(dplyr::all_of(c(sample_col, group_col)))
  df_long <- df_long %>%
    dplyr::left_join(sample_select, by = sample_col)

  # Rename for internal use
  df_long <- df_long %>%
    dplyr::rename(run = !!rlang::sym(sample_col),
                  group = !!rlang::sym(group_col),
                  otu = !!rlang::sym(feature_col))

  # Replace NA with "Unknown" in taxonomic columns
  df_long <- df_long %>% dplyr::mutate(dplyr::across(tidyselect::where(is.character),
    ~ dplyr::na_if(., "")))
  df_long <- df_long %>% tidyr::replace_na(setNames(
    rep(list("Unknown"), ncol(df_long)), colnames(df_long)
  ))

  if (verbose) message("Data prepared: ", nrow(df_long), " observations")

  # ── Iterate over taxonomic levels ────────────────────────────────────────
  all_acc <- NULL
  all_predict <- list()
  all_importance <- NULL

  for (level in tax_levels) {
    if (verbose) message("\n=== Processing: ", level, " ===")

    # Select the taxonomic column
    if (level == "OTU") {
      tax_col <- "otu"
    } else {
      tax_col <- level
      # Map to actual column name
      if (tax_col %in% names(tax_cols)) {
        actual_col <- tax_cols[tax_col]
        if (actual_col %in% colnames(df_long)) {
          df_long$tax <- df_long[[actual_col]]
        } else {
          if (verbose) message("  Skipping: column '", actual_col, "' not found")
          next
        }
      } else {
        if (verbose) message("  Skipping: level '", level, "' not in tax_cols")
        next
      }
    }
    if (level == "OTU") {
      df_long$tax <- df_long$otu
    }

    # Aggregate
    df_tmp <- df_long %>%
      dplyr::select(run, group, tax, value) %>%
      dplyr::group_by(run, tax, group) %>%
      dplyr::summarise(sum_value = sum(value), .groups = "drop")

    # Create code mapping for wide format
    distinct_tax <- df_tmp %>%
      dplyr::select(tax) %>%
      dplyr::distinct() %>%
      dplyr::mutate(tmp = paste0("tmp_", seq_len(dplyr::n())))

    # Wide format for ML
    df_ml <- df_tmp %>%
      dplyr::left_join(distinct_tax, by = "tax") %>%
      dplyr::select(run, group, tmp, sum_value) %>%
      tidyr::pivot_wider(names_from = "tmp", values_from = "sum_value",
                         values_fill = 0) %>%
      dplyr::select(-run)

    # Replace NA with 0
    df_ml[is.na(df_ml)] <- 0

    # mlr3 task
    task_rf <- mlr3::as_task_classif(df_ml, target = "group")

    # Train/test split
    set.seed(seed)
    split <- mlr3::partition(task_rf, ratio = ratio)

    # Learner
    ranger_lrn <- mlr3::lrn("classif.ranger",
                             importance = "impurity",
                             predict_type = "prob")

    # Search space
    search_space <- paradox::ps(
      num.trees = paradox::p_int(
        lower = num_trees_range[1], upper = num_trees_range[2],
        trafo = function(x) 20 * x),
      min.node.size = paradox::p_int(
        lower = min_node_size_range[1], upper = min_node_size_range[2])
    )

    # Auto tuner
    at <- mlr3tuning::auto_tuner(
      tuner = mlr3tuning::tnr("random_search"),
      learner = ranger_lrn,
      resampling = mlr3::rsmp("cv", folds = cv_folds),
      measure = mlr3::msr("classif.acc"),
      search_space = search_space,
      term_evals = term_evals
    )

    # Tune
    set.seed(seed_tune)
    at$train(task_rf, row_ids = split$train)

    # Train with best params
    ranger_lrn$param_set$values <- at$tuning_result$learner_param_vals[[1]]
    ranger_lrn$train(task_rf, row_ids = split$train)

    # Predict
    predictions <- ranger_lrn$predict(task_rf, row_ids = split$test)
    all_predict[[level]] <- predictions

    # Accuracy
    acc <- predictions$score(mlr3::msr("classif.acc"))
    all_acc <- dplyr::bind_rows(all_acc,
      data.frame(tax = level, acc = acc, stringsAsFactors = FALSE))
    if (verbose) message("  Accuracy: ", round(acc * 100, 2), "%")

    # Variable importance
    imp <- ranger_lrn$model$variable.importance
    if (!is.null(imp)) {
      imp_df <- data.frame(
        tmp = names(imp),
        importance = as.numeric(imp),
        stringsAsFactors = FALSE
      ) %>%
        dplyr::arrange(desc(importance)) %>%
        dplyr::left_join(distinct_tax, by = "tmp") %>%
        dplyr::select(tax, importance) %>%
        dplyr::mutate(tax_group = level)

      all_importance <- dplyr::bind_rows(all_importance, imp_df)
    }
  }

  # ── Accuracy plot ────────────────────────────────────────────────────────
  plot_acc <- all_acc %>%
    dplyr::arrange(dplyr::desc(acc)) %>%
    dplyr::mutate(tax = factor(tax, levels = unique(tax))) %>%
    ggplot2::ggplot(ggplot2::aes(x = tax, y = acc * 100)) +
    ggplot2::geom_col(width = 0.6, fill = "#2874A6", alpha = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = round(acc * 100, 2)),
                       vjust = -0.5, size = 4) +
    ggplot2::labs(x = NULL, y = "Prediction accuracy (%)") +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 105)) +
    ggprism::theme_prism()

  # ── Top features boxplot ─────────────────────────────────────────────────
  plot_top <- NULL
  top_features <- all_importance %>%
    dplyr::filter(tax_group == "otu") %>%
    dplyr::arrange(desc(importance)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::pull(tax)

  if (length(top_features) > 0) {
    # Get group col from sample for the boxplot
    sample_for_plot <- sample %>%
      dplyr::select(dplyr::all_of(c(sample_col, group_col)))

    top_df <- df_long %>%
      dplyr::filter(otu %in% top_features) %>%
      dplyr::select(run, group, otu, value)

    plot_top <- top_df %>%
      ggplot2::ggplot(ggplot2::aes(
        x = group, y = value, fill = group
      )) +
      ggplot2::geom_boxplot() +
      ggplot2::facet_wrap(~ otu, scales = "free", ncol = 3) +
      ggplot2::labs(x = NULL, y = "Abundance") +
      ggplot2::theme(
        legend.title = ggplot2::element_blank(),
        legend.position = "top",
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )

    if (fill_palette == "research") {
      plot_top <- plot_top + bioRtools::scale_fill_research()
    }
    plot_top <- plot_top + ggprism::theme_prism()
  }

  # ── Return ───────────────────────────────────────────────────────────────
  list(
    accuracy       = all_acc,
    predictions    = all_predict,
    importance     = all_importance,
    plot_accuracy  = plot_acc,
    plot_top_features = plot_top
  )
}
