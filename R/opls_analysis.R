#' Perform Orthogonal Partial Least Squares Discriminant Analysis (OPLS-DA)
#'
#' This function performs Orthogonal Partial Least Squares Discriminant Analysis
#' (OPLS-DA) for supervised multivariate analysis of high-dimensional data.
#' OPLS-DA is particularly useful for metabolomics, proteomics, and other omics
#' data where the goal is to identify variables that discriminate between groups
#' while filtering out orthogonal variation not related to group differences.
#'
#' @param data Numerical matrix, data frame, SummarizedExperiment, or ExpressionSet
#'   object containing the feature data. Rows represent observations (samples) and
#'   columns represent variables (features/metabolites/genes). Missing values (NA)
#'   are allowed and will be handled by the OPLS algorithm. For metabolomics data,
#'   this typically contains peak intensities or concentrations
#' @param sample Optional data frame containing sample metadata. When provided together
#'   with \code{sample_col} and \code{group_col}, the function filters \code{data}
#'   rows to match samples in \code{sample} and extracts the group vector automatically.
#'   If NULL (default), uses the \code{group} parameter directly.
#' @param sample_col Character string specifying the column in \code{sample} that
#'   contains sample IDs. Used to match rows in \code{data} (matched by row names).
#'   Default is "sample".
#' @param group_col Character string specifying the column in \code{sample} that
#'   contains group labels. Default is "group".
#' @param group Factor vector specifying group membership for discriminant analysis.
#'   Must have the same length as the number of rows in data. Ignored when
#'   \code{sample} is provided. For two-group comparisons, a binary factor is
#'   required. Multi-group OPLS-DA is supported for more than two groups
#' @param vip_threshold Numerical threshold for Variable Importance in Projection
#'   (VIP) scores (default: 1.0). Variables with VIP scores >= this threshold
#'   are considered important for group discrimination. Common thresholds:
#'   \itemize{
#'     \item 1.0: Standard threshold for variable selection
#'     \item 1.5: More stringent selection for highly important variables
#'     \item 0.8: More liberal selection for exploratory analysis
#'   }
#' @param ortho_components Number of orthogonal components to calculate (default: 1).
#'   Orthogonal components capture variation unrelated to group differences.
#'   Typically 1-3 components are sufficient for most datasets
#' @param pred_components Number of predictive components (default: 1 for binary,
#'   automatically determined for multi-group). Usually 1 for two-group comparisons,
#'   more may be needed for complex multi-group analyses
#' @param scaling Method for data scaling (default: "standard"). Options include:
#'   \itemize{
#'     \item "none": No scaling (use when data is already appropriately scaled)
#'     \item "center": Mean centering only
#'     \item "standard": Mean centering and unit variance scaling (recommended)
#'     \item "pareto": Pareto scaling (square root of standard deviation)
#'   }
#' @param validation Cross-validation method (default: "CV"). Options:
#'   \itemize{
#'     \item "CV": Cross-validation for model assessment
#'     \item "none": No validation (faster but no performance metrics)
#'   }
#' @param cv_folds Number of cross-validation folds (default: 7). Only used when
#'   validation = "CV". Should be <= number of samples in smallest group
#' @param ref_group Character string specifying the reference group for differential
#'   analysis. When specified, each non-reference group is compared against this group
#'   using t-tests with log2FC calculation. If NULL (default), no differential analysis
#'   is performed. Similar to \code{rstatix::t_test(ref.group = )}.
#' @param p_threshold Numeric p-value threshold for significance labeling (default: 0.05).
#'   Used in conjunction with log2FC direction to label features as Up/Down.
#' @param verbose Logical, whether to print progress messages during pairwise model fitting
#'   (default: FALSE).
#'
#' @return A named list containing:
#'   \describe{
#'     \item{\code{model}}{OPLS model object from ropls package. In pairwise mode,
#'       this is the first pairwise model.}
#'     \item{\code{scores}}{Data frame with sample scores for visualization.
#'       In pairwise mode, includes a \code{comparison} column.}
#'     \item{\code{vip_scores}}{Data frame with VIP scores, sorted by VIP descending.
#'       In pairwise mode, includes \code{group} and \code{ref_group} columns.}
#'     \item{\code{loadings}}{Data frame with variable loadings.}
#'     \item{\code{differential_analysis}}{Data frame with log2FC, p-values, and
#'       regulation labels for each group vs ref_group. NULL if no ref_group.}
#'     \item{\code{model_summary}}{List with R2Y, Q2Y, and other model metrics.
#'       In pairwise mode, includes per-comparison results.}
#'   }
#'
#' @details
#' \strong{OPLS-DA Method:}
#' OPLS-DA extends PLS-DA by separating predictive variation (correlated with Y)
#' from orthogonal variation (uncorrelated with Y). This provides:
#' \itemize{
#'   \item Better interpretation of group differences
#'   \item Cleaner visualization with reduced model complexity
#'   \item Improved identification of discriminating variables
#' }
#'
#' \strong{Model Interpretation:}
#' \itemize{
#'   \item \strong{R2Y}: Proportion of Y-variance explained by the model
#'   \item \strong{Q2Y}: Cross-validated predictive ability (>0.5 indicates good model)
#'   \item \strong{VIP scores}: Variable importance (>1.0 indicates important variables)
#'   \item \strong{Score plots}: Visualization of sample separation
#'   \item \strong{Loading plots}: Variable contributions to components
#' }
#'
#' \strong{Data Preprocessing Recommendations:}
#' \itemize{
#'   \item Remove variables with >50% missing values
#'   \item Consider log-transformation for skewed data (e.g., metabolomics)
#'   \item Use standard scaling for variables with different units
#'   \item Filter out low-variance variables if appropriate
#' }
#'
#' @note
#' \itemize{
#'   \item OPLS-DA is a supervised method; avoid overfitting with small sample sizes
#'   \item Cross-validation is essential for model validation
#'   \item VIP scores help identify the most discriminating variables
#'   \item For metabolomics: consider sample normalization before analysis
#'   \item Minimum 6-10 samples per group recommended for stable models
#' }
#'
#' @references
#' Trygg, J. and Wold, S. (2002). Orthogonal projections to latent structures (O-PLS).
#' Journal of Chemometrics, 16(3), 119-128.
#'
#' Bylesjö, M., et al. (2006). OPLS discriminant analysis: combining the strengths of
#' PLS‐DA and SIMCA classification. Journal of Chemometrics, 20(8‐10), 341-351.
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{spls_analysis}} for sparse PLS-DA analysis
#'   \item \code{\link{pca_analysis}} for unsupervised dimensionality reduction
#'   \item \code{\link[ropls]{opls}} for underlying OPLS implementation
#' }
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' library(bioRtools)
#' library(dplyr)
#' library(ggplot2)
#'
#' # Example 1: Basic OPLS-DA analysis
#' \dontrun{
#' # Load example metabolomics data
#' data("df.splsda.meta")     # Metabolite intensity matrix
#' data("df.splsda.sample")   # Sample metadata
#'
#' # Prepare data for OPLS-DA (filter to specific timepoint and treatments)
#' sample_filtered <- df.splsda.sample %>%
#'   filter(day == 180, treatment != "Cover") %>%
#'   mutate(treatment = factor(treatment))
#'
#' # Filter metabolite data to match samples
#' meta_filtered <- df.splsda.meta %>%
#'   filter(rownames(.) %in% sample_filtered$sample)
#'
#' # Perform OPLS-DA analysis
#' oplsda_results <- opls_analysis(
#'   data = meta_filtered,
#'   group = sample_filtered$treatment
#' )
#'
#' # View model summary
#' print("OPLS-DA Model Summary:")
#' print(oplsda_results$model_summary)
#'
#' # Check important variables
#' important_vars <- oplsda_results$vip_scores %>%
#'   filter(important == TRUE) %>%
#'   head(10)
#'
#' print("Top 10 important metabolites:")
#' print(important_vars)
#' }
#'
#' # Example 2: Custom VIP threshold and validation
#' \dontrun{
#' # Stringent analysis with higher VIP threshold
#' oplsda_strict <- opls_analysis(
#'   data = meta_filtered,
#'   group = sample_filtered$treatment,
#'   vip_threshold = 1.5,        # More stringent VIP cutoff
#'   scaling = "pareto",          # Pareto scaling for metabolomics
#'   cv_folds = 5                 # 5-fold cross-validation
#' )
#'
#' print("Strict model performance:")
#' print(paste("R2Y:", round(oplsda_strict$model_summary$R2Y, 3)))
#' print(paste("Q2Y:", round(oplsda_strict$model_summary$Q2Y, 3)))
#' print(paste("Important variables:", sum(oplsda_strict$vip_scores$important)))
#' }
#'
#' # Example 3: Multi-group OPLS-DA
#' \dontrun{
#' # Prepare multi-group data
#' multi_group_samples <- df.splsda.sample %>%
#'   filter(day %in% c(90, 180)) %>%
#'   mutate(
#'     group = paste(treatment, day, sep = "_"),
#'     group = factor(group)
#'   )
#'
#' multi_group_meta <- df.splsda.meta %>%
#'   filter(rownames(.) %in% multi_group_samples$sample)
#'
#' # Multi-group OPLS-DA
#' multigroup_oplsda <- opls_analysis(
#'   data = multi_group_meta,
#'   group = multi_group_samples$group,
#'   pred_components = 2,         # May need more components for multiple groups
#'   vip_threshold = 1.2
#' )
#'
#' print("Multi-group model summary:")
#' print(multigroup_oplsda$model_summary)
#'
#' # Visualize group separation
#' score_plot <- multigroup_oplsda$scores %>%
#'   ggplot(aes(x = t1, y = to1, color = group)) +
#'   geom_point(size = 3, alpha = 0.7) +
#'   stat_ellipse(level = 0.95) +
#'   labs(
#'     x = "Predictive Component 1",
#'     y = "Orthogonal Component 1",
#'     title = "OPLS-DA Score Plot",
#'     subtitle = "95% confidence ellipses"
#'   ) +
#'   theme_minimal()
#'
#' print(score_plot)
#' }
#'
#' # Example 4: Variable importance analysis and biomarker discovery
#' \dontrun{
#' # Analyze variable importance patterns
#' vip_analysis <- oplsda_results$vip_scores %>%
#'   mutate(
#'     importance_level = case_when(
#'       vip >= 2.0 ~ "Very High",
#'       vip >= 1.5 ~ "High",
#'       vip >= 1.0 ~ "Moderate",
#'       TRUE ~ "Low"
#'     )
#'   )
#'
#' # Summary by importance level
#' importance_summary <- vip_analysis %>%
#'   count(importance_level) %>%
#'   arrange(desc(n))
#'
#' print("Variable importance distribution:")
#' print(importance_summary)
#'
#' # Create VIP plot
#' vip_plot <- vip_analysis %>%
#'   filter(vip >= 0.8) %>%  # Show moderately important and above
#'   slice_head(n = 20) %>%   # Top 20 variables
#'   ggplot(aes(x = reorder(feature, vip), y = vip, fill = importance_level)) +
#'   geom_col() +
#'   geom_hline(yintercept = 1.0, linetype = "dashed", color = "red") +
#'   coord_flip() +
#'   labs(
#'     x = "Metabolites",
#'     y = "VIP Score",
#'     title = "Variable Importance in Projection (VIP)",
#'     subtitle = "Top 20 discriminating metabolites",
#'     fill = "Importance"
#'   ) +
#'   theme_minimal() +
#'   theme(axis.text.y = element_text(size = 8))
#'
#' print(vip_plot)
#' }
#'
#' # Example 5: Model validation and diagnostics
#' \dontrun{
#' # Extract detailed model metrics
#' model_object <- oplsda_results$model
#'
#' # Cross-validation results
#' cv_metrics <- data.frame(
#'   Metric = c("R2Y", "Q2Y", "RMSEE", "RMSECV"),
#'   Value = c(
#'     model_object@summaryDF$R2Y,
#'     model_object@summaryDF$Q2Y,
#'     model_object@summaryDF$RMSEE,
#'     model_object@summaryDF$RMSECV
#'   )
#' )
#'
#' print("Cross-validation metrics:")
#' print(cv_metrics)
#'
#' # Model interpretation guidelines
#' cat("Model Interpretation Guidelines:\n")
#' cat("R2Y > 0.5: Good explanatory power\n")
#' cat("Q2Y > 0.5: Good predictive ability\n")
#' cat("Q2Y > 0.9: Excellent predictive ability\n")
#' cat("VIP > 1.0: Important for discrimination\n")
#' cat("VIP > 1.5: Highly important variables\n\n")
#'
#' # Check for overfitting
#' if (cv_metrics$Value[cv_metrics$Metric == "Q2Y"] < 0.5) {
#'   warning("Q2Y < 0.5: Potential overfitting. Consider reducing model complexity.")
#' }
#'
#' if ((cv_metrics$Value[cv_metrics$Metric == "R2Y"] -
#'   cv_metrics$Value[cv_metrics$Metric == "Q2Y"]) > 0.3) {
#'   warning("Large R2Y-Q2Y gap: Possible overfitting detected.")
#' }
#' }
#'
#' # Example 6: Comparative analysis with different scaling methods
#' \dontrun{
#' scaling_methods <- c("standard", "pareto", "center")
#' scaling_results <- list()
#'
#' for (method in scaling_methods) {
#'   scaling_results[[method]] <- opls_analysis(
#'     data = meta_filtered,
#'     group = sample_filtered$treatment,
#'     scaling = method,
#'     validation = "CV"
#'   )
#' }
#'
#' # Compare model performance
#' scaling_comparison <- data.frame(
#'   Scaling = scaling_methods,
#'   R2Y = sapply(scaling_results, function(x) x$model_summary$R2Y),
#'   Q2Y = sapply(scaling_results, function(x) x$model_summary$Q2Y),
#'   n_important = sapply(scaling_results, function(x) sum(x$vip_scores$important))
#' )
#'
#' print("Scaling method comparison:")
#' print(scaling_comparison)
#'
#' # Recommend best scaling method
#' best_scaling <- scaling_comparison[which.max(scaling_comparison$Q2Y), "Scaling"]
#' print(paste("Best scaling method based on Q2Y:", best_scaling))
#' }
#'
opls_analysis <- function(data, sample = NULL, sample_col = "sample",
                          group_col = "group", group = NULL,
                          vip_threshold = 1.0, ortho_components = 1,
                          pred_components = NULL, scaling = "standard",
                          validation = "CV", cv_folds = 7,
                          ref_group = NULL, p_threshold = 0.05,
                          verbose = FALSE) {
  # Input validation
  if (!is.matrix(data) && !is.data.frame(data) &&
    !methods::is(data, "SummarizedExperiment") &&
    !methods::is(data, "ExpressionSet")) {
    stop("'data' must be a matrix, data.frame, SummarizedExperiment, or ExpressionSet")
  }

  # If sample data frame is provided, filter data and extract group
  if (!is.null(sample)) {
    if (!is.data.frame(sample)) {
      stop("'sample' must be a data frame")
    }
    if (!sample_col %in% colnames(sample)) {
      stop("'sample_col' (", sample_col, ") not found in 'sample'")
    }
    if (!group_col %in% colnames(sample)) {
      stop("'group_col' (", group_col, ") not found in 'sample'")
    }

    # Convert data to matrix first to get row names
    if (is.data.frame(data)) {
      data <- as.matrix(data)
    } else if (methods::is(data, "SummarizedExperiment")) {
      data <- SummarizedExperiment::assay(data)
    } else if (methods::is(data, "ExpressionSet")) {
      data <- Biobase::exprs(data)
    }

    sample_ids <- sample[[sample_col]]
    data_rownames <- rownames(data)
    if (is.null(data_rownames)) {
      stop("'data' must have row names matching sample IDs when using 'sample' parameter")
    }

    # Match and filter
    common <- intersect(sample_ids, data_rownames)
    if (length(common) == 0) {
      stop("No matching samples between 'data' row names and 'sample'", call. = FALSE)
    }

    sample <- sample[match(common, sample[[sample_col]]), ]
    data <- data[common, , drop = FALSE]
    group <- sample[[group_col]]
    sample_info <- sample

    if (nrow(data) < length(sample_ids)) {
      warning(sprintf("Matched %d/%d samples from 'sample'", nrow(data), length(sample_ids)))
    }
  } else {
    # Use group vector directly
    if (is.null(group)) {
      stop("Either 'sample' or 'group' must be provided")
    }
    if (!is.factor(group) && !is.character(group)) {
      stop("'group' must be a factor or character vector")
    }
    sample_info <- NULL
  }

  # Convert data to matrix if needed
  if (!is.null(sample)) {
    # Already converted to matrix in sample branch above
    data_matrix <- data
  } else if (is.data.frame(data)) {
    data_matrix <- as.matrix(data)
  } else if (methods::is(data, "SummarizedExperiment")) {
    data_matrix <- SummarizedExperiment::assay(data)
  } else if (methods::is(data, "ExpressionSet")) {
    data_matrix <- Biobase::exprs(data)
  } else {
    data_matrix <- data
  }

  # Validate dimensions
  if (nrow(data_matrix) != length(group)) {
    stop("Number of rows in 'data' must equal length of 'group'")
  }

  if (nrow(data_matrix) == 0 || ncol(data_matrix) == 0) {
    stop("'data' cannot be empty")
  }

  # Convert group to factor and validate
  if (is.character(group)) {
    group <- factor(group)
  }

  group_levels <- levels(group)
  n_groups <- length(group_levels)

  if (n_groups < 2) {
    stop("'group' must have at least 2 levels for discriminant analysis")
  }

  # Check group sizes
  group_sizes <- table(group)
  min_group_size <- min(group_sizes)

  if (min_group_size < 3) {
    warning(paste("Small group sizes detected. Minimum group size:", min_group_size,
      ". Results may be unreliable with < 6 samples per group."))
  }

  # Validate parameters
  if (!is.numeric(vip_threshold) || length(vip_threshold) != 1 || vip_threshold < 0) {
    stop("'vip_threshold' must be a single non-negative number")
  }

  if (!is.numeric(ortho_components) || length(ortho_components) != 1 ||
    ortho_components < 0 || ortho_components != round(ortho_components)) {
    stop("'ortho_components' must be a single non-negative integer")
  }

  if (!is.null(pred_components)) {
    if (!is.numeric(pred_components) || length(pred_components) != 1 ||
      pred_components < 1 || pred_components != round(pred_components)) {
      stop("'pred_components' must be a single positive integer")
    }
  } else {
    # Auto-determine predictive components
    pred_components <- min(n_groups - 1, 3)  # Usually sufficient
  }

  valid_scaling <- c("none", "center", "standard", "pareto")
  if (!scaling %in% valid_scaling) {
    stop(paste("'scaling' must be one of:", paste(valid_scaling, collapse = ", ")))
  }

  valid_validation <- c("CV", "none")
  if (!validation %in% valid_validation) {
    stop("'validation' must be 'CV' or 'none'")
  }

  if (!is.null(ref_group)) {
    if (!ref_group %in% group_levels) {
      stop("'ref_group' must be one of the group levels: ", paste(group_levels, collapse = ", "))
    }
  }

  if (!is.numeric(p_threshold) || length(p_threshold) != 1 ||
    p_threshold <= 0 || p_threshold >= 1) {
    stop("'p_threshold' must be a single number between 0 and 1")
  }

  if (validation == "CV") {
    if (!is.numeric(cv_folds) || length(cv_folds) != 1 ||
      cv_folds < 2 || cv_folds != round(cv_folds)) {
      stop("'cv_folds' must be a single integer >= 2")
    }

    if (cv_folds > min_group_size) {
      warning(paste("cv_folds (", cv_folds, ") > smallest group size (", min_group_size,
        "). Reducing to", min_group_size))
      cv_folds <- min_group_size
    }
  }

  # Check for missing values
  na_count <- sum(is.na(data_matrix))
  if (na_count > 0) {
    warning(paste("Data contains", na_count, "missing values.",
      "OPLS will handle these, but consider imputation for better results."))
  }

  # Check for zero variance variables
  var_check <- apply(data_matrix, 2, var, na.rm = TRUE)
  zero_var_count <- sum(var_check == 0, na.rm = TRUE)

  if (zero_var_count > 0) {
    warning(paste("Found", zero_var_count, "variables with zero variance.",
      "Consider removing these before analysis."))
  }

  # Prepare sample names for output
  sample_names <- rownames(data_matrix)
  if (is.null(sample_names)) {
    sample_names <- paste0("Sample_", 1:nrow(data_matrix))
  }

  # Prepare variable names
  variable_names <- colnames(data_matrix)
  if (is.null(variable_names)) {
    variable_names <- paste0("Var_", 1:ncol(data_matrix))
  }

  # Run OPLS-DA analysis
  # When ref_group is specified with >2 groups, fit pairwise binary OPLS-DA models
  pairwise_models <- NULL
  if (!is.null(ref_group) && n_groups > 2) {
    # Pairwise mode: one OPLS-DA per (group vs ref_group)
    other_groups <- setdiff(group_levels, ref_group)
    pairwise_models <- list()
    all_scores <- list()
    all_vip <- list()
    all_summaries <- list()

    for (grp in other_groups) {
      idx <- which(group %in% c(ref_group, grp))
      sub_data <- data_matrix[idx, , drop = FALSE]
      sub_group <- droplevels(factor(as.character(group[idx])))

      grp_cv_folds <- cv_folds
      min_grp_size <- min(table(sub_group))
      if (grp_cv_folds > min_grp_size) {
        grp_cv_folds <- min_grp_size
      }

      if (verbose) {
        message("Fitting OPLS-DA: ", grp, " vs ", ref_group,
          " (n=", nrow(sub_data), ", cv_folds=", grp_cv_folds, ")")
      }

      model <- tryCatch(
        {
          if (validation == "CV") {
            ropls::opls(sub_data, sub_group,
              predI = 1, orthoI = ortho_components, scaleC = scaling,
              crossvalI = grp_cv_folds, fig.pdfC = "none", info.txtC = "none")
          } else {
            ropls::opls(sub_data, sub_group,
              predI = 1, orthoI = ortho_components, scaleC = scaling,
              crossvalI = 0, fig.pdfC = "none", info.txtC = "none")
          }
        },
        error = function(e) {
          warning("OPLS-DA failed for ", grp, " vs ", ref_group, ": ", e$message)
          NULL
        })

      if (!is.null(model)) {
        pairwise_models[[grp]] <- model

        # Extract scores for this pair
        tryCatch(
          {
            pair_scores <- chemhelper::get_scores(model) %>%
              as.data.frame()
            pair_scores$sample_id <- rownames(sub_data)
            pair_scores$group <- as.character(sub_group)
            pair_scores$comparison <- paste0(grp, " vs ", ref_group)
            all_scores[[grp]] <- pair_scores
          },
          error = function(e) {
            pair_scores <- data.frame(
              t1 = model@scoreMN[, 1],
              sample_id = rownames(sub_data),
              group = as.character(sub_group),
              comparison = paste0(grp, " vs ", ref_group)
            )
            if (!is.null(model@orthoScoreMN) && ncol(model@orthoScoreMN) > 0) {
              pair_scores$to1 <- model@orthoScoreMN[, 1]
            }
            all_scores[[grp]] <- pair_scores
          })

        # Extract VIP
        tryCatch(
          {
            vip_vals <- model@vipVn
            all_vip[[grp]] <- data.frame(
              feature = variable_names,
              vip = as.numeric(vip_vals),
              group = grp,
              ref_group = ref_group,
              stringsAsFactors = FALSE
            )
          },
          error = function(e) {
            all_vip[[grp]] <- data.frame(
              feature = variable_names,
              vip = NA_real_,
              group = grp,
              ref_group = ref_group,
              stringsAsFactors = FALSE
            )
          })

        # Extract summary
        tryCatch(
          {
            s <- model@summaryDF
            all_summaries[[grp]] <- list(
              group = grp,
              R2Y = if ("R2Y" %in% names(s)) s$R2Y else NA,
              Q2Y = if ("Q2Y" %in% names(s)) s$Q2Y else NA,
              R2X = if ("R2X(cum)" %in% names(s)) s$`R2X(cum)` else NA
            )
          },
          error = function(e) NULL)
      }
    }

    # Combine results
    opls_model <- pairwise_models[[1]]  # Use first model as representative
    scores_data <- do.call(rbind, all_scores)

    # Merge sample metadata into scores if available
    if (!is.null(sample_info)) {
      scores_data$group <- NULL
      scores_data <- scores_data %>%
        dplyr::left_join(sample_info, by = stats::setNames(sample_col, "sample_id"))
    }

    vip_data <- do.call(rbind, all_vip) %>%
      dplyr::arrange(desc(vip))
    vip_data$important <- vip_data$vip >= vip_threshold

    # Build model_summary from pairwise summaries
    model_summary <- list(
      mode = "pairwise (ref_group)",
      ref_group = ref_group,
      pairwise_results = all_summaries,
      R2Y = mean(sapply(all_summaries, function(x) x$R2Y), na.rm = TRUE),
      Q2Y = mean(sapply(all_summaries, function(x) x$Q2Y), na.rm = TRUE),
      n_groups = n_groups,
      n_variables = ncol(data_matrix),
      n_samples = nrow(data_matrix)
    )

    loadings_data <- data.frame(feature = variable_names)

    # Pairwise differential analysis (log2FC + p-value for each comparison)
    diff_analysis <- NULL
    if (!is.null(ref_group)) {
      diff_list <- list()
      ref_idx <- which(group == ref_group)

      for (grp in other_groups) {
        grp_idx <- which(group == grp)
        ref_data <- data_matrix[ref_idx, , drop = FALSE]
        grp_data <- data_matrix[grp_idx, , drop = FALSE]

        n_ref <- nrow(ref_data)
        n_grp <- nrow(grp_data)
        use_wilcox <- n_ref < 6 || n_grp < 6

        log2fc_vals <- numeric(length(variable_names))
        p_vals <- numeric(length(variable_names))

        for (j in seq_along(variable_names)) {
          ref_vals <- ref_data[, j]
          grp_vals <- grp_data[, j]
          ref_mean <- mean(ref_vals, na.rm = TRUE)
          grp_mean <- mean(grp_vals, na.rm = TRUE)
          log2fc_vals[j] <- log2((grp_mean + 1e-10) / (ref_mean + 1e-10))

          valid_ref <- ref_vals[!is.na(ref_vals)]
          valid_grp <- grp_vals[!is.na(grp_vals)]
          if (length(valid_ref) >= 2 && length(valid_grp) >= 2 &&
            stats::var(valid_ref) > 0 && stats::var(valid_grp) > 0) {
            if (use_wilcox) {
              test_result <- stats::wilcox.test(valid_grp, valid_ref)
            } else {
              test_result <- stats::t.test(valid_grp, valid_ref)
            }
            p_vals[j] <- test_result$p.value
          } else {
            p_vals[j] <- NA
          }
        }

        diff_list[[grp]] <- data.frame(
          feature = variable_names,
          group = grp,
          ref_group = ref_group,
          log2_fc = round(log2fc_vals, 4),
          p_value = p_vals,
          stringsAsFactors = FALSE
        )
      }

      diff_analysis <- do.call(rbind, diff_list)
      rownames(diff_analysis) <- NULL
      diff_analysis$p_adjust <- stats::p.adjust(diff_analysis$p_value, method = "BH")

      # Add VIP from pairwise results
      vip_lookup <- vip_data %>%
        dplyr::select(feature, group, vip)
      diff_analysis <- diff_analysis %>%
        dplyr::left_join(vip_lookup, by = c("feature", "group"))
      diff_analysis$vip <- round(diff_analysis$vip, 4)

      diff_analysis$regulation <- ifelse(
        is.na(diff_analysis$p_adjust), "NS",
        ifelse(diff_analysis$log2_fc > 0 & diff_analysis$p_adjust < p_threshold, "Up",
          ifelse(diff_analysis$log2_fc < 0 & diff_analysis$p_adjust < p_threshold, "Down", "NS")
        )
      )
      reg_order <- c("Up", "Down", "NS")
      diff_analysis$regulation <- factor(diff_analysis$regulation, levels = reg_order)
      diff_analysis <- diff_analysis[order(diff_analysis$regulation, -abs(diff_analysis$log2_fc)), ]
      diff_analysis$regulation <- as.character(diff_analysis$regulation)
    }

  } else {
    # Standard single-model mode (2 groups or no ref_group)
    tryCatch(
      {
        if (validation == "CV") {
          opls_model <- ropls::opls(
            x = data_matrix,
            y = group,
            predI = pred_components,
            orthoI = ortho_components,
            scaleC = scaling,
            crossvalI = cv_folds,
            fig.pdfC = "none",
            info.txtC = "none"
          )
        } else {
          opls_model <- ropls::opls(
            x = data_matrix,
            y = group,
            predI = pred_components,
            orthoI = ortho_components,
            scaleC = scaling,
            crossvalI = 0,
            fig.pdfC = "none",
            info.txtC = "none"
          )
        }
      },
      error = function(e) {
        stop(paste("OPLS model fitting failed:", e$message,
          "\nTry reducing the number of components or checking data quality."))
      })

    # Extract scores for visualization
    tryCatch(
      {
        scores_data <- chemhelper::get_scores(opls_model) %>%
          as.data.frame() %>%
          dplyr::mutate(
            sample_id = sample_names,
            group = group
          )
      },
      error = function(e) {
        warning(paste("Error extracting scores with chemhelper:", e$message))
        scores_data <- data.frame(
          t1 = opls_model@scoreMN[, 1],
          sample_id = sample_names,
          group = group
        )

        if (!is.null(opls_model@orthoScoreMN) && ncol(opls_model@orthoScoreMN) > 0) {
          scores_data$to1 <- opls_model@orthoScoreMN[, 1]
        }
      })

    # Merge sample metadata into scores if available
    if (!is.null(sample_info)) {
      scores_data$group <- NULL
      scores_data <- scores_data %>%
        dplyr::left_join(sample_info, by = stats::setNames(sample_col, "sample_id"))
    }

    # Annotate score columns with variance info
    tryCatch(
      {
        r2x_cum <- opls_model@summaryDF$`R2X(cum)`
        r2y_cum <- opls_model@summaryDF$`R2Y(cum)`
        q2_cum <- opls_model@summaryDF$`Q2(cum)`

        score_cols <- grep("^(t|to|p|o)\\d+$", names(scores_data), value = TRUE)
        if (length(score_cols) > 0) {
          pct_str <- paste0("R2X=", round(r2x_cum * 100, 1), "% R2Y=", round(r2y_cum * 100, 1), "% Q2=", round(q2_cum * 100, 1), "%")
          for (col in score_cols) {
            label <- if (grepl("^(t|p)", col)) {
              paste0(col, " (", pct_str, ")")
            } else {
              paste0(col, " (orthogonal)")
            }
            names(scores_data)[names(scores_data) == col] <- label
          }
        }
      },
      error = function(e) {
        # Keep original column names if annotation fails
      })

    # Extract VIP scores
    vip_data <- tryCatch(
      {
        vip_values <- opls_model@vipVn

        if (is.null(vip_values)) {
          warning("VIP scores not available in model")
          data.frame(
            feature = variable_names,
            vip = rep(NA, length(variable_names)),
            important = rep(FALSE, length(variable_names))
          )
        } else {
          data.frame(
            feature = variable_names,
            vip = as.numeric(vip_values),
            important = as.numeric(vip_values) >= vip_threshold
          ) %>%
            dplyr::arrange(desc(vip))
        }
      },
      error = function(e) {
        warning(paste("Error extracting VIP scores:", e$message))
        data.frame(
          feature = variable_names,
          vip = rep(NA, length(variable_names)),
          important = rep(FALSE, length(variable_names))
        )
      })

    # Extract loadings
    loadings_data <- tryCatch(
      {
        if (!is.null(opls_model@loadingMN)) {
          loading_df <- as.data.frame(opls_model@loadingMN)
          loading_df$feature <- variable_names
          loading_df
        } else {
          data.frame(feature = variable_names)
        }
      },
      error = function(e) {
        warning(paste("Error extracting loadings:", e$message))
        data.frame(feature = variable_names)
      })

    # Extract model summary statistics
    model_summary <- tryCatch(
      {
        summary_df <- opls_model@summaryDF

        list(
          R2Y = if ("R2Y" %in% names(summary_df)) summary_df$R2Y else NA,
          Q2Y = if ("Q2Y" %in% names(summary_df)) summary_df$Q2Y else NA,
          RMSEE = if ("RMSEE" %in% names(summary_df)) summary_df$RMSEE else NA,
          RMSECV = if ("RMSECV" %in% names(summary_df)) summary_df$RMSECV else NA,
          n_pred_components = pred_components,
          n_ortho_components = ortho_components,
          n_variables = ncol(data_matrix),
          n_samples = nrow(data_matrix),
          n_groups = n_groups,
          group_sizes = as.list(group_sizes),
          scaling_method = scaling,
          validation_method = validation
        )
      },
      error = function(e) {
        warning(paste("Error extracting model summary:", e$message))
        list(
          R2Y = NA, Q2Y = NA, RMSEE = NA, RMSECV = NA,
          n_pred_components = pred_components,
          n_ortho_components = ortho_components,
          n_variables = ncol(data_matrix),
          n_samples = nrow(data_matrix),
          n_groups = n_groups
        )
      })
  } # end of else (standard single-model mode)

  # Count important variables
  n_important_vars <- sum(vip_data$important, na.rm = TRUE)

  # Differential analysis against ref_group (single-model mode only)
  diff_analysis <- NULL
  if (!is.null(ref_group) && is.null(pairwise_models)) {
    other_groups <- setdiff(group_levels, ref_group)
    ref_idx <- which(group == ref_group)

    diff_list <- list()
    vip_lookup <- stats::setNames(vip_data$vip, vip_data$feature)

    for (grp in other_groups) {
      grp_idx <- which(group == grp)
      ref_data <- data_matrix[ref_idx, , drop = FALSE]
      grp_data <- data_matrix[grp_idx, , drop = FALSE]

      n_ref <- nrow(ref_data)
      n_grp <- nrow(grp_data)
      use_wilcox <- n_ref < 6 || n_grp < 6

      diff_list[[grp]] <- data.frame(
        feature = variable_names,
        group = grp,
        ref_group = ref_group,
        stringsAsFactors = FALSE
      )

      # Calculate log2FC and p-value per variable
      log2fc_vals <- numeric(length(variable_names))
      p_vals <- numeric(length(variable_names))

      for (j in seq_along(variable_names)) {
        ref_vals <- ref_data[, j]
        grp_vals <- grp_data[, j]

        ref_mean <- mean(ref_vals, na.rm = TRUE)
        grp_mean <- mean(grp_vals, na.rm = TRUE)

        log2fc_vals[j] <- log2((grp_mean + 1e-10) / (ref_mean + 1e-10))

        valid_ref <- ref_vals[!is.na(ref_vals)]
        valid_grp <- grp_vals[!is.na(grp_vals)]

        if (length(valid_ref) >= 2 && length(valid_grp) >= 2 &&
          stats::var(valid_ref) > 0 && stats::var(valid_grp) > 0) {
          if (use_wilcox) {
            test_result <- stats::wilcox.test(valid_grp, valid_ref)
          } else {
            test_result <- stats::t.test(valid_grp, valid_ref)
          }
          p_vals[j] <- test_result$p.value
        } else {
          p_vals[j] <- NA
        }
      }

      diff_list[[grp]]$log2_fc <- round(log2fc_vals, 4)
      diff_list[[grp]]$p_value <- p_vals
    }

    diff_analysis <- do.call(rbind, diff_list)
    rownames(diff_analysis) <- NULL

    # Adjust p-values across all comparisons
    diff_analysis$p_adjust <- stats::p.adjust(diff_analysis$p_value, method = "BH")

    # Add VIP scores
    diff_analysis$vip <- round(as.numeric(vip_lookup[diff_analysis$feature]), 4)

    # Label regulation
    diff_analysis$regulation <- ifelse(
      is.na(diff_analysis$p_adjust), "NS",
      ifelse(diff_analysis$log2_fc > 0 & diff_analysis$p_adjust < p_threshold, "Up",
        ifelse(diff_analysis$log2_fc < 0 & diff_analysis$p_adjust < p_threshold, "Down", "NS")
      )
    )

    # Sort by regulation significance then log2_fc
    reg_order <- c("Up", "Down", "NS")
    diff_analysis$regulation <- factor(diff_analysis$regulation, levels = reg_order)
    diff_analysis <- diff_analysis[order(diff_analysis$regulation, -abs(diff_analysis$log2_fc)), ]
    diff_analysis$regulation <- as.character(diff_analysis$regulation)
  }

  # Prepare final results
  results <- list(
    model = opls_model,
    scores = scores_data,
    vip_scores = vip_data,
    loadings = loadings_data,
    differential_analysis = diff_analysis,
    model_summary = model_summary
  )

  # Add metadata as attributes
  attr(results, "analysis_type") <- "OPLS-DA"
  attr(results, "n_samples") <- nrow(data_matrix)
  attr(results, "n_variables") <- ncol(data_matrix)
  attr(results, "n_groups") <- n_groups
  attr(results, "n_important_vars") <- n_important_vars
  attr(results, "vip_threshold") <- vip_threshold
  attr(results, "scaling_method") <- scaling

  # Print summary if interactive
  if (interactive()) {
    cat("OPLS-DA Analysis Summary:\n")
    cat("========================\n")
    cat("Samples:", nrow(data_matrix), "\n")
    cat("Variables:", ncol(data_matrix), "\n")
    cat("Groups:", n_groups, "(", paste(names(group_sizes), collapse = ", "), ")\n")
    cat("Group sizes:", paste(group_sizes, collapse = ", "), "\n")
    cat("Predictive components:", pred_components, "\n")
    cat("Orthogonal components:", ortho_components, "\n")
    cat("Scaling method:", scaling, "\n\n")

    if (!is.na(model_summary$R2Y) && !is.na(model_summary$Q2Y)) {
      cat("Model Performance:\n")
      cat("R2Y (explained variance):", round(model_summary$R2Y, 3), "\n")
      cat("Q2Y (predictive ability):", round(model_summary$Q2Y, 3), "\n")

      if (model_summary$Q2Y > 0.9) {
        cat("Model quality: Excellent\n")
      } else if (model_summary$Q2Y > 0.5) {
        cat("Model quality: Good\n")
      } else {
        cat("Model quality: Poor (possible overfitting)\n")
      }
    }

    cat("\nVariable Selection:\n")
    cat("VIP threshold:", vip_threshold, "\n")
    cat("Important variables:", n_important_vars,
      paste0("(", round(100 * n_important_vars / ncol(data_matrix), 1), "%)"), "\n")

    if (n_important_vars > 0) {
      cat("\nTop 5 important variables:\n")
      top_vip <- head(vip_data[vip_data$important, ], 5)
      for (i in 1:nrow(top_vip)) {
        cat(sprintf("  %s: VIP = %.3f\n", top_vip$feature[i], top_vip$vip[i]))
      }
    }

    if (!is.null(diff_analysis)) {
      cat("\nDifferential Analysis (ref_group =", ref_group, "):\n")
      cat("------------------------------------------\n")
      other_groups <- unique(diff_analysis$group)
      for (grp in other_groups) {
        subset <- diff_analysis[diff_analysis$group == grp, ]
        n_up <- sum(subset$regulation == "Up", na.rm = TRUE)
        n_down <- sum(subset$regulation == "Down", na.rm = TRUE)
        cat(sprintf("  %s vs %s: %d Up, %d Down, %d NS\n",
          grp, ref_group, n_up, n_down, sum(subset$regulation == "NS", na.rm = TRUE)))
      }
    }
    cat("\n")
  }

  results
}
