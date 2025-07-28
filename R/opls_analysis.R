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
#' @param group Factor vector specifying group membership for discriminant analysis.
#'   Must have the same length as the number of rows in data. For two-group
#'   comparisons, a binary factor is required. Multi-group OPLS-DA is supported
#'   for more than two groups
#' @param vip.threshold Numerical threshold for Variable Importance in Projection
#'   (VIP) scores (default: 1.0). Variables with VIP scores >= this threshold
#'   are considered important for group discrimination. Common thresholds:
#'   \itemize{
#'     \item 1.0: Standard threshold for variable selection
#'     \item 1.5: More stringent selection for highly important variables
#'     \item 0.8: More liberal selection for exploratory analysis
#'   }
#' @param ortho.components Number of orthogonal components to calculate (default: 1).
#'   Orthogonal components capture variation unrelated to group differences.
#'   Typically 1-3 components are sufficient for most datasets
#' @param pred.components Number of predictive components (default: 1 for binary,
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
#' @param cv.folds Number of cross-validation folds (default: 7). Only used when
#'   validation = "CV". Should be <= number of samples in smallest group
#'
#' @return A named list containing:
#'   \describe{
#'     \item{\code{model}}{Complete OPLS model object from ropls package containing:
#'       \itemize{
#'         \item Model parameters and coefficients
#'         \item Cross-validation metrics (R2Y, Q2Y)
#'         \item Model diagnostics and loadings
#'       }}
#'     \item{\code{scores}}{Data frame with sample scores (coordinates) for visualization:
#'       \itemize{
#'         \item \code{t1, t2, ...}: Predictive component scores
#'         \item \code{to1, to2, ...}: Orthogonal component scores  
#'         \item \code{sample_id}: Sample identifiers
#'         \item \code{group}: Group labels for coloring plots
#'       }}
#'     \item{\code{vip_scores}}{Data frame with Variable Importance in Projection:
#'       \itemize{
#'         \item \code{feature}: Variable/feature names
#'         \item \code{vip}: VIP score values
#'         \item \code{important}: Logical indicating if VIP >= threshold
#'       }
#'       Sorted by VIP score in descending order}
#'     \item{\code{loadings}}{Data frame with variable loadings on components}
#'     \item{\code{model_summary}}{List with model performance metrics:
#'       \itemize{
#'         \item \code{R2Y}: Fraction of Y-variance explained
#'         \item \code{Q2Y}: Cross-validated predictive ability
#'         \item \code{n_components}: Number of predictive/orthogonal components
#'         \item \code{n_variables}: Number of variables in model
#'         \item \code{n_samples}: Number of samples
#'       }}
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
#'   vip.threshold = 1.5,        # More stringent VIP cutoff
#'   scaling = "pareto",          # Pareto scaling for metabolomics
#'   cv.folds = 5                 # 5-fold cross-validation
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
#'   pred.components = 2,         # May need more components for multiple groups
#'   vip.threshold = 1.2
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
#'      cv_metrics$Value[cv_metrics$Metric == "Q2Y"]) > 0.3) {
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
opls_analysis <- function(data, group, vip.threshold = 1.0, ortho.components = 1,
                         pred.components = NULL, scaling = "standard", 
                         validation = "CV", cv.folds = 7) {
  
  # Input validation
  if (!is.matrix(data) && !is.data.frame(data) && 
      !methods::is(data, "SummarizedExperiment") && 
      !methods::is(data, "ExpressionSet")) {
    stop("'data' must be a matrix, data.frame, SummarizedExperiment, or ExpressionSet")
  }
  
  if (!is.factor(group) && !is.character(group)) {
    stop("'group' must be a factor or character vector")
  }
  
  # Convert data to matrix if needed
  if (is.data.frame(data)) {
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
  if (!is.numeric(vip.threshold) || length(vip.threshold) != 1 || vip.threshold < 0) {
    stop("'vip.threshold' must be a single non-negative number")
  }
  
  if (!is.numeric(ortho.components) || length(ortho.components) != 1 || 
      ortho.components < 0 || ortho.components != round(ortho.components)) {
    stop("'ortho.components' must be a single non-negative integer")
  }
  
  if (!is.null(pred.components)) {
    if (!is.numeric(pred.components) || length(pred.components) != 1 || 
        pred.components < 1 || pred.components != round(pred.components)) {
      stop("'pred.components' must be a single positive integer")
    }
  } else {
    # Auto-determine predictive components
    pred.components <- min(n_groups - 1, 3)  # Usually sufficient
  }
  
  valid_scaling <- c("none", "center", "standard", "pareto")
  if (!scaling %in% valid_scaling) {
    stop(paste("'scaling' must be one of:", paste(valid_scaling, collapse = ", ")))
  }
  
  valid_validation <- c("CV", "none")
  if (!validation %in% valid_validation) {
    stop("'validation' must be 'CV' or 'none'")
  }
  
  if (validation == "CV") {
    if (!is.numeric(cv.folds) || length(cv.folds) != 1 || 
        cv.folds < 2 || cv.folds != round(cv.folds)) {
      stop("'cv.folds' must be a single integer >= 2")
    }
    
    if (cv.folds > min_group_size) {
      warning(paste("cv.folds (", cv.folds, ") > smallest group size (", min_group_size, 
                    "). Reducing to", min_group_size))
      cv.folds <- min_group_size
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
  tryCatch({
    if (validation == "CV") {
      opls_model <- ropls::opls(
        x = data_matrix,
        y = group,
        predI = pred.components,
        orthoI = ortho.components,
        scaleC = scaling,
        crossvalI = cv.folds,
        fig.pdfC = "none",  # Suppress automatic plots
        info.txtC = "none"  # Suppress verbose output
      )
    } else {
      opls_model <- ropls::opls(
        x = data_matrix,
        y = group,
        predI = pred.components, 
        orthoI = ortho.components,
        scaleC = scaling,
        crossvalI = 0,  # No cross-validation
        fig.pdfC = "none",
        info.txtC = "none"
      )
    }
  }, error = function(e) {
    stop(paste("OPLS model fitting failed:", e$message,
               "\nTry reducing the number of components or checking data quality."))
  })
  
  # Extract scores for visualization
  tryCatch({
    scores_data <- chemhelper::get_scores(opls_model) %>%
      as.data.frame() %>%
      dplyr::mutate(
        sample_id = sample_names,
        group = group
      )
  }, error = function(e) {
    warning(paste("Error extracting scores with chemhelper:", e$message))
    # Fallback to manual extraction
    scores_data <- data.frame(
      t1 = opls_model@scoreMN[, 1],
      sample_id = sample_names,
      group = group
    )
    
    # Add orthogonal scores if available
    if (!is.null(opls_model@orthoScoreMN) && ncol(opls_model@orthoScoreMN) > 0) {
      scores_data$to1 <- opls_model@orthoScoreMN[, 1]
    }
  })
  
  # Extract VIP scores
  vip_data <- tryCatch({
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
        important = as.numeric(vip_values) >= vip.threshold
      ) %>%
        dplyr::arrange(desc(vip))
    }
  }, error = function(e) {
    warning(paste("Error extracting VIP scores:", e$message))
    data.frame(
      feature = variable_names,
      vip = rep(NA, length(variable_names)),
      important = rep(FALSE, length(variable_names))
    )
  })
  
  # Extract loadings
  loadings_data <- tryCatch({
    if (!is.null(opls_model@loadingMN)) {
      loading_df <- as.data.frame(opls_model@loadingMN)
      loading_df$feature <- variable_names
      loading_df
    } else {
      data.frame(feature = variable_names)
    }
  }, error = function(e) {
    warning(paste("Error extracting loadings:", e$message))
    data.frame(feature = variable_names)
  })
  
  # Extract model summary statistics
  model_summary <- tryCatch({
    summary_df <- opls_model@summaryDF
    
    list(
      R2Y = if ("R2Y" %in% names(summary_df)) summary_df$R2Y else NA,
      Q2Y = if ("Q2Y" %in% names(summary_df)) summary_df$Q2Y else NA,
      RMSEE = if ("RMSEE" %in% names(summary_df)) summary_df$RMSEE else NA,
      RMSECV = if ("RMSECV" %in% names(summary_df)) summary_df$RMSECV else NA,
      n_pred_components = pred.components,
      n_ortho_components = ortho.components,
      n_variables = ncol(data_matrix),
      n_samples = nrow(data_matrix),
      n_groups = n_groups,
      group_sizes = as.list(group_sizes),
      scaling_method = scaling,
      validation_method = validation
    )
  }, error = function(e) {
    warning(paste("Error extracting model summary:", e$message))
    list(
      R2Y = NA, Q2Y = NA, RMSEE = NA, RMSECV = NA,
      n_pred_components = pred.components,
      n_ortho_components = ortho.components,
      n_variables = ncol(data_matrix),
      n_samples = nrow(data_matrix),
      n_groups = n_groups
    )
  })
  
  # Count important variables
  n_important_vars <- sum(vip_data$important, na.rm = TRUE)
  
  # Prepare final results
  results <- list(
    model = opls_model,
    scores = scores_data,
    vip_scores = vip_data,
    loadings = loadings_data,
    model_summary = model_summary
  )
  
  # Add metadata as attributes
  attr(results, "analysis_type") <- "OPLS-DA"
  attr(results, "n_samples") <- nrow(data_matrix)
  attr(results, "n_variables") <- ncol(data_matrix)
  attr(results, "n_groups") <- n_groups
  attr(results, "n_important_vars") <- n_important_vars
  attr(results, "vip_threshold") <- vip.threshold
  attr(results, "scaling_method") <- scaling
  
  # Print summary if interactive
  if (interactive()) {
    cat("OPLS-DA Analysis Summary:\n")
    cat("========================\n")
    cat("Samples:", nrow(data_matrix), "\n")
    cat("Variables:", ncol(data_matrix), "\n")
    cat("Groups:", n_groups, "(", paste(names(group_sizes), collapse = ", "), ")\n")
    cat("Group sizes:", paste(group_sizes, collapse = ", "), "\n")
    cat("Predictive components:", pred.components, "\n")
    cat("Orthogonal components:", ortho.components, "\n")
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
    cat("VIP threshold:", vip.threshold, "\n")
    cat("Important variables:", n_important_vars, 
        paste0("(", round(100 * n_important_vars / ncol(data_matrix), 1), "%)"), "\n")
    
    if (n_important_vars > 0) {
      cat("\nTop 5 important variables:\n")
      top_vip <- head(vip_data[vip_data$important, ], 5)
      for (i in 1:nrow(top_vip)) {
        cat(sprintf("  %s: VIP = %.3f\n", top_vip$feature[i], top_vip$vip[i]))
      }
    }
    cat("\n")
  }
  
  return(results)
}