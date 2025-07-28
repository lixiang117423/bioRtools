#' Perform Principal Component Analysis (PCA) with Comprehensive Visualization
#'
#' This function performs Principal Component Analysis (PCA) on high-dimensional
#' data and provides comprehensive results including eigenvalue analysis, sample
#' scores, variable loadings, and multiple visualization options. PCA is an
#' unsupervised dimensionality reduction technique that identifies the directions
#' of maximum variance in the data, making it essential for exploratory data
#' analysis in genomics, metabolomics, and other omics fields.
#'
#' @param data Numeric data frame or matrix where rows represent samples (observations)
#'   and columns represent variables (features/genes/metabolites). Missing values
#'   are allowed and will be handled through imputation or exclusion. All variables
#'   should be numeric; categorical variables should be excluded or encoded separately
#' @param sample Data frame containing sample metadata with at least one column
#'   matching the row names of the data matrix. Must include:
#'   \itemize{
#'     \item Sample identifiers that match data row names
#'     \item Grouping variables for visualization (e.g., treatment, condition, species)
#'     \item Optional covariates for interpretation
#'   }
#' @param n.components Maximum number of principal components to calculate and
#'   return (default: 5). Should be <= min(n_samples-1, n_variables). More
#'   components provide more detailed analysis but increase computation time
#' @param scale.data Logical indicating whether to scale variables to unit variance
#'   (default: TRUE). Recommended when variables have different units or scales.
#'   Set to FALSE if data is already standardized or when preserving original scales
#' @param center.data Logical indicating whether to center variables to zero mean
#'   (default: TRUE). Almost always recommended for PCA analysis
#' @param x.axis Principal component for X-axis in score plot (default: "PC1").
#'   Can be "PC1", "PC2", "PC3", etc., up to the number of calculated components
#' @param y.axis Principal component for Y-axis in score plot (default: "PC2").
#'   Should be different from x.axis for meaningful 2D visualization
#' @param color.by Column name in sample metadata for point colors (default: "group").
#'   Used to visualize sample groupings or experimental conditions
#' @param shape.by Column name in sample metadata for point shapes (default: same as color.by).
#'   Can be the same as color.by or different to show additional sample attributes
#' @param plot.type Type of plots to generate (default: "all"). Options:
#'   \itemize{
#'     \item "all": Generate all available plots
#'     \item "scores": Only sample score plots
#'     \item "loadings": Only variable loading plots
#'     \item "variance": Only variance explanation plots
#'     \item "none": No plots (results only)
#'   }
#' @param conf.ellipses Logical indicating whether to add confidence ellipses
#'   around groups in score plots (default: FALSE). Useful for visualizing
#'   group separation and overlap
#' @param ellipse.level Confidence level for ellipses (default: 0.95). Only
#'   used when conf.ellipses = TRUE
#'
#' @return A named list containing comprehensive PCA results:
#'   \describe{
#'     \item{\code{pca_model}}{Complete PCA model object from FactoMineR package}
#'     \item{\code{eigenvalues}}{Data frame with eigenvalue analysis:
#'       \itemize{
#'         \item \code{component}: Principal component names (PC1, PC2, ...)
#'         \item \code{eigenvalue}: Eigenvalue for each component
#'         \item \code{variance_percent}: Percentage of variance explained
#'         \item \code{cumulative_percent}: Cumulative variance explained
#'       }}
#'     \item{\code{sample_scores}}{Data frame with sample coordinates on principal components:
#'       \itemize{
#'         \item \code{sample_id}: Sample identifiers
#'         \item \code{PC1, PC2, ...}: Principal component scores
#'         \item All columns from sample metadata for grouping/coloring
#'       }}
#'     \item{\code{variable_loadings}}{Data frame with variable contributions:
#'       \itemize{
#'         \item \code{variable}: Variable names
#'         \item \code{PC1, PC2, ...}: Loading values on each component
#'         \item \code{cos2_PC1, cos2_PC2, ...}: Quality of representation
#'       }}
#'     \item{\code{plots}}{List of ggplot2 objects (when plot.type != "none"):
#'       \itemize{
#'         \item \code{score_plot}: Sample scores with grouping
#'         \item \code{scree_plot}: Variance explained by each component
#'         \item \code{loading_plot}: Variable loadings visualization
#'         \item \code{biplot}: Combined scores and loadings plot
#'         \item \code{contribution_plot}: Variable contributions to components
#'       }}
#'     \item{\code{summary_stats}}{List with key summary statistics:
#'       \itemize{
#'         \item \code{total_variance_explained}: Sum of first n components
#'         \item \code{kaiser_criterion}: Components with eigenvalue > 1
#'         \item \code{broken_stick}: Components above broken stick model
#'       }}
#'   }
#'
#' @details
#' \strong{PCA Workflow:}
#' \enumerate{
#'   \item Data preprocessing (centering and optional scaling)
#'   \item Covariance/correlation matrix calculation
#'   \item Eigenvalue decomposition to find principal components
#'   \item Sample projection onto component space
#'   \item Variable loading calculation
#'   \item Variance explanation analysis
#' }
#' 
#' \strong{Interpretation Guidelines:}
#' \itemize{
#'   \item \strong{Eigenvalues}: Components with eigenvalue > 1 are typically retained
#'   \item \strong{Variance explained}: First 2-3 components usually explain 50-80% of variance
#'   \item \strong{Sample scores}: Show relationships between samples in reduced space
#'   \item \strong{Variable loadings}: Indicate which variables contribute to each component
#'   \item \strong{Cos2 values}: Quality of variable representation (0-1, higher is better)
#' }
#' 
#' \strong{Data Requirements:}
#' \itemize{
#'   \item Numeric variables only (continuous or discrete)
#'   \item Sample size >> number of variables (for stable results)
#'   \item Consider removing low-variance variables
#'   \item Handle missing values appropriately
#' }
#'
#' @note
#' \itemize{
#'   \item PCA assumes linear relationships between variables
#'   \item Results are sensitive to outliers; consider robust PCA for outlier-prone data
#'   \item Scaling is crucial when variables have different units or ranges
#'   \item For large datasets, consider approximate PCA methods for computational efficiency
#'   \item Interpretation requires domain knowledge of the variables and samples
#' }
#'
#' @references
#' Jolliffe, I.T. and Cadima, J. (2016). Principal component analysis: a review and
#' recent developments. Philosophical Transactions of the Royal Society A, 374(2065), 20150202.
#' 
#' Husson, F., et al. (2010). FactoMineR: An R Package for Multivariate Exploratory
#' Data Analysis and Data Mining. Journal of Statistical Software, 25(1), 1-18.
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{opls_analysis}} for supervised multivariate analysis
#'   \item \code{\link{cor_analysis}} for correlation-based analysis
#'   \item \code{\link[FactoMineR]{PCA}} for underlying PCA implementation
#'   \item \code{\link[factoextra]{fviz_pca_ind}} for advanced PCA visualization
#' }
#'
#' @export
#'
#' @examples
#' library(bioRtools)
#' library(dplyr)
#' library(ggplot2)
#' 
#' # Example 1: Basic PCA analysis with iris data
#' # Prepare data and sample information
#' iris_data <- iris[, 1:4]  # Numeric variables only
#' iris_samples <- data.frame(
#'   sample_id = paste0("sample_", 1:nrow(iris)),
#'   species = iris$Species,
#'   sepal_size = cut(iris$Sepal.Length, breaks = 3, labels = c("Small", "Medium", "Large"))
#' )
#' rownames(iris_data) <- iris_samples$sample_id
#' 
#' # Perform PCA analysis
#' pca_results <- pca_analysis(
#'   data = iris_data,
#'   sample = iris_samples,
#'   color.by = "species"
#' )
#' 
#' # View variance explained
#' print("Variance explained by each component:")
#' print(pca_results$eigenvalues)
#' 
#' # Display score plot
#' print(pca_results$plots$score_plot)
#' 
#' # Check summary statistics
#' print(paste("Total variance in first 2 PCs:", 
#'             round(pca_results$summary_stats$total_variance_explained[2], 1), "%"))
#' 
#' # Example 2: Customized PCA with different visualization options
#' pca_custom <- pca_analysis(
#'   data = iris_data,
#'   sample = iris_samples,
#'   n.components = 3,           # Calculate 3 components
#'   x.axis = "PC1",             # X-axis component
#'   y.axis = "PC3",             # Y-axis component (skip PC2)
#'   color.by = "species",       # Color by species
#'   shape.by = "sepal_size",    # Shape by sepal size
#'   conf.ellipses = TRUE,       # Add confidence ellipses
#'   ellipse.level = 0.95        # 95% confidence ellipses
#' )
#' 
#' print("Customized PCA plot (PC1 vs PC3):")
#' print(pca_custom$plots$score_plot)
#' 
#' # Example 3: Metabolomics-style analysis with scaling options
#' \dontrun{
#' # Simulate metabolomics data
#' set.seed(123)
#' metabolite_data <- matrix(
#'   rlnorm(50 * 20, meanlog = 5, sdlog = 1),  # Log-normal distribution typical of metabolomics
#'   nrow = 50, ncol = 20,
#'   dimnames = list(
#'     paste0("Sample_", 1:50),
#'     paste0("Metabolite_", 1:20)
#'   )
#' )
#' 
#' metabolite_samples <- data.frame(
#'   sample_id = rownames(metabolite_data),
#'   treatment = rep(c("Control", "Treatment_A", "Treatment_B"), length.out = 50),
#'   batch = rep(c("Batch1", "Batch2"), length.out = 50),
#'   timepoint = rep(c("T0", "T1", "T2", "T3", "T4"), each = 10)
#' )
#' 
#' # PCA with scaling (important for metabolomics)
#' metabolomics_pca <- pca_analysis(
#'   data = metabolite_data,
#'   sample = metabolite_samples,
#'   scale.data = TRUE,          # Scale for different metabolite ranges
#'   color.by = "treatment",
#'   shape.by = "timepoint",
#'   plot.type = "all"
#' )
#' 
#' # Check for batch effects
#' batch_plot <- metabolomics_pca$sample_scores %>%
#'   ggplot(aes(x = PC1, y = PC2, color = batch, shape = treatment)) +
#'   geom_point(size = 3) +
#'   labs(title = "PCA: Batch Effect Assessment") +
#'   theme_minimal()
#' 
#' print(batch_plot)
#' }
#' 
#' # Example 4: Variable importance and loading analysis
#' # Examine which variables contribute most to PC1 and PC2
#' variable_importance <- pca_results$variable_loadings %>%
#'   mutate(
#'     PC1_importance = abs(PC1),
#'     PC2_importance = abs(PC2),
#'     total_importance = sqrt(PC1^2 + PC2^2)
#'   ) %>%
#'   arrange(desc(total_importance))
#' 
#' print("Variables most important for PC1-PC2 space:")
#' print(head(variable_importance, 10))
#' 
#' # Create custom loading plot
#' loading_plot <- variable_importance %>%
#'   ggplot(aes(x = PC1, y = PC2)) +
#'   geom_segment(aes(xend = 0, yend = 0), 
#'                arrow = arrow(length = unit(0.3, "cm"))) +
#'   geom_text(aes(label = variable), vjust = -1) +
#'   labs(
#'     title = "PCA Loading Plot",
#'     subtitle = "Variable contributions to PC1 and PC2",
#'     x = "PC1 Loading",
#'     y = "PC2 Loading"
#'   ) +
#'   theme_minimal() +
#'   coord_equal()
#' 
#' print(loading_plot)
#' 
#' # Example 5: Quality assessment and component selection
#' # Scree plot analysis for component selection
#' scree_analysis <- pca_results$eigenvalues %>%
#'   mutate(
#'     above_kaiser = eigenvalue > 1,  # Kaiser criterion
#'     component_num = row_number()
#'   )
#' 
#' print("Component selection criteria:")
#' print(scree_analysis)
#' 
#' # Enhanced scree plot with selection criteria
#' enhanced_scree <- scree_analysis %>%
#'   ggplot(aes(x = component_num, y = variance_percent)) +
#'   geom_line(color = "blue", size = 1) +
#'   geom_point(aes(color = above_kaiser), size = 3) +
#'   geom_hline(yintercept = 100/ncol(iris_data), 
#'              linetype = "dashed", color = "red",
#'              alpha = 0.7) +
#'   labs(
#'     title = "Scree Plot with Selection Criteria",
#'     subtitle = "Red line: Average variance per component",
#'     x = "Principal Component",
#'     y = "Variance Explained (%)",
#'     color = "Above Kaiser\nCriterion"
#'   ) +
#'   theme_minimal()
#' 
#' print(enhanced_scree)
#' 
#' # Example 6: Sample quality assessment using PCA
#' # Identify potential outliers based on PCA scores
#' sample_distances <- pca_results$sample_scores %>%
#'   mutate(
#'     pc_distance = sqrt(PC1^2 + PC2^2),  # Distance from origin
#'     is_outlier = pc_distance > quantile(pc_distance, 0.95)  # Top 5% as outliers
#'   )
#' 
#' outlier_plot <- sample_distances %>%
#'   ggplot(aes(x = PC1, y = PC2, color = species, shape = is_outlier)) +
#'   geom_point(size = 3) +
#'   scale_shape_manual(values = c(16, 17), name = "Potential\nOutlier") +
#'   labs(title = "PCA Outlier Detection") +
#'   theme_minimal()
#' 
#' print(outlier_plot)
#' 
#' # Summary of outliers
#' outlier_summary <- sample_distances %>%
#'   filter(is_outlier) %>%
#'   select(sample_id, species, PC1, PC2, pc_distance)
#' 
#' if (nrow(outlier_summary) > 0) {
#'   print("Potential outlier samples:")
#'   print(outlier_summary)
#' }
#' 
#' # Example 7: Comparative PCA (with and without scaling)
#' \dontrun{
#' # Compare scaled vs unscaled PCA
#' pca_unscaled <- pca_analysis(
#'   data = iris_data,
#'   sample = iris_samples,
#'   scale.data = FALSE,
#'   color.by = "species",
#'   plot.type = "scores"
#' )
#' 
#' pca_scaled <- pca_analysis(
#'   data = iris_data,
#'   sample = iris_samples,
#'   scale.data = TRUE,
#'   color.by = "species",
#'   plot.type = "scores"
#' )
#' 
#' # Compare variance explained
#' scaling_comparison <- data.frame(
#'   Component = 1:2,
#'   Unscaled_PC1_PC2 = pca_unscaled$eigenvalues$variance_percent[1:2],
#'   Scaled_PC1_PC2 = pca_scaled$eigenvalues$variance_percent[1:2]
#' )
#' 
#' print("Scaling effect on variance explained:")
#' print(scaling_comparison)
#' 
#' # Side-by-side plot comparison
#' library(gridExtra)
#' grid.arrange(
#'   pca_unscaled$plots$score_plot + ggtitle("Unscaled PCA"),
#'   pca_scaled$plots$score_plot + ggtitle("Scaled PCA"),
#'   ncol = 2
#' )
#' }
#'
pca_analysis <- function(data, sample, n.components = 5, scale.data = TRUE, 
                        center.data = TRUE, x.axis = "PC1", y.axis = "PC2",
                        color.by = "group", shape.by = NULL, plot.type = "all",
                        conf.ellipses = FALSE, ellipse.level = 0.95) {
  
  # Input validation
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("'data' must be a data frame or matrix")
  }
  
  if (!is.data.frame(sample)) {
    stop("'sample' must be a data frame")
  }
  
  if (nrow(data) == 0 || ncol(data) == 0) {
    stop("'data' cannot be empty")
  }
  
  # Convert data to numeric matrix
  data_matrix <- as.matrix(data)
  if (!is.numeric(data_matrix)) {
    stop("'data' must contain only numeric values")
  }
  
  # Check for sample matching
  data_samples <- rownames(data)
  if (is.null(data_samples)) {
    data_samples <- paste0("Sample_", 1:nrow(data))
    rownames(data) <- data_samples
  }
  
  # Validate sample metadata
  if (!"sample_id" %in% names(sample)) {
    # Try to find a sample identifier column
    sample_id_cols <- c("sample", "Sample", "ID", "id", "sample_name")
    found_col <- intersect(sample_id_cols, names(sample))
    
    if (length(found_col) > 0) {
      sample$sample_id <- sample[[found_col[1]]]
    } else {
      stop("Sample metadata must contain a 'sample_id' column or similar identifier")
    }
  }
  
  # Check sample matching between data and metadata
  if (!all(data_samples %in% sample$sample_id)) {
    missing_samples <- setdiff(data_samples, sample$sample_id)
    stop(paste("Samples missing from metadata:", paste(head(missing_samples, 5), collapse = ", ")))
  }
  
  # Align sample metadata with data
  sample_aligned <- sample[match(data_samples, sample$sample_id), ]
  
  # Validate parameters
  if (!is.numeric(n.components) || length(n.components) != 1 || 
      n.components < 1 || n.components != round(n.components)) {
    stop("'n.components' must be a positive integer")
  }
  
  max_components <- min(nrow(data) - 1, ncol(data))
  if (n.components > max_components) {
    warning(paste("n.components reduced from", n.components, "to", max_components))
    n.components <- max_components
  }
  
  if (!is.logical(scale.data) || length(scale.data) != 1) {
    stop("'scale.data' must be a single logical value")
  }
  
  if (!is.logical(center.data) || length(center.data) != 1) {
    stop("'center.data' must be a single logical value")
  }
  
  # Validate axis specifications
  valid_axes <- paste0("PC", 1:n.components)
  if (!x.axis %in% valid_axes) {
    stop(paste("'x.axis' must be one of:", paste(valid_axes, collapse = ", ")))
  }
  
  if (!y.axis %in% valid_axes) {
    stop(paste("'y.axis' must be one of:", paste(valid_axes, collapse = ", ")))
  }
  
  if (x.axis == y.axis) {
    stop("'x.axis' and 'y.axis' must be different")
  }
  
  # Validate grouping variables
  if (!color.by %in% names(sample_aligned)) {
    stop(paste("'color.by' column '", color.by, "' not found in sample metadata"))
  }
  
  if (is.null(shape.by)) {
    shape.by <- color.by
  } else if (!shape.by %in% names(sample_aligned)) {
    stop(paste("'shape.by' column '", shape.by, "' not found in sample metadata"))
  }
  
  valid_plot_types <- c("all", "scores", "loadings", "variance", "none")
  if (!plot.type %in% valid_plot_types) {
    stop(paste("'plot.type' must be one of:", paste(valid_plot_types, collapse = ", ")))
  }
  
  # Data quality checks
  if (any(is.na(data_matrix))) {
    na_count <- sum(is.na(data_matrix))
    warning(paste("Data contains", na_count, "missing values. PCA will handle these via imputation."))
  }
  
  # Check for zero variance variables
  var_check <- apply(data_matrix, 2, var, na.rm = TRUE)
  zero_var_vars <- sum(var_check == 0, na.rm = TRUE)
  
  if (zero_var_vars > 0) {
    warning(paste("Found", zero_var_vars, "variables with zero variance. Consider removing these."))
  }
  
  # Perform PCA analysis
  tryCatch({
    pca_model <- FactoMineR::PCA(
      X = data_matrix,
      ncp = n.components,
      scale.unit = scale.data,
      graph = FALSE,
      ind.sup = NULL,
      quanti.sup = NULL
    )
  }, error = function(e) {
    stop(paste("PCA analysis failed:", e$message))
  })
  
  # Extract eigenvalues and variance explanation
  eigenvalue_data <- tryCatch({
    eig_values <- factoextra::get_eigenvalue(pca_model)
    
    data.frame(
      component = paste0("PC", 1:nrow(eig_values)),
      eigenvalue = eig_values$eigenvalue,
      variance_percent = round(eig_values$variance.percent, 2),
      cumulative_percent = round(eig_values$cumulative.variance.percent, 2)
    )
  }, error = function(e) {
    warning(paste("Error extracting eigenvalues:", e$message))
    data.frame(
      component = paste0("PC", 1:n.components),
      eigenvalue = rep(NA, n.components),
      variance_percent = rep(NA, n.components),
      cumulative_percent = rep(NA, n.components)
    )
  })
  
  # Extract sample scores
  sample_scores <- tryCatch({
    scores_matrix <- pca_model$ind$coord
    scores_df <- as.data.frame(scores_matrix)
    names(scores_df) <- paste0("PC", 1:ncol(scores_df))
    
    scores_df$sample_id <- data_samples
    scores_df <- scores_df %>%
      dplyr::left_join(sample_aligned, by = "sample_id")
    
    scores_df
  }, error = function(e) {
    warning(paste("Error extracting sample scores:", e$message))
    data.frame(sample_id = data_samples)
  })
  
  # Extract variable loadings
  variable_loadings <- tryCatch({
    loadings_matrix <- pca_model$var$coord
    cos2_matrix <- pca_model$var$cos2
    
    loadings_df <- as.data.frame(loadings_matrix)
    cos2_df <- as.data.frame(cos2_matrix)
    
    names(loadings_df) <- paste0("PC", 1:ncol(loadings_df))
    names(cos2_df) <- paste0("cos2_PC", 1:ncol(cos2_df))
    
    result_df <- data.frame(
      variable = rownames(loadings_matrix),
      loadings_df,
      cos2_df
    )
    
    result_df
  }, error = function(e) {
    warning(paste("Error extracting variable loadings:", e$message))
    data.frame(variable = colnames(data))
  })
  
  # Generate plots based on plot.type
  plots <- list()
  
  if (plot.type %in% c("all", "scores")) {
    # Score plot
    x_pc_num <- as.numeric(gsub("PC", "", x.axis))
    y_pc_num <- as.numeric(gsub("PC", "", y.axis))
    
    x_variance <- eigenvalue_data$variance_percent[x_pc_num]
    y_variance <- eigenvalue_data$variance_percent[y_pc_num]
    
    score_plot <- sample_scores %>%
      ggplot2::ggplot(ggplot2::aes(x = !!rlang::sym(x.axis), y = !!rlang::sym(y.axis),
                                   color = !!rlang::sym(color.by), 
                                   shape = !!rlang::sym(shape.by))) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", alpha = 0.7) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", alpha = 0.7) +
      ggplot2::geom_point(size = 3, alpha = 0.8) +
      ggplot2::labs(
        x = paste0(x.axis, " (", x_variance, "%)"),
        y = paste0(y.axis, " (", y_variance, "%)"),
        title = "PCA Score Plot",
        subtitle = paste("First", n.components, "components explain", 
                        round(eigenvalue_data$cumulative_percent[n.components], 1), "% of variance"),
        color = tools::toTitleCase(gsub("[._]", " ", color.by)),
        shape = tools::toTitleCase(gsub("[._]", " ", shape.by))
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 14, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 11),
        axis.title = ggplot2::element_text(size = 12),
        legend.title = ggplot2::element_text(size = 11)
      )
    
    # Add confidence ellipses if requested
    if (conf.ellipses) {
      score_plot <- score_plot +
        ggplot2::stat_ellipse(ggplot2::aes(color = !!rlang::sym(color.by)), 
                             level = ellipse.level, type = "norm", alpha = 0.3)
    }
    
    plots$score_plot <- score_plot
  }
  
  if (plot.type %in% c("all", "variance")) {
    # Scree plot
    scree_plot <- eigenvalue_data %>%
      ggplot2::ggplot(ggplot2::aes(x = factor(1:nrow(eigenvalue_data)), y = variance_percent)) +
      ggplot2::geom_line(group = 1, color = "blue", size = 1) +
      ggplot2::geom_point(color = "blue", size = 3) +
      ggplot2::geom_hline(yintercept = 100/ncol(data), linetype = "dashed", 
                         color = "red", alpha = 0.7) +
      ggplot2::labs(
        x = "Principal Component",
        y = "Variance Explained (%)",
        title = "Scree Plot",
        subtitle = "Red line: Average variance per variable"
      ) +
      ggplot2::theme_minimal()
    
    plots$scree_plot <- scree_plot
  }
  
  if (plot.type %in% c("all", "loadings")) {
    # Loading plot
    x_loadings <- paste0("PC", x_pc_num)
    y_loadings <- paste0("PC", y_pc_num)
    
    loading_plot <- variable_loadings %>%
      ggplot2::ggplot(ggplot2::aes(x = !!rlang::sym(x_loadings), y = !!rlang::sym(y_loadings))) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
      ggplot2::geom_segment(ggplot2::aes(xend = 0, yend = 0), 
                           arrow = ggplot2::arrow(length = ggplot2::unit(0.3, "cm")),
                           alpha = 0.7) +
      ggplot2::geom_text(ggplot2::aes(label = variable), 
                        vjust = -0.5, hjust = 0.5, size = 3, alpha = 0.8) +
      ggplot2::labs(
        x = paste0(x.axis, " Loadings"),
        y = paste0(y.axis, " Loadings"),
        title = "PCA Loading Plot",
        subtitle = "Variable contributions to principal components"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::coord_equal()
    
    plots$loading_plot <- loading_plot
  }
  
  if (plot.type == "all") {
    # Biplot (combined scores and loadings)
    # Scale loadings for better visualization
    loading_scale <- 3
    scaled_loadings <- variable_loadings
    scaled_loadings[[x_loadings]] <- scaled_loadings[[x_loadings]] * loading_scale
    scaled_loadings[[y_loadings]] <- scaled_loadings[[y_loadings]] * loading_scale
    
    biplot <- sample_scores %>%
      ggplot2::ggplot(ggplot2::aes(x = !!rlang::sym(x.axis), y = !!rlang::sym(y.axis))) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", alpha = 0.7) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", alpha = 0.7) +
      ggplot2::geom_point(ggplot2::aes(color = !!rlang::sym(color.by)), size = 2, alpha = 0.7) +
      ggplot2::geom_segment(data = scaled_loadings,
                           ggplot2::aes(x = 0, y = 0, 
                                       xend = !!rlang::sym(x_loadings), 
                                       yend = !!rlang::sym(y_loadings)),
                           arrow = ggplot2::arrow(length = ggplot2::unit(0.2, "cm")),
                           color = "red", alpha = 0.6) +
      ggplot2::geom_text(data = scaled_loadings,
                        ggplot2::aes(x = !!rlang::sym(x_loadings), 
                                    y = !!rlang::sym(y_loadings),
                                    label = variable),
                        color = "red", size = 2.5, alpha = 0.8) +
      ggplot2::labs(
        x = paste0(x.axis, " (", x_variance, "%)"),
        y = paste0(y.axis, " (", y_variance, "%)"),
        title = "PCA Biplot",
        subtitle = "Samples (points) and variables (arrows)",
        color = tools::toTitleCase(gsub("[._]", " ", color.by))
      ) +
      ggplot2::theme_minimal()
    
    plots$biplot <- biplot
    
    # Contribution plot
    contrib_data <- variable_loadings %>%
      dplyr::select(variable, dplyr::all_of(paste0("cos2_PC", 1:min(3, n.components)))) %>%
      tidyr::pivot_longer(cols = -variable, names_to = "component", values_to = "contribution") %>%
      dplyr::mutate(component = gsub("cos2_", "", component))
    
    contribution_plot <- contrib_data %>%
      ggplot2::ggplot(ggplot2::aes(x = reorder(variable, contribution), 
                                   y = contribution, fill = component)) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::coord_flip() +
      ggplot2::labs(
        x = "Variables",
        y = "Quality of Representation (cos²)",
        title = "Variable Contribution to Principal Components",
        fill = "Component"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
    
    plots$contribution_plot <- contribution_plot
  }
  
  # Calculate summary statistics
  summary_stats <- list(
    total_variance_explained = eigenvalue_data$cumulative_percent,
    kaiser_criterion = eigenvalue_data$component[eigenvalue_data$eigenvalue > 1],
    broken_stick = {
      # Broken stick model for component selection
      n_vars <- ncol(data)
      broken_stick_values <- sapply(1:n.components, function(i) {
        100 * sum(1/((i:n_vars))) / n_vars
      })
      eigenvalue_data$component[eigenvalue_data$variance_percent > broken_stick_values]
    },
    average_variance_per_variable = 100 / ncol(data),
    n_components_kaiser = sum(eigenvalue_data$eigenvalue > 1, na.rm = TRUE),
    effective_dimensionality = sum(eigenvalue_data$eigenvalue^2, na.rm = TRUE) / sum(eigenvalue_data$eigenvalue, na.rm = TRUE)^2
  )
  
  # Prepare final results
  results <- list(
    pca_model = pca_model,
    eigenvalues = eigenvalue_data,
    sample_scores = sample_scores,
    variable_loadings = variable_loadings,
    plots = if (plot.type != "none") plots else NULL,
    summary_stats = summary_stats
  )
  
  # Add metadata as attributes
  attr(results, "analysis_type") <- "PCA"
  attr(results, "n_samples") <- nrow(data)
  attr(results, "n_variables") <- ncol(data)
  attr(results, "n_components") <- n.components
  attr(results, "scaled") <- scale.data
  attr(results, "centered") <- center.data
  
  # Interactive summary
  if (interactive()) {
    cat("Principal Component Analysis Summary:\n")
    cat("====================================\n")
    cat("Samples:", nrow(data), "\n")
    cat("Variables:", ncol(data), "\n")
    cat("Components calculated:", n.components, "\n")
    cat("Data scaling:", ifelse(scale.data, "Applied", "Not applied"), "\n")
    cat("Data centering:", ifelse(center.data, "Applied", "Not applied"), "\n\n")
    
    cat("Variance Explanation:\n")
    for (i in 1:min(3, nrow(eigenvalue_data))) {
      cat(sprintf("  %s: %.1f%% (eigenvalue = %.2f)\n",
                  eigenvalue_data$component[i],
                  eigenvalue_data$variance_percent[i],
                  eigenvalue_data$eigenvalue[i]))
    }
    
    cat(sprintf("  First %d components: %.1f%% total variance\n",
                min(3, n.components),
                eigenvalue_data$cumulative_percent[min(3, n.components)]))
    
    cat("\nComponent Selection Criteria:\n")
    cat("  Kaiser criterion (eigenvalue > 1):", length(summary_stats$kaiser_criterion), "components\n")
    cat("  Broken stick model suggests:", length(summary_stats$broken_stick), "components\n")
    cat("  Average variance per variable:", round(summary_stats$average_variance_per_variable, 1), "%\n")
    
    # Variable importance
    if (nrow(variable_loadings) > 0) {
      important_vars_pc1 <- variable_loadings %>%
        dplyr::arrange(desc(abs(PC1))) %>%
        head(3)
      
      cat("\nTop 3 variables contributing to PC1:\n")
      for (i in 1:nrow(important_vars_pc1)) {
        cat(sprintf("  %s (loading = %.3f)\n",
                    important_vars_pc1$variable[i],
                    important_vars_pc1$PC1[i]))
      }
    }
    
    # Data quality assessment
    if (any(is.na(data_matrix))) {
      cat("\nData Quality Notes:\n")
      cat("  Missing values detected - handled via imputation\n")
    }
    
    if (zero_var_vars > 0) {
      cat("  Zero variance variables detected - consider removal\n")
    }
    
    cat("\n")
  }
  
  return(results)
}