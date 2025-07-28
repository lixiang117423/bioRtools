#' Calculate Correlation Coefficients and Statistical Significance
#'
#' This function calculates pairwise correlations between variables within a single
#' dataset or between variables from two different datasets. It provides correlation
#' coefficients, p-values, and allows filtering based on correlation strength and
#' statistical significance thresholds.
#'
#' @param data.1 First data frame where each column represents a feature/variable.
#'   All columns should be numeric.
#' @param data.2 Optional second data frame for cross-dataset correlation analysis.
#'   If NULL (default), correlations are calculated within data.1. Each column
#'   should represent a feature/variable and be numeric.
#' @param method Correlation method to use. Options are:
#'   \itemize{
#'     \item \code{"pearson"} (default): Pearson product-moment correlation
#'     \item \code{"kendall"}: Kendall's tau correlation (non-parametric)
#'     \item \code{"spearman"}: Spearman's rank correlation (non-parametric)
#'   }
#' @param cor Minimum absolute correlation coefficient threshold for filtering
#'   results (default: 0.6). Only correlations with |r| >= cor will be returned.
#' @param pvalue Maximum p-value threshold for statistical significance
#'   (default: 0.05). Only correlations with p-value <= pvalue will be returned.
#'
#' @return A data frame containing significant correlations with columns:
#'   \itemize{
#'     \item \code{from}: Variable name from the first dataset
#'     \item \code{to}: Variable name from the second dataset (or same dataset)
#'     \item \code{cor}: Correlation coefficient
#'     \item \code{pvalue}: Statistical significance (p-value)
#'   }
#'   
#' @details
#' The function uses WGCNA::corAndPvalue() for efficient correlation calculation
#' with significance testing. When data.2 is NULL, it performs within-dataset
#' correlations and excludes self-correlations (from != to). When data.2 is
#' provided, it performs between-dataset correlations.
#' 
#' The results are filtered to include only correlations that meet both the
#' correlation strength threshold (|cor| >= cor) and significance threshold
#' (p-value <= pvalue).
#'
#' @note
#' \itemize{
#'   \item All input data should be numeric
#'   \item Missing values are handled by the underlying correlation function
#'   \item For large datasets, consider adjusting p-value thresholds for multiple testing
#' }
#'
#' @export
#'
#' @examples
#' library(bioRtools)
#' library(dplyr)
#' 
#' # Example 1: Within-dataset correlation (iris features)
#' iris_numeric <- iris[, 1:4]
#' 
#' # Basic correlation analysis
#' cor_result1 <- cor_analysis(data.1 = iris_numeric)
#' print(cor_result1)
#' 
#' # Example 2: Custom thresholds
#' cor_result2 <- cor_analysis(
#'   data.1 = iris_numeric,
#'   method = "spearman",
#'   cor = 0.7,        # Higher correlation threshold
#'   pvalue = 0.01     # More stringent significance
#' )
#' print(cor_result2)
#' 
#' # Example 3: Between-dataset correlation
#' # Create two subsets of iris data
#' morphology <- iris[, 1:2]  # Sepal measurements  
#' anatomy <- iris[, 3:4]     # Petal measurements
#' 
#' cor_result3 <- cor_analysis(
#'   data.1 = morphology,
#'   data.2 = anatomy,
#'   method = "pearson"
#' )
#' print(cor_result3)
#' 
#' # Example 4: Using different correlation methods
#' # Pearson (for linear relationships)
#' pearson_result <- cor_analysis(iris_numeric, method = "pearson")
#' 
#' # Spearman (for monotonic relationships)
#' spearman_result <- cor_analysis(iris_numeric, method = "spearman")
#' 
#' # Kendall (robust to outliers)
#' kendall_result <- cor_analysis(iris_numeric, method = "kendall")
#' 
#' # Compare results
#' print("Pearson correlations:")
#' print(pearson_result)
#' print("Spearman correlations:")
#' print(spearman_result)
#'
cor_analysis <- function(data.1, data.2 = NULL, method = "pearson", cor = 0.6, pvalue = 0.05) {
  
  # Input validation
  if (!is.data.frame(data.1)) {
    stop("'data.1' must be a data frame")
  }
  
  if (!is.null(data.2) && !is.data.frame(data.2)) {
    stop("'data.2' must be a data frame or NULL")
  }
  
  if (!method %in% c("pearson", "kendall", "spearman")) {
    stop("'method' must be one of: 'pearson', 'kendall', 'spearman'")
  }
  
  if (!is.numeric(cor) || cor < 0 || cor > 1) {
    stop("'cor' must be a numeric value between 0 and 1")
  }
  
  if (!is.numeric(pvalue) || pvalue <= 0 || pvalue > 1) {
    stop("'pvalue' must be a numeric value between 0 and 1")
  }
  
  # Check if data.1 has numeric columns
  numeric_cols_1 <- sapply(data.1, is.numeric)
  if (!any(numeric_cols_1)) {
    stop("'data.1' must contain at least one numeric column")
  }
  
  # Filter to numeric columns only and warn if non-numeric columns removed
  if (!all(numeric_cols_1)) {
    data.1 <- data.1[, numeric_cols_1, drop = FALSE]
    warning("Non-numeric columns removed from data.1")
  }
  
  # Check data.2 if provided
  if (!is.null(data.2)) {
    numeric_cols_2 <- sapply(data.2, is.numeric)
    if (!any(numeric_cols_2)) {
      stop("'data.2' must contain at least one numeric column")
    }
    if (!all(numeric_cols_2)) {
      data.2 <- data.2[, numeric_cols_2, drop = FALSE]
      warning("Non-numeric columns removed from data.2")
    }
  }
  
  # Check for minimum number of observations
  if (nrow(data.1) < 3) {
    stop("'data.1' must have at least 3 observations for correlation analysis")
  }
  
  if (!is.null(data.2) && nrow(data.2) < 3) {
    stop("'data.2' must have at least 3 observations for correlation analysis")
  }
  
  # Calculate correlations based on whether data.2 is provided
  if (is.null(data.2)) {
    # Within-dataset correlation
    cor_result <- WGCNA::corAndPvalue(x = data.1, method = method)
    
    # Convert correlation matrix to long format
    df_cor <- cor_result$cor %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "from") %>%
      tidyr::pivot_longer(
        cols = -from, 
        names_to = "to", 
        values_to = "cor"
      )
    
    # Convert p-value matrix to long format
    df_pvalue <- cor_result$p %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "from") %>%
      tidyr::pivot_longer(
        cols = -from, 
        names_to = "to", 
        values_to = "pvalue"
      )
    
    # Combine and filter results
    result_cor <- dplyr::left_join(df_cor, df_pvalue, by = c("from", "to")) %>%
      dplyr::filter(
        pvalue <= !!pvalue,           # Significance threshold
        abs(cor) >= !!cor,            # Correlation strength threshold
        from != to                    # Exclude self-correlations
      ) %>%
      dplyr::arrange(desc(abs(cor)))  # Order by correlation strength
      
  } else {
    # Between-dataset correlation
    if (nrow(data.1) != nrow(data.2)) {
      stop("'data.1' and 'data.2' must have the same number of rows for between-dataset correlation")
    }
    
    cor_result <- WGCNA::corAndPvalue(x = data.1, y = data.2, method = method)
    
    # Convert correlation matrix to long format
    df_cor <- cor_result$cor %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "from") %>%
      tidyr::pivot_longer(
        cols = -from, 
        names_to = "to", 
        values_to = "cor"
      )
    
    # Convert p-value matrix to long format  
    df_pvalue <- cor_result$p %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "from") %>%
      tidyr::pivot_longer(
        cols = -from, 
        names_to = "to", 
        values_to = "pvalue"
      )
    
    # Combine and filter results
    result_cor <- dplyr::left_join(df_cor, df_pvalue, by = c("from", "to")) %>%
      dplyr::filter(
        pvalue <= !!pvalue,           # Significance threshold
        abs(cor) >= !!cor             # Correlation strength threshold
      ) %>%
      dplyr::arrange(desc(abs(cor)))  # Order by correlation strength
  }
  
  # Add correlation interpretation
  result_cor <- result_cor %>%
    dplyr::mutate(
      correlation_strength = dplyr::case_when(
        abs(cor) >= 0.9 ~ "Very strong",
        abs(cor) >= 0.7 ~ "Strong", 
        abs(cor) >= 0.5 ~ "Moderate",
        abs(cor) >= 0.3 ~ "Weak",
        TRUE ~ "Very weak"
      ),
      correlation_direction = ifelse(cor > 0, "Positive", "Negative")
    )
  
  # Add attributes to result for method and thresholds used
  attr(result_cor, "method") <- method
  attr(result_cor, "cor_threshold") <- cor
  attr(result_cor, "pvalue_threshold") <- pvalue
  attr(result_cor, "analysis_type") <- if (is.null(data.2)) "within-dataset" else "between-dataset"
  
  # Print summary if interactive
  if (interactive()) {
    cat("Correlation Analysis Summary:\n")
    cat("Method:", method, "\n")
    cat("Correlation threshold: |r| >=", cor, "\n")
    cat("P-value threshold: p <=", pvalue, "\n")
    cat("Analysis type:", attr(result_cor, "analysis_type"), "\n")
    cat("Significant correlations found:", nrow(result_cor), "\n\n")
  }
  
  return(result_cor)
}