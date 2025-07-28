#' Identify Statistical Outliers Using Interquartile Range (IQR) Method
#'
#' This function identifies statistical outliers in numeric data using the
#' standard interquartile range (IQR) method, also known as Tukey's fences.
#' Data points that fall more than 1.5 × IQR beyond the first or third quartile
#' are classified as outliers. This is the same method used in box plots to
#' identify and visualize outliers.
#'
#' @param x Numeric vector containing the data to be analyzed for outliers.
#'   Missing values (NA) are handled gracefully and will not be classified
#'   as outliers
#' @param method Character string specifying the outlier detection method
#'   (default: "iqr"). Currently supported methods:
#'   \itemize{
#'     \item \code{"iqr"}: Interquartile range method (Tukey's fences)
#'     \item \code{"modified_zscore"}: Modified Z-score using median absolute deviation
#'     \item \code{"zscore"}: Standard Z-score method (mean ± k*SD)
#'   }
#' @param k Multiplier for the outlier threshold (default: 1.5 for IQR method).
#'   Common values:
#'   \itemize{
#'     \item 1.5: Standard outlier threshold (used in box plots)
#'     \item 3.0: More conservative threshold for extreme outliers
#'     \item 2.5: Intermediate threshold
#'   }
#' @param return.logical Logical indicating return format (default: FALSE).
#'   If TRUE, returns logical vector (TRUE/FALSE). If FALSE, returns character
#'   vector ("yes"/"no") for backward compatibility
#' @param na.rm Logical indicating whether to remove missing values before
#'   calculation (default: TRUE). If FALSE and NAs are present, the entire
#'   result will be NA
#'
#' @return Vector indicating outlier status for each element in x:
#'   \itemize{
#'     \item If \code{return.logical = TRUE}: Logical vector where TRUE indicates outlier
#'     \item If \code{return.logical = FALSE}: Character vector with "yes" for outliers, "no" for normal values
#'     \item Missing values in input are preserved as NA in output
#'   }
#'
#' @details
#' \strong{IQR Method (Tukey's Fences):}
#' \enumerate{
#'   \item Calculate first quartile (Q1) and third quartile (Q3)
#'   \item Compute IQR = Q3 - Q1
#'   \item Define lower fence = Q1 - k × IQR
#'   \item Define upper fence = Q3 + k × IQR
#'   \item Values outside these fences are outliers
#' }
#' 
#' \strong{Modified Z-Score Method:}
#' Uses median and median absolute deviation (MAD) instead of mean and standard
#' deviation, making it more robust to outliers:
#' \deqn{Modified Z-score = 0.6745 × (x - median) / MAD}
#' 
#' \strong{Standard Z-Score Method:}
#' \deqn{Z-score = (x - mean) / SD}
#' 
#' The IQR method is recommended for most applications as it:
#' \itemize{
#'   \item Is robust to extreme values
#'   \item Works well with skewed distributions
#'   \item Provides interpretable thresholds
#'   \item Is widely used and accepted
#' }
#'
#' @note
#' \itemize{
#'   \item The IQR method assumes the data has sufficient spread (IQR > 0)
#'   \item For small sample sizes (n < 10), outlier detection may be unreliable
#'   \item Consider the biological/scientific context when interpreting outliers
#'   \item Outliers are not necessarily erroneous; they may represent important biological variation
#' }
#'
#' @references
#' Tukey, J.W. (1977). Exploratory Data Analysis. Addison-Wesley.
#' 
#' Iglewicz, B. and Hoaglin, D.C. (1993). How to Detect and Handle Outliers. 
#' ASQC Quality Press.
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[stats]{quantile}} for quantile calculation
#'   \item \code{\link[stats]{IQR}} for interquartile range
#'   \item \code{\link[graphics]{boxplot}} for visual outlier identification
#' }
#'
#' @export
#'
#' @examples
#' library(bioRtools)
#' library(dplyr)
#' 
#' # Example 1: Basic outlier detection
#' data(iris)
#' 
#' # Identify outliers in sepal length
#' sepal_outliers <- find_outliers(iris$Sepal.Length)
#' print(table(sepal_outliers))
#' 
#' # Show outlier values
#' outlier_values <- iris$Sepal.Length[sepal_outliers == "yes"]
#' print("Outlier values:")
#' print(sort(outlier_values))
#' 
#' # Example 2: Group-wise outlier detection (original use case)
#' iris_with_outliers <- iris %>%
#'   group_by(Species) %>%
#'   mutate(
#'     sepal_length_outlier = find_outliers(Sepal.Length),
#'     sepal_width_outlier = find_outliers(Sepal.Width),
#'     petal_length_outlier = find_outliers(Petal.Length),
#'     petal_width_outlier = find_outliers(Petal.Width)
#'   ) %>%
#'   ungroup()
#' 
#' # Summary of outliers by species
#' outlier_summary <- iris_with_outliers %>%
#'   group_by(Species) %>%
#'   summarise(
#'     n_samples = n(),
#'     sepal_length_outliers = sum(sepal_length_outlier == "yes"),
#'     sepal_width_outliers = sum(sepal_width_outlier == "yes"),
#'     petal_length_outliers = sum(petal_length_outlier == "yes"),
#'     petal_width_outliers = sum(petal_width_outlier == "yes"),
#'     .groups = 'drop'
#'   )
#' 
#' print("Outliers by species and measurement:")
#' print(outlier_summary)
#' 
#' # Example 3: Different outlier detection methods
#' set.seed(123)
#' test_data <- c(rnorm(95, mean = 10, sd = 2), c(1, 25, 30))  # Normal data + outliers
#' 
#' # Compare different methods
#' iqr_outliers <- find_outliers(test_data, method = "iqr", return.logical = TRUE)
#' zscore_outliers <- find_outliers(test_data, method = "zscore", k = 2, return.logical = TRUE)
#' modified_zscore_outliers <- find_outliers(test_data, method = "modified_zscore", 
#'                                           k = 3.5, return.logical = TRUE)
#' 
#' # Create comparison data frame
#' comparison <- data.frame(
#'   value = test_data,
#'   iqr_method = iqr_outliers,
#'   zscore_method = zscore_outliers,
#'   modified_zscore_method = modified_zscore_outliers
#' )
#' 
#' # Show detected outliers
#' outlier_comparison <- comparison[apply(comparison[,2:4], 1, any), ]
#' print("Outliers detected by different methods:")
#' print(outlier_comparison)
#' 
#' # Method agreement
#' print("Method comparison:")
#' print(paste("IQR method:", sum(iqr_outliers), "outliers"))
#' print(paste("Z-score method:", sum(zscore_outliers), "outliers"))
#' print(paste("Modified Z-score method:", sum(modified_zscore_outliers), "outliers"))
#' 
#' # Example 4: Custom threshold values
#' # Conservative outlier detection (fewer outliers)
#' conservative_outliers <- find_outliers(iris$Sepal.Length, k = 3.0)
#' 
#' # Liberal outlier detection (more outliers)  
#' liberal_outliers <- find_outliers(iris$Sepal.Length, k = 1.0)
#' 
#' print("Threshold comparison:")
#' print(paste("Standard (k=1.5):", sum(find_outliers(iris$Sepal.Length) == "yes"), "outliers"))
#' print(paste("Conservative (k=3.0):", sum(conservative_outliers == "yes"), "outliers"))
#' print(paste("Liberal (k=1.0):", sum(liberal_outliers == "yes"), "outliers"))
#' 
#' # Example 5: Handling missing values
#' data_with_na <- c(1, 2, 3, NA, 4, 5, 100, NA, 6, 7)
#' 
#' # Default behavior (na.rm = TRUE)
#' outliers_default <- find_outliers(data_with_na)
#' print("With missing values (na.rm = TRUE):")
#' print(data.frame(value = data_with_na, outlier = outliers_default))
#' 
#' # Keep NAs in calculation (na.rm = FALSE)
#' outliers_with_na <- find_outliers(data_with_na, na.rm = FALSE)
#' print("With missing values (na.rm = FALSE):")
#' print(data.frame(value = data_with_na, outlier = outliers_with_na))
#' 
#' # Example 6: Biological interpretation
#' # Gene expression data simulation
#' set.seed(456)
#' gene_expression <- c(
#'   rnorm(50, mean = 5, sd = 0.5),    # Normal expression
#'   c(2, 8.5, 9.2)                    # Potential outliers
#' )
#' 
#' expression_analysis <- data.frame(
#'   sample_id = paste0("Sample_", 1:length(gene_expression)),
#'   expression = gene_expression,
#'   is_outlier = find_outliers(gene_expression),
#'   outlier_logical = find_outliers(gene_expression, return.logical = TRUE)
#' )
#' 
#' # Identify samples with outlier expression
#' outlier_samples <- expression_analysis[expression_analysis$is_outlier == "yes", ]
#' print("Samples with outlier gene expression:")
#' print(outlier_samples)
#' 
#' # Statistical summary
#' print("Expression data summary:")
#' print(summary(gene_expression))
#' print(paste("IQR:", round(IQR(gene_expression), 2)))
#' print(paste("Outlier threshold (lower):", 
#'             round(quantile(gene_expression, 0.25) - 1.5 * IQR(gene_expression), 2)))
#' print(paste("Outlier threshold (upper):", 
#'             round(quantile(gene_expression, 0.75) + 1.5 * IQR(gene_expression), 2)))
#'
find_outliers <- function(x, method = "iqr", k = 1.5, return.logical = FALSE, na.rm = TRUE) {
  
  # Input validation
  if (!is.numeric(x)) {
    stop("'x' must be a numeric vector")
  }
  
  if (length(x) == 0) {
    stop("'x' cannot be empty")
  }
  
  if (!is.character(method) || length(method) != 1) {
    stop("'method' must be a single character string")
  }
  
  valid_methods <- c("iqr", "zscore", "modified_zscore")
  if (!method %in% valid_methods) {
    stop(paste("'method' must be one of:", paste(valid_methods, collapse = ", ")))
  }
  
  if (!is.numeric(k) || length(k) != 1 || k <= 0) {
    stop("'k' must be a single positive number")
  }
  
  if (!is.logical(return.logical) || length(return.logical) != 1) {
    stop("'return.logical' must be a single logical value")
  }
  
  if (!is.logical(na.rm) || length(na.rm) != 1) {
    stop("'na.rm' must be a single logical value")
  }
  
  # Handle all-NA case
  if (all(is.na(x))) {
    result <- rep(NA, length(x))
    if (!return.logical) {
      result <- as.character(result)
    }
    return(result)
  }
  
  # Check for sufficient data
  valid_data <- x[!is.na(x)]
  if (length(valid_data) < 4) {
    warning("Very few data points available. Outlier detection may be unreliable.")
  }
  
  # Detect outliers based on method
  if (method == "iqr") {
    # IQR method (Tukey's fences)
    if (na.rm) {
      q1 <- stats::quantile(x, 0.25, na.rm = TRUE)
      q3 <- stats::quantile(x, 0.75, na.rm = TRUE)
      iqr_value <- stats::IQR(x, na.rm = TRUE)
    } else {
      if (any(is.na(x))) {
        # If na.rm = FALSE and NAs present, return all NAs
        result <- rep(NA, length(x))
        if (!return.logical) {
          result <- as.character(result)
        }
        return(result)
      }
      q1 <- stats::quantile(x, 0.25)
      q3 <- stats::quantile(x, 0.75)
      iqr_value <- stats::IQR(x)
    }
    
    # Check for zero IQR (no variation in data)
    if (iqr_value == 0) {
      warning("IQR is zero (no variation in data). No outliers will be detected.")
      is_outlier <- rep(FALSE, length(x))
    } else {
      # Calculate outlier boundaries
      lower_fence <- q1 - k * iqr_value
      upper_fence <- q3 + k * iqr_value
      
      # Identify outliers
      is_outlier <- (x < lower_fence | x > upper_fence) & !is.na(x)
    }
    
  } else if (method == "zscore") {
    # Standard Z-score method
    if (na.rm) {
      mean_x <- mean(x, na.rm = TRUE)
      sd_x <- sd(x, na.rm = TRUE)
    } else {
      if (any(is.na(x))) {
        result <- rep(NA, length(x))
        if (!return.logical) {
          result <- as.character(result)
        }
        return(result)
      }
      mean_x <- mean(x)
      sd_x <- sd(x)
    }
    
    if (sd_x == 0) {
      warning("Standard deviation is zero (no variation in data). No outliers will be detected.")
      is_outlier <- rep(FALSE, length(x))
    } else {
      z_scores <- abs((x - mean_x) / sd_x)
      is_outlier <- (z_scores > k) & !is.na(z_scores)
    }
    
  } else if (method == "modified_zscore") {
    # Modified Z-score using median and MAD
    if (na.rm) {
      median_x <- median(x, na.rm = TRUE)
      mad_x <- mad(x, na.rm = TRUE)
    } else {
      if (any(is.na(x))) {
        result <- rep(NA, length(x))
        if (!return.logical) {
          result <- as.character(result)
        }
        return(result)
      }
      median_x <- median(x)
      mad_x <- mad(x)
    }
    
    if (mad_x == 0) {
      warning("MAD is zero (no variation in data). No outliers will be detected.")
      is_outlier <- rep(FALSE, length(x))
    } else {
      # Modified Z-score formula
      modified_z_scores <- abs(0.6745 * (x - median_x) / mad_x)
      is_outlier <- (modified_z_scores > k) & !is.na(modified_z_scores)
    }
  }
  
  # Handle NAs in the original data
  is_outlier[is.na(x)] <- NA
  
  # Return in requested format
  if (return.logical) {
    return(is_outlier)
  } else {
    # Convert to character for backward compatibility
    result <- ifelse(is_outlier, "yes", "no")
    result[is.na(is_outlier)] <- NA
    return(result)
  }
}