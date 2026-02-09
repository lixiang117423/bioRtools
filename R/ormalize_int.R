#' Apply Rank-based Inverse Normal Transformation (INT)
#'
#' This function transforms each row of a data frame (representing features like
#' OTUs or genes) to a normal distribution using the Rank-based Inverse Normal
#' Transformation (INT). It's a non-parametric method robust to outliers and
#' suitable for preparing data for statistical models that assume normality,
#' such as GWAS.
#'
#' @description
#' The transformation follows the formula described by Blom:
#' `y = qnorm((r - c) / (N - 2*c + 1))`, where:
#' - `r` is the rank of each observation within a feature.
#' - `N` is the total number of non-missing observations for that feature.
#' - `c` is the Blom offset constant.
#' - `qnorm` is the standard normal quantile function.
#'
#' The function is designed to work within a `tidyverse` pipeline.
#'
#' @param data A data frame where rows are features (e.g., OTUs, genes) and
#'   columns are samples. It must contain one column with feature identifiers.
#' @param feature_col A string specifying the name of the column that contains
#'   the feature identifiers (e.g., `"otu_id"`).
#' @param c The Blom offset constant. Defaults to `3/8` as recommended.
#'
#' @return A data frame with the same dimensions as the input `data`, where the
#'   numeric values have been replaced by their INT-transformed values.
#'
#' @importFrom dplyr group_by mutate ungroup select
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom rlang sym
#' @importFrom stats qnorm
#' @importFrom tidyselect where
#'
#' @export
#'
#' @examples
#' library(bioRtools)
#' # 1. Create a sample dataset with non-normal data
#' set.seed(123)
#' sample_data <- data.frame(
#'   otu_id = paste0("OTU_", 1:5),
#'   sample_A = rpois(5, lambda = 2)^2,  # Skewed data
#'   sample_B = rpois(5, lambda = 1),
#'   sample_C = c(100, 0, 0, 0, 50),     # With outliers and zeros
#'   sample_D = rpois(5, lambda = 5)
#' )
#'
#' print("Original Data:")
#' print(sample_data)
#'
#' # 2. Apply the INT normalization
#' normalized_data <- normalize_int(sample_data, feature_col = "otu_id")
#'
#' print("Normalized Data:")
#' print(normalized_data)
#'
#' # 3. Verify the transformation for one feature (row)
#' # The result should be perfectly normal (or as close as ranks allow)
#' # and centered around 0.
#' transformed_row <- as.numeric(normalized_data[1, -1])
#' hist(transformed_row, main = "Histogram of a Transformed OTU", xlab = "INT Value")
#'
#' # Check for normality using a Q-Q plot
#' qqnorm(transformed_row)
#' qqline(transformed_row, col = "red")
#'
normalize_int <- function(data, feature_col, c = 3 / 8) {
  # --- Input Validation ---
  if (!is.data.frame(data)) {
    stop("Input 'data' must be a data frame.")
  }
  if (!is.character(feature_col) || length(feature_col) != 1) {
    stop("'feature_col' must be a single string.")
  }
  if (!feature_col %in% names(data)) {
    stop(paste0("Column '", feature_col, "' not found in the data frame."))
  }

  # Convert feature column name to a symbol for tidy evaluation
  feature_col_sym <- rlang::sym(feature_col)

  # --- Transformation Pipeline ---
  transformed_data <- data %>%
    # Reshape data from wide to long format for row-wise operations
    tidyr::pivot_longer(
      cols = tidyselect::where(is.numeric),
      names_to = "sample_id",
      values_to = "abundance"
    ) %>%
    # Group by each feature to perform calculations independently per row
    dplyr::group_by(!!feature_col_sym) %>%
    # Apply the INT formula within each group
    dplyr::mutate(
      # N is the count of non-NA observations for the current feature
      N = sum(!is.na(abundance)),
      # Rank observations, handling ties by averaging. NAs are kept as NAs.
      r = rank(abundance, na.last = "keep", ties.method = "average"),
      # Apply the full INT formula, overwriting the original abundance
      abundance = stats::qnorm((r - c) / (N - 2 * c + 1))
    ) %>%
    # Clean up intermediate columns before reshaping
    dplyr::select(!!feature_col_sym, sample_id, abundance) %>%
    # Reshape the data back to the original wide format
    tidyr::pivot_wider(
      names_from = "sample_id",
      values_from = "abundance"
    ) %>%
    # It's good practice to ungroup after all grouped operations are done
    dplyr::ungroup()

  return(transformed_data)
}
