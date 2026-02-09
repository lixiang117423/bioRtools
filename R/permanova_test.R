#' Perform Permutational Multivariate Analysis of Variance (PERMANOVA)
#'
#' @param formula A formula specifying the response data and explanatory variables,
#'   such as \code{data ~ factor1 * factor2 * factor3}. The left-hand side should be
#'   a data frame or matrix containing community data (e.g., OTU table), and the
#'   right-hand side should specify the explanatory variables.
#' @param sample Sample information table containing the explanatory variables
#'   referenced in the formula. Column names must match variable names in the formula.
#' @param method The distance method used in \code{vegan::vegdist()} to calculate
#'   pairwise distances. Default is "bray" (Bray-Curtis dissimilarity). Other options
#'   include "jaccard", "euclidean", "manhattan", etc.
#' @param permutations The number of permutations for the permutation test, or a
#'   list of control values for permutations as returned by \code{vegan::how()},
#'   or a permutation matrix where each row gives the permuted indices. Default is 999.
#'
#' @return A list containing three components:
#' \describe{
#'   \item{result.permanova}{A data frame containing the PERMANOVA results table
#'     with degrees of freedom, sum of squares, mean squares, F values, R-squared,
#'     and p-values for each factor.}
#'   \item{raw.result}{The original \code{vegan::adonis2()} output object, which
#'     can be used for further analysis or custom formatting.}
#'   \item{summary.stats}{A simplified summary table with key statistics including
#'     R-squared and p-values for each factor.}
#' }
#'
#' @details
#' PERMANOVA tests the null hypothesis that the centroids and dispersion of groups
#' as defined by the explanatory variables are equivalent for all groups. It is a
#' non-parametric method suitable for ecological community data. The function uses
#' \code{vegan::adonis2()} which is the updated version of \code{vegan::adonis()}.
#'
#' @note
#' This function uses \code{vegan::adonis2()} instead of the deprecated
#' \code{vegan::adonis()} for better performance and accuracy. If you need to use
#' the old \code{adonis()} function for compatibility reasons, please specify this
#' in your analysis requirements.
#'
#' @export
#'
#' @examples
#' library(dplyr)
#' library(bioRtools)
#'
#' data(df.permanova.otu)
#' data(df.permanova.sample)
#'
#' # Basic PERMANOVA test
#' permanova_test(
#'   df.permanova.otu ~ location * variety * group,
#'   df.permanova.sample
#' ) -> result.permanova
#'
#' # Access different components
#' result.permanova$result.permanova  # Main results table
#' result.permanova$summary.stats     # Summary statistics
#'
#' # Using different distance method and permutations
#' permanova_test(
#'   df.permanova.otu ~ location * variety,
#'   df.permanova.sample,
#'   method = "jaccard",
#'   permutations = 9999
#' ) -> result.jaccard
#'
permanova_test <- function(formula, sample, method = "bray", permutations = 999) {
  # Input validation
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a valid formula object")
  }

  if (!is.data.frame(sample)) {
    stop("'sample' must be a data frame")
  }

  if (!is.character(method) || length(method) != 1) {
    stop("'method' must be a single character string")
  }

  if (!is.numeric(permutations) || length(permutations) != 1 || permutations < 1) {
    if (!inherits(permutations, c("how", "matrix"))) {
      stop("'permutations' must be a positive integer, a 'how' object, or a permutation matrix")
    }
  }

  # Perform PERMANOVA using the updated adonis2 function
  tryCatch(
    {
      permanova.result <- vegan::adonis2(
        formula = formula,
        data = sample,
        permutations = permutations,
        method = method,
        by = "terms"  # Analyze terms sequentially
      )
    },
    error = function(e) {
      stop("PERMANOVA analysis failed: ", e$message, call. = FALSE)
    })

  # Extract main results table
  result.permanova <- permanova.result %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Factor") %>%
    dplyr::mutate(
      Factor = ifelse(Factor == "Residual", "Residuals", Factor),
      Factor = ifelse(Factor == "Total", "Total", Factor)
    )

  # Create summary statistics table for significant factors only
  summary.stats <- result.permanova %>%
    dplyr::filter(
      !Factor %in% c("Residuals", "Total"),
      !is.na(`Pr(>F)`)
    ) %>%
    dplyr::select(
      Factor,
      R2,
      `Pr(>F)`
    ) %>%
    dplyr::rename(
      P_value = `Pr(>F)`,
      R_squared = R2
    ) %>%
    dplyr::mutate(
      Significance = dplyr::case_when(
        P_value < 0.001 ~ "***",
        P_value < 0.01 ~ "**",
        P_value < 0.05 ~ "*",
        P_value < 0.1 ~ ".",
        TRUE ~ ""
      )
    ) %>%
    dplyr::arrange(P_value)

  # Return comprehensive results following project conventions
  return(list(
    result.permanova = result.permanova,
    raw.result = permanova.result,
    summary.stats = summary.stats
  ))
}
