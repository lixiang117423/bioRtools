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
#'   \item{result_permanova}{A data frame containing the PERMANOVA results table
#'     with degrees of freedom, sum of squares, mean squares, F values, R-squared,
#'     and p-values for each factor.}
#'   \item{raw_result}{The original \code{vegan::adonis2()} output object, which
#'     can be used for further analysis or custom formatting.}
#'   \item{summary_stats}{A simplified summary table with key statistics including
#'     R-squared and p-values for each factor.}
#' }
#'
#' @details
#' PERMANOVA tests the null hypothesis that the centroids and dispersion of groups
#' as defined by the explanatory variables are equivalent for all groups. It is a
#' non-parametric method suitable for ecological community data. The function uses
#' \code{vegan::adonis2()} which is the updated version of \code{vegan::adonis()}.
#' Response row names and metadata sample IDs must identify the same complete,
#' unique sample set; metadata are reordered automatically before fitting.
#'
#' @note
#' This function uses \code{vegan::adonis2()} instead of the deprecated
#' \code{vegan::adonis()} for better performance and accuracy. If you need to use
#' the old \code{adonis()} function for compatibility reasons, please specify this
#' in your analysis requirements.
#'
#' @author Xiang LI <lixiang117423@gmail.com>
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
#' ) -> result_permanova
#'
#' # Access different components
#' result_permanova$result_permanova  # Main results table
#' result_permanova$summary_stats     # Summary statistics
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

  # The response is supplied through the formula. Align metadata by row names
  # before fitting so a reordered metadata table cannot silently change results.
  response <- tryCatch(
    eval(formula[[2]], envir = environment(formula)),
    error = function(e) {
      stop("Could not evaluate the response in 'formula': ", e$message, call. = FALSE)
    }
  )
  response_ids <- rownames(response)
  if (is.null(response_ids) || has_default_rownames(response)) {
    stop(
      "The response in 'formula' must have non-default sample IDs in row names",
      call. = FALSE
    )
  }
  sample <- align_sample_metadata(sample, response_ids)$sample

  # Rebuild the formula around the evaluated response. This keeps adonis2 from
  # re-evaluating a caller-local expression such as `dataset$matrix` inside the
  # metadata data frame.
  fit_formula <- formula
  fit_formula[[2]] <- as.name(".bioRtools_response")
  fit_environment <- new.env(parent = environment(formula))
  fit_environment$.bioRtools_response <- response
  environment(fit_formula) <- fit_environment
  .bioRtools_response <- response

  # Perform PERMANOVA using the updated adonis2 function
  tryCatch(
    {
      permanova_result <- vegan::adonis2(
        formula = fit_formula,
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
  result_permanova <- permanova_result %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "factor") %>%
    dplyr::mutate(
      factor = ifelse(factor == "Residual", "Residuals", factor),
      factor = ifelse(factor == "Total", "Total", factor)
    )

  # Create summary statistics table for significant factors only
  summary_stats <- result_permanova %>%
    dplyr::filter(
      !factor %in% c("Residuals", "Total"),
      !is.na(`Pr(>F)`)
    ) %>%
    dplyr::select(
      factor,
      R2,
      `Pr(>F)`
    ) %>%
    dplyr::rename(
      p_value = `Pr(>F)`,
      r_squared = R2
    ) %>%
    dplyr::mutate(
      significance = dplyr::case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        p_value < 0.1 ~ ".",
        TRUE ~ ""
      )
    ) %>%
    dplyr::arrange(p_value)

  # Return comprehensive results using the shared result contract
  results <- list(
    result_permanova = result_permanova,
    raw_result = permanova_result,
    summary_stats = summary_stats
  )

  new_bioRtools_fit(
    results,
    model = permanova_result,
    tables = list(
      permanova = result_permanova,
      summary = summary_stats
    ),
    plots = list(),
    analysis_type = "PERMANOVA",
    engine = "vegan::adonis2",
    engine_versions = c(vegan = as.character(utils::packageVersion("vegan"))),
    data_info = list(
      n_samples = nrow(sample),
      distance = method,
      permutations = permutations
    ),
    call = match.call()
  )
}
