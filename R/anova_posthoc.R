#' Perform ANOVA and Post-hoc Multiple Comparison Tests
#'
#' This function performs one-way ANOVA followed by post-hoc multiple comparison
#' tests using either Tukey HSD or Duncan's test. It provides both ANOVA p-values
#' and post-hoc group comparisons with significance letters.
#'
#' @param data A data frame containing the variables for analysis
#' @param group Column name (as string) containing group information for comparison
#' @param value Column name (as string) containing the numeric values to compare
#' @param method Method for post-hoc comparison. Either "Tukey" (default) for 
#'   Tukey's Honest Significant Difference test or "Duncan" for Duncan's Multiple 
#'   Range Test
#' @param level Confidence level for the comparisons (default: 0.95, i.e., 95%)
#'
#' @return A data frame containing:
#'   \itemize{
#'     \item \code{group}: Group levels
#'     \item \code{anova.pvalue}: Overall ANOVA p-value
#'     \item \code{anova.signif}: ANOVA significance levels (NS, *, **, ***)
#'     \item \code{Tukey.signif} or \code{Duncan.signif}: Post-hoc significance letters
#'   }
#'
#' @details
#' The function first performs a one-way ANOVA to test for overall differences
#' between groups. If significant differences are detected, post-hoc tests are
#' performed to identify which specific groups differ from each other.
#' 
#' Significance levels:
#' \itemize{
#'   \item NS: p > 0.05 (not significant)
#'   \item *: 0.01 < p ≤ 0.05
#'   \item **: 0.001 < p ≤ 0.01  
#'   \item ***: p ≤ 0.001
#' }
#'
#' @export
#'
#' @examples
#' library(bioRtools)
#' library(dplyr)
#'
#' # Load iris dataset
#' data(iris)
#'
#' # Perform ANOVA with Tukey post-hoc test
#' tukey_result <- anova_posthoc(
#'   data = iris,
#'   group = "Species", 
#'   value = "Sepal.Length"
#' )
#' print(tukey_result)
#'
#' # Perform ANOVA with Duncan post-hoc test
#' duncan_result <- anova_posthoc(
#'   data = iris,
#'   group = "Species",
#'   value = "Sepal.Length", 
#'   method = "Duncan"
#' )
#' print(duncan_result)
#'
#' # Using different confidence level
#' result_99 <- anova_posthoc(
#'   data = iris,
#'   group = "Species",
#'   value = "Sepal.Width",
#'   level = 0.99
#' )
#' print(result_99)
#'
anova_posthoc <- function(data, group, value, method = "Tukey", level = 0.95) {
  
  # Validate inputs
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  
  if (!group %in% names(data)) {
    stop(paste("Column '", group, "' not found in data"))
  }
  
  if (!value %in% names(data)) {
    stop(paste("Column '", value, "' not found in data"))
  }
  
  if (!method %in% c("Tukey", "Duncan")) {
    stop("'method' must be either 'Tukey' or 'Duncan'")
  }
  
  if (level <= 0 || level >= 1) {
    stop("'level' must be between 0 and 1")
  }
  
  # Prepare data for analysis
  {{ data }} %>%
    dplyr::select({{ group }}, {{ value }}) %>%
    magrittr::set_names(c("group.anova", "value")) %>%
    dplyr::mutate(group.anova = factor(group.anova, levels = unique(group.anova))) %>%
    # Remove missing values
    tidyr::drop_na() -> data.new
  
  # Check if we have enough data
  if (nrow(data.new) < 2) {
    stop("Not enough valid observations for analysis")
  }
  
  # Check if we have at least 2 groups
  if (length(unique(data.new$group.anova)) < 2) {
    stop("Need at least 2 groups for ANOVA")
  }
  
  # Perform ANOVA
  fit <- stats::aov(value ~ group.anova, data = data.new)
  anova_pvalue <- broom::tidy(fit)$p.value[1]
  
  # Perform post-hoc test based on method
  if (method == "Tukey") {
    multcomp::glht(fit, linfct = multcomp::mcp(group.anova = "Tukey")) %>%
      multcomp::cld(level = level, decreasing = TRUE) -> posthoc_result
    
    res <- posthoc_result[["mcletters"]][["Letters"]] %>%
      as.data.frame() %>%
      tibble::rownames_to_column() %>%
      magrittr::set_names(c("group", "Tukey.signif")) %>%
      dplyr::mutate(anova.pvalue = anova_pvalue) %>%
      dplyr::select(group, anova.pvalue, Tukey.signif)
      
  } else if (method == "Duncan") {
    duncan_result <- agricolae::duncan.test(
      fit, 
      "group.anova", 
      alpha = 1 - level, 
      console = FALSE  # Set to FALSE to avoid console output
    )
    
    res <- duncan_result[["groups"]] %>%
      as.data.frame() %>%
      tibble::rownames_to_column() %>%
      dplyr::mutate(anova.pvalue = anova_pvalue) %>%
      dplyr::select(rowname, anova.pvalue, groups) %>%
      magrittr::set_names(c("group", "anova.pvalue", "Duncan.signif"))
  }
  
  # Add ANOVA significance levels
  anova_result <- res %>%
    dplyr::mutate(
      anova.signif = dplyr::case_when(
        anova.pvalue > 0.05 ~ "NS",
        anova.pvalue > 0.01 ~ "*", 
        anova.pvalue > 0.001 ~ "**",
        TRUE ~ "***"
      )
    ) %>%
    # Reorder columns: group, anova.pvalue, anova.signif, posthoc results
    dplyr::select(group, anova.pvalue, anova.signif, dplyr::everything())
  
  return(anova_result)
}