#' Perform ANOVA and Post-hoc Multiple Comparison Tests
#'
#' This function performs one-way ANOVA followed by post-hoc multiple comparison
#' tests using either Tukey HSD or Duncan's test. It provides both ANOVA p-values
#' and post-hoc group comparisons with significance letters. Supports grouped data
#' from dplyr::group_by().
#'
#' @param data A data frame or grouped data frame containing the variables for analysis
#' @param group Column name (as string) containing group information for comparison
#' @param value Column name (as string) containing the numeric values to compare
#' @param method Method for post-hoc comparison. Either "Tukey" (default) for
#'   Tukey's Honest Significant Difference test or "Duncan" for Duncan's Multiple
#'   Range Test
#' @param level Confidence level for the comparisons (default: 0.95, i.e., 95%)
#'
#' @return A data frame containing:
#'   \itemize{
#'     \item Grouping columns (if data was grouped)
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
#' When used with grouped data (from dplyr::group_by()), the function will
#' perform separate ANOVA analyses for each group combination.
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
#' # Basic usage: Perform ANOVA with Tukey post-hoc test
#' tukey_result <- anova_posthoc(
#'   data = iris,
#'   group = "Species",
#'   value = "Sepal.Length"
#' )
#' print(tukey_result)
#'
#' # With Duncan post-hoc test
#' duncan_result <- anova_posthoc(
#'   data = iris,
#'   group = "Species",
#'   value = "Sepal.Length",
#'   method = "Duncan"
#' )
#' print(duncan_result)
#'
#' # Using grouped data
#' library(tidyr)
#'
#' # Create example data with multiple factors
#' iris_long <- iris %>%
#'   mutate(Size = ifelse(Sepal.Width > median(Sepal.Width), "Large", "Small")) %>%
#'   group_by(Size)
#'
#' # Perform ANOVA for each Size group
#' grouped_result <- anova_posthoc(
#'   data = iris_long,
#'   group = "Species",
#'   value = "Sepal.Length"
#' )
#' print(grouped_result)
#'
#' # Multiple grouping variables
#' mtcars_grouped <- mtcars %>%
#'   mutate(cyl = as.factor(cyl),
#'     am = as.factor(am),
#'     gear_group = as.factor(gear)) %>%
#'   group_by(cyl, am)
#'
#' result <- anova_posthoc(
#'   data = mtcars_grouped,
#'   group = "gear_group",
#'   value = "mpg",
#'   method = "Tukey"
#' )
#' print(result)
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

  # Check if data is grouped
  is_grouped <- dplyr::is_grouped_df(data)

  if (is_grouped) {
    # Get grouping variables
    group_vars <- dplyr::group_vars(data)

    # Process each group separately
    result <- data %>%
      dplyr::group_modify(~ {
        anova_single_group(.x, group, value, method, level)
      }) %>%
      dplyr::ungroup()

    return(result)
  } else {
    # Process ungrouped data
    return(anova_single_group(data, group, value, method, level))
  }
}

#' Internal function to perform ANOVA on a single group
#' @noRd
anova_single_group <- function(data, group, value, method, level) {
  # Prepare data for analysis
  data.new <- data %>%
    dplyr::select(dplyr::all_of(c(group, value))) %>%
    magrittr::set_names(c("group.anova", "value")) %>%
    dplyr::mutate(group.anova = factor(group.anova, levels = unique(group.anova))) %>%
    # Remove missing values
    tidyr::drop_na()

  # Check if we have enough data
  if (nrow(data.new) < 2) {
    warning("Not enough valid observations for analysis in this group")
    return(tibble::tibble(
      group = character(),
      anova.pvalue = numeric(),
      anova.signif = character()
    ))
  }

  # Check if we have at least 2 groups
  n_groups <- length(unique(data.new$group.anova))
  if (n_groups < 2) {
    warning("Need at least 2 groups for ANOVA in this group")
    return(tibble::tibble(
      group = as.character(unique(data.new$group.anova)),
      anova.pvalue = NA_real_,
      anova.signif = "NS"
    ))
  }

  # Perform ANOVA
  fit <- tryCatch(
    {
      stats::aov(value ~ group.anova, data = data.new)
    },
    error = function(e) {
      warning(paste("ANOVA failed:", e$message))
      return(NULL)
    })

  if (is.null(fit)) {
    return(tibble::tibble(
      group = character(),
      anova.pvalue = numeric(),
      anova.signif = character()
    ))
  }

  anova_pvalue <- broom::tidy(fit)$p.value[1]

  # Perform post-hoc test based on method
  if (method == "Tukey") {
    posthoc_result <- tryCatch(
      {
        multcomp::glht(fit, linfct = multcomp::mcp(group.anova = "Tukey")) %>%
          multcomp::cld(level = level, decreasing = TRUE)
      },
      error = function(e) {
        warning(paste("Tukey post-hoc test failed:", e$message))
        return(NULL)
      })

    if (is.null(posthoc_result)) {
      return(tibble::tibble(
        group = as.character(unique(data.new$group.anova)),
        anova.pvalue = anova_pvalue,
        anova.signif = get_anova_signif(anova_pvalue),
        Tukey.signif = NA_character_
      ))
    }

    res <- posthoc_result[["mcletters"]][["Letters"]] %>%
      as.data.frame() %>%
      tibble::rownames_to_column() %>%
      magrittr::set_names(c("group", "Tukey.signif")) %>%
      dplyr::mutate(anova.pvalue = anova_pvalue) %>%
      dplyr::select(group, anova.pvalue, Tukey.signif)

  } else if (method == "Duncan") {
    duncan_result <- tryCatch(
      {
        agricolae::duncan.test(
          fit,
          "group.anova",
          alpha = 1 - level,
          console = FALSE
        )
      },
      error = function(e) {
        warning(paste("Duncan post-hoc test failed:", e$message))
        return(NULL)
      })

    if (is.null(duncan_result)) {
      return(tibble::tibble(
        group = as.character(unique(data.new$group.anova)),
        anova.pvalue = anova_pvalue,
        anova.signif = get_anova_signif(anova_pvalue),
        Duncan.signif = NA_character_
      ))
    }

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
      anova.signif = get_anova_signif(anova.pvalue)
    ) %>%
    # Reorder columns: group, anova.pvalue, anova.signif, posthoc results
    dplyr::select(group, anova.pvalue, anova.signif, dplyr::everything())

  return(anova_result)
}

#' Helper function to get ANOVA significance label
#' @noRd
get_anova_signif <- function(pvalue) {
  dplyr::case_when(
    pvalue > 0.05 ~ "NS",
    pvalue > 0.01 ~ "*",
    pvalue > 0.001 ~ "**",
    TRUE ~ "***"
  )
}
