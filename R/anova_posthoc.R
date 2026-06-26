#' Perform ANOVA and Post-hoc Multiple Comparison Tests
#'
#' This function performs one-way ANOVA followed by post-hoc multiple comparison
#' tests using either Tukey HSD or Duncan's test. It provides both ANOVA p-values
#' and post-hoc group comparisons with significance letters. Supports grouped data
#' from dplyr::group_by().
#'
#' @param data A data frame or grouped data frame containing the variables for analysis
#' @param group Either a column name (as string) containing group information,
#'   or a formula like \code{value ~ group} (e.g., \code{Sepal.Length ~ Species}).
#'   When a formula is passed, the \code{value} parameter is ignored.
#' @param value Column name (as string) containing the numeric values to compare.
#'   Ignored when \code{group} is a formula.
#' @param method Method for post-hoc comparison. Either "Tukey" (default) for
#'   Tukey's Honest Significant Difference test or "Duncan" for Duncan's Multiple
#'   Range Test
#' @param level Confidence level for the comparisons (default: 0.95, i.e., 95%)
#'
#' @return A data frame containing:
#'   - Grouping columns (if data was grouped)
#'   - \code{group}: Group levels
#'   - \code{anova_pvalue}: Overall ANOVA p-value
#'   - \code{anova_signif}: ANOVA significance levels (NS, *, **, ***)
#'   - \code{tukey_signif} or \code{duncan_signif}: Post-hoc significance letters
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
#'   - NS: p > 0.05 (not significant)
#'   - *: 0.01 < p ≤ 0.05
#'   - **: 0.001 < p ≤ 0.01
#'   - ***: p ≤ 0.001
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' library(bioRtools)
#' library(dplyr)
#'
#' # Load iris dataset
#' data(iris)
#'
#' # Basic usage: formula interface
#' tukey_result <- anova_posthoc(iris, Sepal.Length ~ Species)
#' print(tukey_result)
#'
#' # Traditional string interface (still supported)
#' tukey_result <- anova_posthoc(
#'   data = iris,
#'   group = "Species",
#'   value = "Sepal.Length"
#' )
#'
#' # With Duncan post-hoc test
#' duncan_result <- anova_posthoc(iris, Sepal.Length ~ Species, method = "Duncan")
#'
#' # Using grouped data
#' iris %>%
#'   dplyr::mutate(size = ifelse(Sepal.Width > median(Sepal.Width), "Large", "Small")) %>%
#'   dplyr::group_by(size) %>%
#'   anova_posthoc(Sepal.Length ~ Species)
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
anova_posthoc <- function(data, group, value = NULL, method = "Tukey", level = 0.95) {

  # Formula interface: group can be a formula like value ~ treatment
  if (inherits(group, "formula")) {
    if (length(group) != 3) {
      stop("Formula must be in the form 'value ~ group'")
    }
    value <- as.character(group[[2]])
    group <- as.character(group[[3]])
  }

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

  # Rename the output column "group" to the original input name on the way out.
  finalize <- function(df) {
    if ("group" %in% names(df)) {
      names(df)[names(df) == "group"] <- group
    }
    df
  }

  # Prepare data for analysis
  data_new <- data %>%
    dplyr::select(dplyr::all_of(c(group, value))) %>%
    magrittr::set_names(c("group_anova", "value")) %>%
    # Remove missing values
    tidyr::drop_na()

  # Order factor levels by ascending mean so cld(decreasing=TRUE) assigns
  # 'a' to the group with the highest mean
  mean_order <- data_new %>%
    dplyr::group_by(group_anova) %>%
    dplyr::summarise(mean_val = mean(value, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(mean_val) %>%
    dplyr::pull(group_anova)

  data_new$group_anova <- factor(data_new$group_anova, levels = mean_order)

  # Check if we have enough data
  if (nrow(data_new) < 2) {
    warning("Not enough valid observations for analysis in this group")
    return(finalize(tibble::tibble(
      group = character(),
      anova_pvalue = numeric(),
      anova_signif = character()
    )))
  }

  # Check if we have at least 2 groups
  n_groups <- length(unique(data_new$group_anova))
  if (n_groups < 2) {
    warning("Need at least 2 groups for ANOVA in this group")
    return(finalize(tibble::tibble(
      group = as.character(unique(data_new$group_anova)),
      anova_pvalue = NA_real_,
      anova_signif = "NS"
    )))
  }

  # Perform ANOVA
  fit <- tryCatch(
    {
      stats::aov(value ~ group_anova, data = data_new)
    },
    error = function(e) {
      warning(paste("ANOVA failed:", e$message))
      NULL
    })

  if (is.null(fit)) {
    return(finalize(tibble::tibble(
      group = character(),
      anova_pvalue = numeric(),
      anova_signif = character()
    )))
  }

  anova_pvalue <- broom::tidy(fit)$p.value[1]

  # Perform post-hoc test based on method
  if (method == "Tukey") {
    posthoc_result <- tryCatch(
      {
        multcomp::glht(fit, linfct = multcomp::mcp(group_anova = "Tukey")) %>%
          multcomp::cld(level = level, decreasing = TRUE)
      },
      error = function(e) {
        warning(paste("Tukey post-hoc test failed:", e$message))
        NULL
      })

    if (is.null(posthoc_result)) {
      return(finalize(tibble::tibble(
        group = as.character(unique(data_new$group_anova)),
        anova_pvalue = anova_pvalue,
        anova_signif = get_anova_signif(anova_pvalue),
        tukey_signif = NA_character_
      )))
    }

    res <- posthoc_result[["mcletters"]][["Letters"]] %>%
      as.data.frame() %>%
      tibble::rownames_to_column() %>%
      magrittr::set_names(c("group", "tukey_signif")) %>%
      dplyr::mutate(anova_pvalue = anova_pvalue) %>%
      dplyr::select(group, anova_pvalue, tukey_signif)

  } else if (method == "Duncan") {
    duncan_result <- tryCatch(
      {
        agricolae::duncan.test(
          fit,
          "group_anova",
          alpha = 1 - level,
          console = FALSE
        )
      },
      error = function(e) {
        warning(paste("Duncan post-hoc test failed:", e$message))
        NULL
      })

    if (is.null(duncan_result)) {
      return(finalize(tibble::tibble(
        group = as.character(unique(data_new$group_anova)),
        anova_pvalue = anova_pvalue,
        anova_signif = get_anova_signif(anova_pvalue),
        duncan_signif = NA_character_
      )))
    }

    res <- duncan_result[["groups"]] %>%
      as.data.frame() %>%
      tibble::rownames_to_column() %>%
      dplyr::mutate(anova_pvalue = anova_pvalue) %>%
      dplyr::select(rowname, anova_pvalue, groups) %>%
      magrittr::set_names(c("group", "anova_pvalue", "duncan_signif"))
  }

  # Add ANOVA significance levels
  anova_result <- res %>%
    dplyr::mutate(
      anova_signif = get_anova_signif(anova_pvalue)
    ) %>%
    # Reorder columns: group, anova_pvalue, anova_signif, posthoc results
    dplyr::select(group, anova_pvalue, anova_signif, dplyr::everything())

  finalize(anova_result)
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
