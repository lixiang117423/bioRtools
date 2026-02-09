#' Extract R-squared and P-value from Linear Model
#'
#' @description
#' Extract adjusted R-squared and overall model p-value from a fitted linear
#' model object. This is useful for adding statistics to plots or creating
#' summary tables of multiple models.
#'
#' @param model A fitted linear model object from \code{lm()} or compatible
#'   modeling functions that have \code{summary()} and \code{anova()} methods.
#' @param r_squared_type Character string specifying which R-squared to extract.
#'   Options are "adjusted" (default) or "multiple". The adjusted R-squared
#'   is generally preferred as it accounts for the number of predictors.
#' @param digits Integer. Number of decimal places to round R-squared value.
#'   Default is 3. Use NULL for no rounding.
#'
#' @details
#' The function extracts:
#' \itemize{
#'   \item \strong{R-squared}: By default, adjusted R-squared which adjusts
#'     for the number of predictors. Can also extract multiple R-squared.
#'   \item \strong{P-value}: Overall model significance from the F-statistic
#'     in the ANOVA table. This tests whether at least one predictor is
#'     significantly related to the outcome.
#' }
#'
#' For models without an F-statistic (e.g., some GLMs), the p-value will be NA.
#'
#' @return A data frame (tibble) with one row and two columns:
#'   \item{r_squared}{Adjusted (or multiple) R-squared value, optionally rounded}
#'   \item{p_value}{Overall model p-value from F-test}
#'
#' @author Xiang LI \email{lixiang117423@@foxmail.com}
#'
#' @importFrom tibble tibble
#' @importFrom stats anova
#'
#' @examples
#' # Simple linear regression
#' model <- lm(mpg ~ wt, data = mtcars)
#' get_lm_stats(model)
#'
#' # Multiple regression
#' model2 <- lm(mpg ~ wt + hp + cyl, data = mtcars)
#' get_lm_stats(model2)
#'
#' # Get multiple R-squared instead of adjusted
#' get_lm_stats(model, r_squared_type = "multiple")
#'
#' # No rounding
#' get_lm_stats(model, digits = NULL)
#'
#' # More decimal places
#' get_lm_stats(model, digits = 5)
#'
#' # Use in a plot annotation
#' library(ggplot2)
#' stats <- get_lm_stats(model)
#'
#' ggplot(mtcars, aes(x = wt, y = mpg)) +
#'   geom_point() +
#'   geom_smooth(method = "lm", se = FALSE) +
#'   annotate(
#'     "text",
#'     x = 4.5, y = 30,
#'     label = sprintf(
#'       "R² = %.3f\np = %.4f",
#'       stats$r_squared,
#'       stats$p_value
#'     )
#'   )
#'
#' # Compare multiple models
#' library(dplyr)
#' library(purrr)
#'
#' models <- list(
#'   model1 = lm(mpg ~ wt, data = mtcars),
#'   model2 = lm(mpg ~ wt + hp, data = mtcars),
#'   model3 = lm(mpg ~ wt + hp + cyl, data = mtcars)
#' )
#'
#' map_dfr(models, get_lm_stats, .id = "model")
#'
#' @author Xiang LI \email{lixiang117423@@foxmail.com}
#' @export
get_lm_stats <- function(model,
                         r_squared_type = c("adjusted", "multiple"),
                         digits = 3) {
  # Input validation
  validate_lm_model(model)
  r_squared_type <- match.arg(r_squared_type)
  validate_digits(digits)

  # Extract statistics
  r_squared <- extract_r_squared(model, r_squared_type, digits)
  p_value <- extract_model_p_value(model)

  # Return as tibble
  tibble::tibble(
    r_squared = r_squared,
    p_value = p_value
  )
}


#' Validate Linear Model Object
#'
#' @param model Object to validate
#' @keywords internal
#' @noRd
validate_lm_model <- function(model) {
  # Check if model has required methods
  if (!inherits(model, "lm")) {
    warning(
      "Model is not of class 'lm'. ",
      "Results may be unreliable for class: ",
      paste(class(model), collapse = ", "),
      call. = FALSE
    )
  }

  # Check if summary method exists
  if (!hasMethod("summary", class(model))) {
    stop(
      "Model object does not have a summary() method",
      call. = FALSE
    )
  }

  # Check if anova method exists
  if (!hasMethod("anova", class(model))) {
    stop(
      "Model object does not have an anova() method",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Validate Digits Parameter
#'
#' @param digits Value to validate
#' @keywords internal
#' @noRd
validate_digits <- function(digits) {
  if (!is.null(digits)) {
    if (!is.numeric(digits) || length(digits) != 1) {
      stop(
        "digits must be a single numeric value or NULL",
        call. = FALSE
      )
    }

    if (digits < 0 || digits != round(digits)) {
      stop(
        "digits must be a non-negative integer",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}


#' Extract R-squared from Model
#'
#' @param model Fitted model object
#' @param type Type of R-squared to extract
#' @param digits Number of digits to round to
#' @return R-squared value
#' @keywords internal
#' @noRd
extract_r_squared <- function(model, type, digits) {
  model_summary <- summary(model)

  # Extract appropriate R-squared
  r_squared <- switch(
    type,
    adjusted = model_summary$adj.r.squared,
    multiple = model_summary$r.squared
  )

  # Handle missing R-squared
  if (is.null(r_squared)) {
    warning(
      "R-squared not available for this model type",
      call. = FALSE
    )
    return(NA_real_)
  }

  # Round if requested
  if (!is.null(digits)) {
    r_squared <- round(r_squared, digits)
  }

  return(r_squared)
}


#' Extract Model P-value from ANOVA
#'
#' @param model Fitted model object
#' @return P-value from F-test
#' @keywords internal
#' @noRd
extract_model_p_value <- function(model) {
  # Try to get ANOVA table
  anova_result <- tryCatch(
    stats::anova(model),
    error = function(e) {
      warning(
        "Could not extract ANOVA table: ",
        e$message,
        call. = FALSE
      )
      return(NULL)
    }
  )

  # Return NA if ANOVA failed
  if (is.null(anova_result)) {
    return(NA_real_)
  }

  # Extract p-value from first row (model effect)
  p_col_name <- get_p_value_column(anova_result)

  if (is.null(p_col_name)) {
    warning(
      "P-value column not found in ANOVA table",
      call. = FALSE
    )
    return(NA_real_)
  }

  p_value <- anova_result[[p_col_name]][1]

  # Handle missing p-value
  if (is.null(p_value) || is.na(p_value)) {
    warning(
      "P-value not available for this model",
      call. = FALSE
    )
    return(NA_real_)
  }

  return(p_value)
}


#' Get P-value Column Name from ANOVA Table
#'
#' @param anova_result ANOVA table
#' @return Name of p-value column or NULL
#' @keywords internal
#' @noRd
get_p_value_column <- function(anova_result) {
  # Common p-value column names
  possible_names <- c(
    "Pr(>F)",
    "Pr(>Chi)",
    "Pr(>|t|)",
    "p.value",
    "P.value",
    "pvalue",
    "Pvalue"
  )

  col_names <- names(anova_result)

  for (name in possible_names) {
    if (name %in% col_names) {
      return(name)
    }
  }

  return(NULL)
}


#' Format Linear Model Statistics as Text
#'
#' @description
#' Format R-squared and p-value from a linear model as text suitable for
#' plot annotations or reports.
#'
#' @param model A fitted linear model object from \code{lm()}.
#' @param r_squared_type Character string. Type of R-squared ("adjusted" or "multiple").
#' @param r_squared_digits Integer. Decimal places for R-squared. Default is 3.
#' @param p_digits Integer. Decimal places for p-value. Default is 4.
#' @param format Character string. Output format:
#'   \itemize{
#'     \item "expression" (default): R expression for ggplot2 annotations
#'     \item "text": Plain text with newline separator
#'     \item "markdown": Markdown formatted text
#'   }
#'
#' @return Formatted text string or expression
#'
#' @examples
#' model <- lm(mpg ~ wt, data = mtcars)
#'
#' # For ggplot2 annotation
#' format_lm_stats(model)
#'
#' # Plain text
#' format_lm_stats(model, format = "text")
#'
#' # Markdown
#' format_lm_stats(model, format = "markdown")
#'
#' @author Xiang LI \email{lixiang117423@@foxmail.com}
#' @export
format_lm_stats <- function(model,
                            r_squared_type = c("adjusted", "multiple"),
                            r_squared_digits = 3,
                            p_digits = 4,
                            format = c("expression", "text", "markdown")) {
  r_squared_type <- match.arg(r_squared_type)
  format <- match.arg(format)

  stats <- get_lm_stats(model, r_squared_type, r_squared_digits)

  r2_symbol <- if (r_squared_type == "adjusted") "R²" else "R²"
  r2_value <- sprintf(paste0("%.", r_squared_digits, "f"), stats$r_squared)
  p_value <- sprintf(paste0("%.", p_digits, "f"), stats$p_value)

  switch(
    format,
    expression = bquote(
      .(r2_symbol) ~ "=" ~ .(r2_value) * "," ~ italic(p) ~ "=" ~ .(p_value)
    ),
    text = sprintf("%s = %s\np = %s", r2_symbol, r2_value, p_value),
    markdown = sprintf("**%s** = %s  \n*p* = %s", r2_symbol, r2_value, p_value)
  )
}


#' @rdname get_lm_stats
#' @author Xiang LI \email{lixiang117423@@foxmail.com}
#' @export
extract_lm_stats <- get_lm_stats
