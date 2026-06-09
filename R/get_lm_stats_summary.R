#' Extract Key Statistics from Linear Model Fits
#'
#' Extract and format essential statistics from linear model (lm) objects,
#' including R-squared, adjusted R-squared, p-values, and model fit metrics.
#' This is a simplified alternative to \code{\link[stats]{summary.lm}} for
#' quick model assessment and reporting.
#'
#' @param fit A fitted linear model object of class "lm" or "glm"
#' @param digits Number of decimal places to round the statistics (default: 3)
#'
#' @return A data frame with columns:
#'   \code{r_squared}: Multiple R-squared value
#'
#'   \code{adj_r_squared}: Adjusted R-squared value
#'
#'   \code{p_value}: P-value for the overall model F-test
#'
#'   \code{f_statistic}: F-statistic value
#'
#'   \code{df_model}: Degrees of freedom for the model (numerator)
#'
#'   \code{df_residual}: Degrees of freedom for residuals (denominator)
#'
#'   \code{residual_se}: Residual standard error
#'
#'   \code{n_obs}: Number of observations
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export?
#'
#' @examples
#' \dontrun{
#' library(bioRtools)
#'
#' # Fit a linear model
#' fit <- lm(mpg ~ hp + wt + qsec, data = mtcars)
#'
#' # Get summary statistics
#' stats <- get_lm_stats_summary(fit)
#' print(stats)
#'
#' # Use with custom precision
#' stats_detailed <- get_lm_stats_summary(fit, digits = 5)
#' print(stats_detailed)
#'
#' # Example with multiple models
#' models <- list(
#'   model1 = lm(mpg ~ hp, data = mtcars),
#'   model2 = lm(mpg ~ hp + wt, data = mtcars),
#'   model3 = lm(mpg ~ hp + wt + qsec, data = mtcars)
#' )
#'
#' # Compare models
#' model_stats <- dplyr::bind_rows(
#'   lapply(names(models), function(name) {
#'     stats <- get_lm_stats_summary(models[[name]])
#'     stats$model <- name
#'     stats
#'   })
#' )
#'
#' print(model_stats)
#'
#' # Format for publication
#' cat(sprintf(
#'   "Model: R² = %.3f, Adj. R² = %.3f, p = %.2e\n",
#'   stats$r_squared, stats$adj_r_squared, stats$p_value
#' ))
#' }
get_lm_stats_summary <- function(fit, digits = 3) {

  # Input validation
  if (missing(fit) || is.null(fit)) {
    stop("'fit' is required")
  }

  if (!inherits(fit, c("lm", "glm"))) {
    stop("'fit' must be an object of class 'lm' or 'glm'")
  }

  if (!is.numeric(digits) || length(digits) != 1 || digits < 0) {
    stop("'digits' must be a non-negative integer")
  }

  digits <- as.integer(digits)

  # Extract model summary
  tryCatch({
    model_summary <- summary(fit)
  }, error = function(e) {
    stop(paste("Failed to extract model summary:", e$message))
  })

  # Extract R-squared values
  r_squared <- model_summary$r.squared

  adj_r_squared <- if (!is.null(model_summary$adj.r.squared)) {
    model_summary$adj.r.squared
  } else {
    NA
  }

  # Extract F-test information
  if (!is.null(model_summary$fstatistic)) {
    f_stat <- model_summary$fstatistic[1]
    df1 <- model_summary$fstatistic[2]
    df2 <- model_summary$fstatistic[3]

    # Calculate p-value from F-statistic
    p_value <- stats::pf(f_stat, df1, df2, lower.tail = FALSE)
  } else {
    # For models without F-test (e.g., some glm models)
    f_stat <- NA
    df1 <- NA
    df2 <- NA
    p_value <- NA
  }

  # Extract residual standard error (sigma)
  sigma <- model_summary$sigma

  # Extract number of observations
  n_observations <- model_summary$df[2] + model_summary$df[1] + 1

  # Handle edge cases for glm objects
  if (inherits(fit, "glm")) {
    # For glm, use deviance-based pseudo R-squared if available
    if (is.null(r_squared)) {
      # Calculate McFadden's pseudo R-squared
      null_deviance <- fit$null.deviance
      residual_deviance <- fit$deviance
      r_squared <- 1 - (residual_deviance / null_deviance)
      adj_r_squared <- NA  # No adjusted R-squared for pseudo R-squared
    }
  }

  # Create output data frame
  result <- data.frame(
    r_squared = round(r_squared, digits),
    adj_r_squared = round(adj_r_squared, digits),
    p_value = if (!is.na(p_value)) format(p_value, scientific = TRUE, digits = 3) else NA_character_,
    f_statistic = round(f_stat, digits),
    df1 = as.integer(df1),
    df2 = as.integer(df2),
    sigma = round(sigma, digits),
    n_observations = as.integer(n_observations),
    stringsAsFactors = FALSE
  )

  # Add informative column names
  names(result) <- c(
    "r_squared",
    "adj_r_squared",
    "p_value",
    "f_statistic",
    "df_model",
    "df_residual",
    "residual_se",
    "n_obs"
  )

  result
}
