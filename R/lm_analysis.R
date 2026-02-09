#' Perform Group-wise Linear Regression Analysis with Visualization
#'
#' This function performs linear regression analysis between two continuous variables
#' within each level of a grouping variable. It calculates regression statistics
#' (R-squared, p-values, coefficients) for each group separately and creates an
#' integrated scatter plot with fitted regression lines. This is particularly
#' useful for comparing relationships across different experimental conditions,
#' species, treatments, or other categorical groupings.
#'
#' @param data Data frame containing all variables needed for the analysis.
#'   Must include the variables specified in x, y, group, and color parameters
#' @param x Character string specifying the column name for the independent
#'   variable (predictor). This variable will be plotted on the X-axis and
#'   used as the predictor in linear regression models
#' @param y Character string specifying the column name for the dependent
#'   variable (response). This variable will be plotted on the Y-axis and
#'   used as the response in linear regression models
#' @param group Character string specifying the column name for the grouping
#'   variable. Separate linear regression models will be fitted for each
#'   unique level of this variable
#' @param color Character string specifying the column name for point colors
#'   in the scatter plot. Often the same as the group variable but can be
#'   different to show additional data dimensions
#' @param alpha Significance level for regression p-values (default: 0.05).
#'   Used to determine statistical significance in the results summary
#' @param conf.level Confidence level for regression confidence intervals
#'   (default: 0.95). Used in both statistical calculations and plot confidence bands
#' @param method Regression method for line fitting (default: "lm" for linear model).
#'   Can also be "glm", "gam", "loess", or other methods supported by ggplot2::geom_smooth
#' @param se Logical indicating whether to display confidence intervals around
#'   regression lines in the plot (default: TRUE)
#'
#' @return A named list containing two elements:
#'   \describe{
#'     \item{\code{regression_results}}{Data frame with group-wise regression statistics:
#'       \itemize{
#'         \item \code{group}: Group levels from the grouping variable
#'         \item \code{n_observations}: Number of observations in each group
#'         \item \code{r_squared}: Coefficient of determination (R²)
#'         \item \code{adj_r_squared}: Adjusted R-squared
#'         \item \code{slope}: Regression slope (effect of x on y)
#'         \item \code{intercept}: Y-intercept of regression line
#'         \item \code{slope_se}: Standard error of slope estimate
#'         \item \code{slope_pvalue}: P-value for slope significance test
#'         \item \code{overall_pvalue}: Overall model F-test p-value
#'         \item \code{significance}: Significance classification based on alpha
#'         \item \code{slope_ci_lower}: Lower bound of slope confidence interval
#'         \item \code{slope_ci_upper}: Upper bound of slope confidence interval
#'       }}
#'     \item{\code{scatter_plot}}{ggplot2 object showing scatter plot with:
#'       \itemize{
#'         \item Individual data points colored by group
#'         \item Fitted regression lines for each group
#'         \item Confidence intervals (optional)
#'         \item Professional theme and labeling
#'       }}
#'   }
#'
#' @details
#' For each group, the function fits a linear regression model: y ~ x and extracts:
#' \itemize{
#'   \item \strong{R-squared}: Proportion of variance in y explained by x
#'   \item \strong{Adjusted R-squared}: R² adjusted for degrees of freedom
#'   \item \strong{Slope}: Rate of change in y per unit change in x
#'   \item \strong{P-values}: Statistical significance of slope and overall model
#'   \item \strong{Confidence intervals}: Uncertainty bounds for slope estimates
#' }
#'
#' The visualization combines all groups in a single plot for easy comparison
#' of regression relationships across groups.
#'
#' @note
#' \itemize{
#'   \item Both x and y variables should be continuous (numeric)
#'   \item Each group should have at least 3 observations for reliable regression
#'   \item The function assumes linear relationships; check residual plots for model validity
#'   \item Groups with insufficient data will generate warnings
#'   \item Consider data transformation if relationships appear non-linear
#' }
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[stats]{lm}} for underlying linear regression implementation
#'   \item \code{\link{cor_analysis}} for correlation analysis
#'   \item \code{\link{anova_posthoc}} for group comparisons
#' }
#'
#' @export
#'
#' @examples
#' library(bioRtools)
#' library(dplyr)
#' library(ggplot2)
#'
#' # Example 1: Basic regression analysis with iris data
#' data(iris)
#'
#' # Analyze relationship between sepal length and width by species
#' results <- lm_analysis(
#'   data = iris,
#'   x = "Sepal.Length",
#'   y = "Sepal.Width",
#'   group = "Species",
#'   color = "Species"
#' )
#'
#' # View regression statistics
#' print("Regression results by species:")
#' print(results$regression_results)
#'
#' # Display the plot
#' print(results$scatter_plot)
#'
#' # Summary of significant relationships
#' significant_groups <- results$regression_results %>%
#'   filter(significance == "Significant")
#'
#' print(paste("Significant relationships found in", nrow(significant_groups), "groups"))
#'
lm_analysis <- function(data, x, y, group, color, alpha = 0.05, conf.level = 0.95,
                        method = "lm", se = TRUE) {
  # Input validation
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }

  if (nrow(data) == 0) {
    stop("'data' cannot be empty")
  }

  # Validate column names
  required_cols <- c(x, y, group, color)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing columns in data:", paste(missing_cols, collapse = ", ")))
  }

  # Validate parameter types
  if (!is.character(x) || length(x) != 1) {
    stop("'x' must be a single character string (column name)")
  }

  if (!is.character(y) || length(y) != 1) {
    stop("'y' must be a single character string (column name)")
  }

  if (!is.character(group) || length(group) != 1) {
    stop("'group' must be a single character string (column name)")
  }

  if (!is.character(color) || length(color) != 1) {
    stop("'color' must be a single character string (column name)")
  }

  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single number between 0 and 1")
  }

  if (!is.numeric(conf.level) || length(conf.level) != 1 || conf.level <= 0 || conf.level >= 1) {
    stop("'conf.level' must be a single number between 0 and 1")
  }

  if (!is.logical(se) || length(se) != 1) {
    stop("'se' must be a single logical value")
  }

  # Check if x and y are numeric
  if (!is.numeric(data[[x]])) {
    stop(paste("Variable", x, "must be numeric for linear regression"))
  }

  if (!is.numeric(data[[y]])) {
    stop(paste("Variable", y, "must be numeric for linear regression"))
  }

  # Remove rows with missing values in key variables
  clean_data <- data %>%
    dplyr::filter(!is.na(!!rlang::sym(x)), !is.na(!!rlang::sym(y)),
      !is.na(!!rlang::sym(group)), !is.na(!!rlang::sym(color)))

  if (nrow(clean_data) == 0) {
    stop("No complete cases remain after removing missing values")
  }

  if (nrow(clean_data) < nrow(data)) {
    warning(paste("Removed", nrow(data) - nrow(clean_data),
      "rows due to missing values"))
  }

  # Check group sizes
  group_sizes <- clean_data %>%
    dplyr::count(!!rlang::sym(group), name = "n_obs") %>%
    dplyr::arrange(n_obs)

  min_group_size <- min(group_sizes$n_obs)
  if (min_group_size < 3) {
    warning(paste("Some groups have fewer than 3 observations. Smallest group:",
      min_group_size, "observations. Results may be unreliable."))
  }

  # Check for groups with no variation in x or y
  variance_check <- clean_data %>%
    dplyr::group_by(!!rlang::sym(group)) %>%
    dplyr::summarise(
      x_var = var(!!rlang::sym(x), na.rm = TRUE),
      y_var = var(!!rlang::sym(y), na.rm = TRUE),
      .groups = "drop"
    )

  zero_var_groups <- variance_check %>%
    dplyr::filter(x_var == 0 | y_var == 0)

  if (nrow(zero_var_groups) > 0) {
    warning(paste("Groups with zero variance detected:",
      paste(zero_var_groups[[group]], collapse = ", "),
      ". These groups will have undefined regression statistics."))
  }

  # Define helper functions for regression statistics
  get_regression_stats <- function(fit, conf_level = 0.95) {
    if (is.null(fit) || !inherits(fit, "lm")) {
      return(data.frame(
        n_observations = NA, r_squared = NA, adj_r_squared = NA,
        slope = NA, intercept = NA, slope_se = NA, slope_pvalue = NA,
        overall_pvalue = NA, slope_ci_lower = NA, slope_ci_upper = NA
      ))
    }

    tryCatch(
      {
        fit_summary <- summary(fit)
        anova_result <- stats::anova(fit)
        confint_result <- confint(fit, level = conf_level)

        data.frame(
          n_observations = nobs(fit),
          r_squared = round(fit_summary$r.squared, 4),
          adj_r_squared = round(fit_summary$adj.r.squared, 4),
          slope = round(coef(fit)[2], 4),
          intercept = round(coef(fit)[1], 4),
          slope_se = round(fit_summary$coefficients[2, "Std. Error"], 4),
          slope_pvalue = fit_summary$coefficients[2, "Pr(>|t|)"],
          overall_pvalue = anova_result$`Pr(>F)`[1],
          slope_ci_lower = round(confint_result[2, 1], 4),
          slope_ci_upper = round(confint_result[2, 2], 4)
        )
      },
      error = function(e) {
        warning(paste("Error extracting regression statistics:", e$message))
        data.frame(
          n_observations = NA, r_squared = NA, adj_r_squared = NA,
          slope = NA, intercept = NA, slope_se = NA, slope_pvalue = NA,
          overall_pvalue = NA, slope_ci_lower = NA, slope_ci_upper = NA
        )
      })
  }

  # Perform group-wise regression analysis
  regression_results <- clean_data %>%
    dplyr::group_by(!!rlang::sym(group)) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      # Fit linear models
      model = purrr::map(data, ~ {
        tryCatch(
          {
            lm(as.formula(paste(y, "~", x)), data = .x)
          },
          error = function(e) {
            warning(paste("Model fitting failed for group",
              unique(.x[[group]]), ":", e$message))
            NULL
          })
      }),
      # Extract statistics
      stats = purrr::map(model, ~ get_regression_stats(.x, conf.level))
    ) %>%
    dplyr::select(-data, -model) %>%
    tidyr::unnest(cols = stats) %>%
    dplyr::mutate(
      # Add significance classification
      significance = dplyr::case_when(
        is.na(overall_pvalue) ~ "Unable to determine",
        overall_pvalue < alpha ~ "Significant",
        TRUE ~ "Not significant"
      ),
      # Format p-values for display
      overall_pvalue_formatted = dplyr::case_when(
        is.na(overall_pvalue) ~ "NA",
        overall_pvalue < 0.001 ~ "< 0.001",
        TRUE ~ sprintf("%.3f", overall_pvalue)
      ),
      slope_pvalue_formatted = dplyr::case_when(
        is.na(slope_pvalue) ~ "NA",
        slope_pvalue < 0.001 ~ "< 0.001",
        TRUE ~ sprintf("%.3f", slope_pvalue)
      )
    ) %>%
    dplyr::arrange(overall_pvalue)

  # Create enhanced scatter plot
  tryCatch(
    {
      scatter_plot <- clean_data %>%
        ggplot2::ggplot(ggplot2::aes(x = !!rlang::sym(x), y = !!rlang::sym(y),
          color = !!rlang::sym(color))) +
        ggplot2::geom_point(alpha = 0.7, size = 2) +
        ggplot2::geom_smooth(method = method, se = se, level = conf.level) +
        ggplot2::labs(
          title = paste("Linear Regression Analysis by", group),
          subtitle = paste("Method:", method, "| Confidence level:",
            paste0(round(conf.level * 100, 1), "%")),
          x = tools::toTitleCase(gsub("[._]", " ", x)),
          y = tools::toTitleCase(gsub("[._]", " ", y)),
          color = tools::toTitleCase(gsub("[._]", " ", color))
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 14, face = "bold"),
          plot.subtitle = ggplot2::element_text(size = 11),
          axis.title = ggplot2::element_text(size = 12),
          legend.title = ggplot2::element_text(size = 11),
          legend.text = ggplot2::element_text(size = 10),
          panel.grid.minor = ggplot2::element_blank()
        )
    },
    error = function(e) {
      warning(paste("Error creating plot:", e$message))
      scatter_plot <- ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5,
          label = "Plot creation failed", size = 6) +
        ggplot2::theme_void()
    })

  # Add summary statistics
  n_total_obs <- nrow(clean_data)
  n_groups <- nrow(regression_results)
  n_significant <- sum(regression_results$significance == "Significant", na.rm = TRUE)

  # Print summary if interactive
  if (interactive()) {
    cat("Linear Regression Analysis Summary:\n")
    cat("==================================\n")
    cat("Total observations:", n_total_obs, "\n")
    cat("Number of groups:", n_groups, "\n")
    cat("Groups with significant relationships:", n_significant,
      paste0("(", round(100 * n_significant / n_groups, 1), "%)"), "\n")
    cat("Significance level:", alpha, "\n")
    cat("Confidence level:", conf.level, "\n\n")

    if (n_significant > 0) {
      cat("Significant relationships:\n")
      sig_results <- regression_results %>%
        dplyr::filter(significance == "Significant") %>%
        dplyr::select(!!rlang::sym(group), slope, r_squared, overall_pvalue_formatted)

      for (i in 1:nrow(sig_results)) {
        cat(sprintf("  %s: slope = %.3f, R² = %.3f, p = %s\n",
          sig_results[[group]][i],
          sig_results$slope[i],
          sig_results$r_squared[i],
          sig_results$overall_pvalue_formatted[i]))
      }
    }
    cat("\n")
  }

  # Prepare final result
  result <- list(
    regression_results = regression_results,
    scatter_plot = scatter_plot
  )

  # Add metadata as attributes
  attr(result, "analysis_type") <- "group_wise_linear_regression"
  attr(result, "x_variable") <- x
  attr(result, "y_variable") <- y
  attr(result, "group_variable") <- group
  attr(result, "color_variable") <- color
  attr(result, "n_observations") <- n_total_obs
  attr(result, "n_groups") <- n_groups
  attr(result, "n_significant") <- n_significant
  attr(result, "alpha_level") <- alpha
  attr(result, "confidence_level") <- conf.level
  attr(result, "regression_method") <- method

  return(result)
}
