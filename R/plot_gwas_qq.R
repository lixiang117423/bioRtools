#' Create QQ Plot for GWAS P-values
#'
#' This function creates a Quantile-Quantile (QQ) plot for visualizing the
#' distribution of p-values from Genome-Wide Association Studies (GWAS) or
#' other statistical analyses. It helps assess goodness-of-fit and identify
#' potential population stratification or systematic biases.
#'
#' @param df A data frame containing p-values
#' @param p_column Character string specifying the column name containing p-values
#'   (default: "P")
#' @param title Character string for the plot title
#'   (default: "QQ Plot")
#' @param point_color Color for points (default: "#2874A6")
#' @param point_size Numeric value for point size (default: 0.8)
#' @param point_alpha Numeric value for point transparency (default: 0.7)
#' @param ribbon_fill Color for confidence interval ribbon (default: "#AED6F1")
#' @param ribbon_alpha Numeric value for ribbon transparency (default: 0.4)
#' @param line_color Color for the reference line (default: "red")
#' @param line_width Numeric value for line width (default: 0.8)
#' @param line_type Character string for line type (default: "dashed")
#' @param base_size Numeric value for base font size (default: 13)
#'
#' @return A ggplot object representing the QQ plot. The function invisibly
#'   returns the calculated genomic inflation factor (lambda) as an attribute.
#'
#' @details
#' The QQ plot compares the observed distribution of p-values against the
#' expected uniform distribution under the null hypothesis. Key features:
#'
#' \itemize{
#'   \item \strong{Reference line}: Diagonal line representing perfect
#'     agreement between observed and expected p-values
#'   \item \strong{Confidence interval}: Shaded region showing the 95%
#'     confidence bounds calculated using the beta distribution
#'   \item \strong{Lambda (genomic inflation factor)}: Calculated as
#'     median(chisq_observed) / chisq_expected, where values > 1 indicate
#'     potential population stratification or systematic bias
#' }
#'
#' Lambda is calculated but not displayed on the plot by default. You can
#' retrieve it using \code{attr(plot, "lambda")}.
#'
#' @note
#' \itemize{
#'   \item P-values should be between 0 and 1
#'   \item Missing values (NA) in the p-value column will be removed
#'   \item The function uses -log10 transformation for both axes
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(bioRtools)
#'
#' # Generate example p-value data
#' set.seed(123)
#' n <- 1000
#' df <- data.frame(
#'   P = c(runif(900), runif(100, 0, 0.01))  # 90% null, 10% significant
#' )
#'
#' # Create basic QQ plot
#' p1 <- plot_gwas_qq(df)
#' print(p1)
#'
#' # Get lambda value
#' lambda <- attr(p1, "lambda")
#' print(paste0("Genomic inflation factor: ", round(lambda, 3)))
#'
#' # Custom QQ plot with specific column name
#' df2 <- data.frame(pvalue = runif(500))
#' p2 <- plot_gwas_qq(df2, p_column = "pvalue", title = "My Custom QQ Plot")
#' print(p2)
#'
#' # QQ plot with custom styling
#' p3 <- plot_gwas_qq(
#'   df,
#'   title = "GWAS Results",
#'   point_color = "blue",
#'   point_size = 1,
#'   ribbon_fill = "lightgray"
#' )
#' print(p3)
#' }
#'
plot_gwas_qq <- function(df,
                         p_column = "P",
                         title = "QQ Plot",
                         point_color = "#2874A6",
                         point_size = 0.8,
                         point_alpha = 0.7,
                         ribbon_fill = "#AED6F1",
                         ribbon_alpha = 0.4,
                         line_color = "red",
                         line_width = 0.8,
                         line_type = "dashed",
                         base_size = 13) {

  # Input validation
  if (!is.data.frame(df)) {
    stop("'df' must be a data frame")
  }

  if (!is.character(p_column) || length(p_column) != 1) {
    stop("'p_column' must be a single character string")
  }

  if (!p_column %in% names(df)) {
    stop(paste0("Column '", p_column, "' not found in data frame"))
  }

  if (!is.character(title) || length(title) != 1) {
    stop("'title' must be a single character string")
  }

  # Extract p-values and handle missing values
  p_values <- df[[p_column]]
  p_values <- p_values[!is.na(p_values)]

  # Check if we have enough valid p-values
  if (length(p_values) < 2) {
    stop("At least 2 valid p-values are required")
  }

  # Check p-value range
  if (any(p_values < 0 | p_values > 1)) {
    stop("P-values must be between 0 and 1")
  }

  # Calculate QQ statistics
  n <- length(p_values)
  obs <- sort(p_values)
  exp <- (1:n) / (n + 1)

  # Create QQ data frame
  qq_df <- data.frame(
    expected = -log10(exp),
    observed = -log10(obs),
    ci_upper = -log10(qbeta(0.025, 1:n, n - 1:n + 1)),
    ci_lower = -log10(qbeta(0.975, 1:n, n - 1:n + 1))
  )

  # Calculate genomic inflation factor (lambda)
  lambda <- median(qchisq(1 - obs, 1)) / qchisq(0.5, 1)

  # Create the plot
  p <- ggplot2::ggplot(qq_df, ggplot2::aes(x = expected, y = observed)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
      fill = ribbon_fill,
      alpha = ribbon_alpha
    ) +
    ggplot2::geom_abline(
      intercept = 0,
      slope = 1,
      color = line_color,
      linewidth = line_width,
      linetype = line_type
    ) +
    ggplot2::geom_point(
      color = point_color,
      size = point_size,
      alpha = point_alpha
    ) +
    ggplot2::labs(
      title = title,
      x = expression(Expected ~ -log[10](italic(p))),
      y = expression(Observed ~ -log[10](italic(p)))
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )

  # Attach lambda as an attribute
  attr(p, "lambda") <- lambda
  attr(p, "n_tests") <- n

  p
}
