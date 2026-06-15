#' Create QQ Plot for P-values or General Numeric Data
#'
#' @description
#' Creates a Quantile-Quantile (QQ) plot with two modes:
#'
#' \itemize{
#'   \item \code{type = "pvalue"}: GWAS p-value QQ plot. Compares observed
#'     p-values against the expected uniform distribution under the null,
#'     using -log10 transformation. Useful for assessing goodness-of-fit
#'     and detecting population stratification.
#'   \item \code{type = "normal"}: Normal QQ plot. Compares sample values
#'     against theoretical normal quantiles. Useful for checking normality
#'     of any numeric variable (e.g., transformed phenotypes, expression
#'     values).
#' }
#'
#' When \code{group_column} is set, QQ statistics are computed per group and
#' the returned plot preserves the group column in every layer, so you can
#' append \code{+ facet_wrap(...)} directly.
#'
#' @param df A data frame.
#' @param value_column Character string specifying the column to plot.
#'   Default: "P" for \code{type = "pvalue"}, "value" for \code{type = "normal"}.
#' @param group_column Optional character string naming a grouping column.
#'   When set, QQ stats are computed per group and the group column is
#'   preserved for downstream faceting (default: NULL).
#' @param type Character string specifying the QQ plot type. One of
#'   \code{"pvalue"} (default) or \code{"normal"}.
#' @param title Plot title (default: "QQ Plot").
#' @param point_color Color for points (default: "#2874A6").
#' @param point_size Numeric value for point size (default: 0.8).
#' @param point_alpha Numeric value for point transparency (default: 0.7).
#' @param ribbon_fill Color for confidence interval ribbon (default: "#AED6F1").
#' @param ribbon_alpha Numeric value for ribbon transparency (default: 0.4).
#' @param line_color Color for the reference line (default: "red").
#' @param line_width Numeric value for line width (default: 0.8).
#' @param line_type Character string for line type (default: "dashed").
#' @param base_size Numeric value for base font size (default: 13).
#' @param ci_level Numeric value for confidence interval level (default: 0.95).
#'
#' @return A ggplot object.
#'   For \code{type = "pvalue"} (single group): the genomic inflation factor
#'   (lambda) is attached as attribute \code{"lambda"}, number of tests as
#'   attribute \code{"n_tests"}.
#'   For \code{type = "normal"} (single group): number of observations is
#'   attached as attribute \code{"n_obs"}.
#'   For grouped input: attribute \code{"group_stats"} holds a data frame of
#'   per-group statistics (group value, n, and lambda for pvalue mode).
#'
#' @details
#' For \code{type = "pvalue"}, values must be in [0, 1]. The reference line
#' is the identity line in -log10 space. The confidence ribbon is computed
#' from the beta distribution. Lambda (genomic inflation factor) is
#' \code{median(chisq_observed) / qchisq(0.5, 1)}; values > 1 suggest
#' population stratification.
#'
#' For \code{type = "normal"}, values can be any numeric. The reference line
#' passes through the sample mean and uses the sample standard deviation as
#' slope. The confidence ribbon is based on the standard error of order
#' statistics under normality: \code{SE(x_(i)) ~ sigma * sqrt(p*(1-p)/n) / phi(z_i)},
#' so intervals widen at the tails.
#'
#' Missing values and non-finite values are removed before plotting.
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' \dontrun{
#' library(bioRtools)
#'
#' # p-value QQ plot (GWAS)
#' set.seed(123)
#' df_p <- data.frame(P = c(runif(900), runif(100, 0, 0.01)))
#' p1 <- plot_qq(df_p, type = "pvalue")
#' print(p1)
#' attr(p1, "lambda")
#'
#' # Normal QQ plot for any numeric data
#' df_v <- data.frame(value = rnorm(500))
#' plot_qq(df_v, type = "normal")
#'
#' # Compare phenotype transformations side-by-side via facet
#' phen <- data.frame(raw = rgamma(300, 2))
#' phen$log <- log(phen$raw)
#' phen$rint <- qnorm((rank(phen$raw) - 0.5) / nrow(phen))
#' phen_long <- tidyr::pivot_longer(phen, dplyr::everything())
#' plot_qq(phen_long, value_column = "value", group_column = "name",
#'         type = "normal") +
#'   ggplot2::facet_wrap(~ name, scales = "free")
#' }
#'
plot_qq <- function(df,
                    value_column = NULL,
                    group_column = NULL,
                    type = c("pvalue", "normal"),
                    title = "QQ Plot",
                    point_color = "#2874A6",
                    point_size = 0.8,
                    point_alpha = 0.7,
                    ribbon_fill = "#AED6F1",
                    ribbon_alpha = 0.4,
                    line_color = "red",
                    line_width = 0.8,
                    line_type = "dashed",
                    base_size = 13,
                    ci_level = 0.95) {

  type <- match.arg(type)

  if (!is.data.frame(df)) {
    stop("'df' must be a data frame")
  }

  if (is.null(value_column)) {
    value_column <- if (type == "pvalue") "P" else "value"
  }

  if (!is.character(value_column) || length(value_column) != 1) {
    stop("'value_column' must be a single character string")
  }
  if (!value_column %in% names(df)) {
    stop(paste0("Column '", value_column, "' not found in data frame"))
  }
  if (!is.numeric(df[[value_column]])) {
    stop(paste0("Column '", value_column, "' must be numeric"))
  }

  if (!is.null(group_column)) {
    if (!is.character(group_column) || length(group_column) != 1) {
      stop("'group_column' must be a single character string")
    }
    if (!group_column %in% names(df)) {
      stop(paste0("Column '", group_column, "' not found in data frame"))
    }
  }

  if (!is.character(title) || length(title) != 1) {
    stop("'title' must be a single character string")
  }
  if (!is.numeric(ci_level) || length(ci_level) != 1 ||
      ci_level <= 0 || ci_level >= 1) {
    stop("'ci_level' must be a single numeric value strictly between 0 and 1")
  }

  is_pvalue <- type == "pvalue"
  x_lab <- if (is_pvalue) {
    expression(Expected ~ -log[10](italic(p)))
  } else {
    "Theoretical Quantiles"
  }
  y_lab <- if (is_pvalue) {
    expression(Observed ~ -log[10](italic(p)))
  } else {
    "Sample Quantiles"
  }

  if (is.null(group_column)) {
    values <- df[[value_column]]
    values <- values[!is.na(values)]
    if (length(values) < 2) {
      stop("At least 2 valid values are required")
    }

    qq <- qq_compute(values, type, ci_level)
    line_df <- data.frame(intercept = qq$intercept, slope = qq$slope)

    p <- qq_build_single(
      qq$data, line_df,
      title = title, x_lab = x_lab, y_lab = y_lab,
      point_color = point_color, point_size = point_size, point_alpha = point_alpha,
      ribbon_fill = ribbon_fill, ribbon_alpha = ribbon_alpha,
      line_color = line_color, line_width = line_width, line_type = line_type,
      base_size = base_size
    )

    if (is_pvalue) {
      attr(p, "lambda") <- qq$lambda
      attr(p, "n_tests") <- qq$n
    } else {
      attr(p, "n_obs") <- qq$n
    }

    return(p)
  }

  # Grouped path: compute per group, preserve group column for faceting
  group_values <- df[[group_column]]
  if (any(is.na(group_values))) {
    df <- df[!is.na(group_values), ]
    group_values <- df[[group_column]]
  }

  splits <- split(df[[value_column]], group_values)
  splits <- lapply(splits, function(v) v[!is.na(v)])
  if (any(lengths(splits) < 2)) {
    stop("Each group must contain at least 2 valid values")
  }

  qq_list <- lapply(splits, qq_compute, type = type, ci_level = ci_level)
  group_names <- names(splits)

  qq_df <- do.call(
    rbind,
    lapply(seq_along(qq_list), function(i) {
      cbind(qq_list[[i]]$data, setNames(list(group_names[[i]]), group_column))
    })
  )

  line_df <- data.frame(
    setNames(list(group_names), group_column),
    intercept = vapply(qq_list, function(x) x$intercept, numeric(1)),
    slope = vapply(qq_list, function(x) x$slope, numeric(1)),
    check.names = FALSE
  )

  group_stats <- data.frame(
    setNames(list(group_names), group_column),
    n = vapply(qq_list, function(x) x$n, integer(1)),
    check.names = FALSE
  )
  if (is_pvalue) {
    group_stats$lambda <- vapply(qq_list, function(x) x$lambda, numeric(1))
  }

  p <- qq_build_grouped(
    qq_df, line_df, group_column = group_column,
    title = title, x_lab = x_lab, y_lab = y_lab,
    point_color = point_color, point_size = point_size, point_alpha = point_alpha,
    ribbon_fill = ribbon_fill, ribbon_alpha = ribbon_alpha,
    line_color = line_color, line_width = line_width, line_type = line_type,
    base_size = base_size
  )

  attr(p, "group_stats") <- group_stats

  p
}

#' Compute QQ stats for a single group
#' @keywords internal
qq_compute <- function(values, type, ci_level) {

  alpha <- 1 - ci_level

  if (type == "pvalue") {
    if (any(values < 0 | values > 1)) {
      stop("For type = 'pvalue', values must be between 0 and 1")
    }

    n <- length(values)
    obs <- sort(values)
    i_seq <- seq_len(n)
    exp_p <- i_seq / (n + 1)

    data <- data.frame(
      expected = -log10(exp_p),
      observed = -log10(obs),
      ci_lower = -log10(stats::qbeta(1 - alpha / 2, i_seq, n - i_seq + 1)),
      ci_upper = -log10(stats::qbeta(alpha / 2, i_seq, n - i_seq + 1))
    )

    lambda <- median(stats::qchisq(1 - obs, 1)) / stats::qchisq(0.5, 1)

    list(data = data, slope = 1, intercept = 0, lambda = lambda, n = n)
  } else {
    finite <- is.finite(values)
    if (any(!finite)) values <- values[finite]
    if (length(values) < 2) stop("At least 2 finite values are required")

    n <- length(values)
    obs <- sort(values)

    p_pos <- (seq_len(n) - 0.5) / n
    exp_z <- stats::qnorm(p_pos)

    mu <- mean(obs)
    sigma <- stats::sd(obs)

    dz <- stats::dnorm(exp_z)
    se <- sigma * sqrt(p_pos * (1 - p_pos) / n) / dz
    z_crit <- stats::qnorm(1 - (1 - ci_level) / 2)
    expected_fit <- mu + sigma * exp_z

    data <- data.frame(
      expected = exp_z,
      observed = obs,
      ci_lower = expected_fit - z_crit * se,
      ci_upper = expected_fit + z_crit * se
    )

    list(data = data, slope = sigma, intercept = mu, n = n)
  }
}

#' Build single-group QQ plot
#' @keywords internal
qq_build_single <- function(qq_df, line_df, title, x_lab, y_lab,
                            point_color, point_size, point_alpha,
                            ribbon_fill, ribbon_alpha,
                            line_color, line_width, line_type, base_size) {

  ggplot2::ggplot(qq_df, ggplot2::aes(x = .data$expected, y = .data$observed)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$ci_lower, ymax = .data$ci_upper),
      fill = ribbon_fill, alpha = ribbon_alpha
    ) +
    ggplot2::geom_abline(
      intercept = line_df$intercept[1], slope = line_df$slope[1],
      color = line_color, linewidth = line_width, linetype = line_type
    ) +
    ggplot2::geom_point(
      color = point_color, size = point_size, alpha = point_alpha
    ) +
    ggplot2::labs(title = title, x = x_lab, y = y_lab) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )
}

#' Build grouped QQ plot (preserves group column for downstream faceting)
#' @keywords internal
qq_build_grouped <- function(qq_df, line_df, group_column, title, x_lab, y_lab,
                             point_color, point_size, point_alpha,
                             ribbon_fill, ribbon_alpha,
                             line_color, line_width, line_type, base_size) {

  ggplot2::ggplot(qq_df, ggplot2::aes(x = .data$expected, y = .data$observed)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$ci_lower, ymax = .data$ci_upper),
      fill = ribbon_fill, alpha = ribbon_alpha
    ) +
    ggplot2::geom_abline(
      data = line_df,
      mapping = ggplot2::aes(intercept = .data$intercept, slope = .data$slope),
      color = line_color, linewidth = line_width, linetype = line_type
    ) +
    ggplot2::geom_point(
      color = point_color, size = point_size, alpha = point_alpha
    ) +
    ggplot2::labs(title = title, x = x_lab, y = y_lab) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )
}
