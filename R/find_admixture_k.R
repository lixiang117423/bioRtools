#' Find the Best K from ADMIXTURE CV Errors
#'
#' Analyzes the cross-validation (CV) error curve of ADMIXTURE runs and
#' picks the best number of ancestral populations (K) with three classic
#' criteria: the Kneedle elbow (largest perpendicular distance from the
#' first-last chord in normalized coordinates), the maximal second-order
#' difference (local curvature), and first-K-that-drops-below relative-
#' change thresholds. Returns all intermediate tables plus two diagnostic
#' plots.
#'
#' @param data A data frame with the per-K CV results (the ADMIXTURE
#'   \code{cv_results.csv} layout), or a file path readable by
#'   \code{\link{read_data}} (e.g. \code{"cv_results.csv"}).
#' @param k_col Column name holding K values (default: \code{"K"}).
#' @param cv_col Column name holding the CV errors (default:
#'   \code{"CV_error"}).
#' @param rel_thresholds Numeric vector of relative-drop percentages; for
#'   each threshold the first K whose CV drop is below it is reported
#'   (default: \code{c(3, 2)}).
#' @param plot Logical; build the CV curve and Kneedle distance plots
#'   (default: \code{TRUE}).
#'
#' @return An object of class \code{admixture_k_results} (a list) with:
#'   \item{best_k}{K selected by the Kneedle elbow}
#'   \item{max_curvature_k}{K with maximal second-order difference}
#'   \item{threshold_k}{Named vector, first K below each relative-drop
#'   threshold (\code{NA} when never reached)}
#'   \item{kneedle_table}{Data frame: K, CV_error, distance}
#'   \item{diff_table}{Data frame: K_from, K_to, d1 (first difference)}
#'   \item{curvature_table}{Data frame: K, d2 (second difference)}
#'   \item{rel_change_table}{Data frame: K_from, K_to, rel_pct}
#'   \item{plots}{List with \code{cv_plot}, \code{distance_plot}, and the
#'   dual-axis \code{combined_plot} (CV error on the left, Kneedle distance
#'   on the right); \code{NULL} when \code{plot = FALSE}}
#'   A \code{print} method shows the key K picks.
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' # Typical ADMIXTURE cv_results.csv shape: CV drops fast then flattens
#' cv_df <- data.frame(
#'   K = 1:10,
#'   CV_error = c(0.58, 0.52, 0.50, 0.49, 0.485, 0.482, 0.481, 0.4805,
#'                0.4802, 0.4801)
#' )
#'
#' res <- find_admixture_k(cv_df)
#' res
#' res$plots$cv_plot
#' res$plots$combined_plot  # CV error + Kneedle distance on dual axes
#'
#' # Read the file directly
#' \dontrun{
#' res <- find_admixture_k("~/path/to/cv_results.csv")
#' res$plots$distance_plot
#' }
find_admixture_k <- function(data, k_col = "K", cv_col = "CV_error",
                             rel_thresholds = c(3, 2), plot = TRUE) {
  # ---- input handling ----
  if (is.character(data) && length(data) == 1L) {
    data <- read_data(data)
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame, or a file path readable by read_data()")
  }
  if (!is.character(k_col) || length(k_col) != 1L ||
      !k_col %in% names(data)) {
    stop("'k_col' must name a column of 'data' (default: 'K')")
  }
  if (!is.character(cv_col) || length(cv_col) != 1L ||
      !cv_col %in% names(data)) {
    stop("'cv_col' must name a column of 'data' (default: 'CV_error')")
  }

  K <- data[[k_col]]
  cv <- data[[cv_col]]
  if (!is.numeric(K) || !is.numeric(cv)) {
    stop("'k_col' and 'cv_col' must both be numeric")
  }
  if (length(K) != length(cv)) {
    stop("'k_col' and 'cv_col' must have the same length")
  }
  if (length(K) < 3) {
    stop("At least 3 K values are required (found ", length(K), ")")
  }
  if (!is.numeric(rel_thresholds) || length(rel_thresholds) == 0 ||
      any(!is.finite(rel_thresholds)) || any(rel_thresholds <= 0)) {
    stop("'rel_thresholds' must be a non-empty vector of positive numbers")
  }
  if (!is.logical(plot) || length(plot) != 1L || is.na(plot)) {
    stop("'plot' must be TRUE or FALSE")
  }

  # sort by K so differences follow the K axis
  ord <- order(K)
  K <- K[ord]
  cv <- cv[ord]
  if (any(duplicated(K))) {
    stop("Duplicate K values are not allowed")
  }
  if (diff(range(cv)) == 0) {
    stop("CV errors are constant across K; cannot locate an elbow")
  }

  n <- length(K)

  # ---- 1. Kneedle: perpendicular distance to the first-last chord ----
  K_norm <- (K - min(K)) / (max(K) - min(K))
  cv_norm <- (cv - min(cv)) / (max(cv) - min(cv))

  x1 <- K_norm[1]
  y1 <- cv_norm[1]
  x2 <- K_norm[n]
  y2 <- cv_norm[n]
  A <- y2 - y1
  B <- -(x2 - x1)
  C <- x2 * y1 - y2 * x1
  dist_to_line <- abs(A * K_norm + B * cv_norm + C) / sqrt(A^2 + B^2)

  kneedle_table <- tibble::tibble(
    K = K, CV_error = cv, distance = dist_to_line
  )
  best_k <- kneedle_table$K[which.max(kneedle_table$distance)]

  # ---- 2. First and second differences ----
  d1 <- diff(cv)
  diff_table <- tibble::tibble(K_from = K[-n], K_to = K[-1], d1 = d1)

  d2 <- diff(d1)
  curvature_table <- tibble::tibble(K = K[-c(1, n)], d2 = d2)
  max_curvature_k <- curvature_table$K[which.max(curvature_table$d2)]

  # ---- 3. Relative drop percentages + thresholds ----
  rel_change_table <- tibble::tibble(
    K_from = K[-n],
    K_to   = K[-1],
    rel_pct = abs(d1 / cv[-n]) * 100
  )

  threshold_k <- vapply(rel_thresholds, function(t) {
    hit <- rel_change_table$K_to[rel_change_table$rel_pct < t]
    if (length(hit) == 0) NA_real_ else hit[[1]]
  }, numeric(1))
  names(threshold_k) <- paste0("rel_pct_lt_", rel_thresholds)
  if (anyNA(threshold_k)) {
    warning("Relative drop never fell below ",
            paste(rel_thresholds[is.na(threshold_k)], collapse = ", "),
            "% within the tested K range")
  }

  # ---- 4. Diagnostic plots ----
  plots <- NULL
  if (plot) {
    cv_plot <- ggplot2::ggplot(kneedle_table,
                               ggplot2::aes(x = .data$K, y = .data$CV_error)) +
      ggplot2::geom_line() +
      ggplot2::geom_point(size = 2) +
      ggplot2::geom_vline(xintercept = best_k, linetype = "dashed",
                          color = "red") +
      ggplot2::annotate("text", x = best_k, y = max(cv),
                        label = paste0("Kneedle K=", best_k),
                        color = "red", vjust = -0.5, hjust = -0.1) +
      ggplot2::scale_x_continuous(breaks = K) +
      ggplot2::labs(x = "K", y = "CV error",
                    title = "ADMIXTURE CV error and Kneedle elbow") +
      ggplot2::theme_bw()

    distance_plot <- ggplot2::ggplot(kneedle_table,
                                     ggplot2::aes(x = .data$K, y = .data$distance)) +
      ggplot2::geom_line() +
      ggplot2::geom_point(size = 2) +
      ggplot2::geom_vline(xintercept = best_k, linetype = "dashed",
                          color = "red") +
      ggplot2::scale_x_continuous(breaks = K) +
      ggplot2::labs(x = "K", y = "Distance to first-last chord",
                    title = "Kneedle distance curve") +
      ggplot2::theme_bw()

    # dual-axis view: CV error (left) + Kneedle distance (right); the
    # distance is scaled so both curves share the left axis, and the right
    # axis shows the original distance scale via sec_axis()
    if (max(dist_to_line) > 0) {
      scale_factor <- max(cv) / max(dist_to_line)
    } else {
      scale_factor <- 1
    }

    combined_plot <- ggplot2::ggplot(kneedle_table,
                                     ggplot2::aes(x = .data$K)) +
      ggplot2::geom_line(ggplot2::aes(y = .data$CV_error),
                         color = "black") +
      ggplot2::geom_point(ggplot2::aes(y = .data$CV_error), size = 2) +
      ggplot2::geom_line(ggplot2::aes(y = .data$distance * scale_factor),
                         color = "blue", linetype = "dashed") +
      ggplot2::geom_point(ggplot2::aes(y = .data$distance * scale_factor),
                          color = "blue", size = 2) +
      ggplot2::geom_vline(xintercept = best_k, linetype = "dashed",
                          color = "red") +
      ggplot2::scale_x_continuous(breaks = K) +
      ggplot2::scale_y_continuous(
        name = "CV error",
        sec.axis = ggplot2::sec_axis(~ . / scale_factor,
                                     name = "Kneedle distance")
      ) +
      ggplot2::labs(x = "K",
                    title = "ADMIXTURE CV error and Kneedle distance") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.title.y.right = ggplot2::element_text(color = "blue"),
        axis.text.y.right  = ggplot2::element_text(color = "blue")
      )

    plots <- list(cv_plot = cv_plot,
                  distance_plot = distance_plot,
                  combined_plot = combined_plot)
  }

  structure(
    list(
      best_k           = best_k,
      max_curvature_k  = max_curvature_k,
      threshold_k      = threshold_k,
      kneedle_table    = kneedle_table,
      diff_table       = diff_table,
      curvature_table  = curvature_table,
      rel_change_table = rel_change_table,
      plots            = plots
    ),
    class = "admixture_k_results"
  )
}

#' Print ADMIXTURE K Selection Results
#'
#' @param x An object of class \code{admixture_k_results} from
#'   \code{\link{find_admixture_k}}.
#' @param ... Ignored.
#'
#' @return Invisibly returns \code{x}.
#' @export
print.admixture_k_results <- function(x, ...) {
  cat("ADMIXTURE CV-error K selection
")
  cat("  Best K (Kneedle elbow):      ", x$best_k, "
")
  cat("  Max curvature K (2nd diff):  ", x$max_curvature_k, "
")
  for (nm in names(x$threshold_k)) {
    t <- sub("^rel_pct_lt_", "", nm)
    cat("  First K with rel drop <", t, "%:", x$threshold_k[[nm]], "
")
  }
  invisible(x)
}
