#' Estimate LD Decay and Automatically Recommend an Appropriate r² Threshold
#'
#' Fits a Hill & Weir (1988) nonlinear model to linkage disequilibrium (LD)
#' decay data and recommends a suitable r² threshold based on the ratio of
#' candidate thresholds to the background LD level estimated from distant
#' marker pairs.
#'
#' The Hill & Weir expectation formula used for fitting is:
#' \deqn{
#'   E(r^2) = \frac{10 + \rho d}{22 + 13\rho d + \rho^2 d^2}
#'            \left(1 + \frac{3 + \rho d}{n(22 + 13\rho d + \rho^2 d^2)}\right)
#' }
#' where \eqn{\rho} is the population recombination rate parameter and
#' \eqn{d} is physical distance in kb.
#'
#' Threshold recommendation logic:
#' \itemize{
#'   \item Thresholds below the background LD level are flagged as unsuitable.
#'   \item Thresholds within 1.5x of background LD are flagged as borderline.
#'   \item Thresholds between 1.5x and 4x background LD are recommended.
#'   \item Thresholds above 4x background LD are flagged as overly stringent.
#'   \item Among recommended thresholds, the one whose ratio to background LD
#'         is closest to 2 is selected as the best single recommendation.
#' }
#'
#' @param data A data frame containing at least two columns: physical distance
#'   and an LD statistic (e.g., r2).
#' @param n Integer. Number of individuals used when computing LD. Passed
#'   directly into the Hill & Weir expectation formula.
#' @param dist_col Character. Name of the column in \code{data} that contains
#'   physical distances. Default: \code{"dist"}.
#' @param value_col Character. Name of the column in \code{data} that contains
#'   the LD statistic. Default: \code{"r2"}.
#' @param dist_unit Character. Unit of the distance column; either \code{"bp"}
#'   (converted internally to kb) or \code{"kb"} (used as-is).
#'   Default: \code{"bp"}.
#' @param max_dist Numeric or \code{NULL}. If provided, only marker pairs
#'   within this distance (in kb, after unit conversion) are used for fitting.
#'   Default: \code{NULL} (no filtering).
#' @param bg_quantile Numeric in (0, 1). Marker pairs whose distance exceeds
#'   this quantile of all pairwise distances are used to estimate background
#'   LD. Default: \code{0.9}.
#' @param start_rho Numeric vector. Candidate starting values for the rho
#'   parameter tried sequentially until \code{\link[stats]{nls}} converges.
#'   Default: \code{c(0.1, 0.01, 1, 0.001, 10)}.
#' @param candidate_thresholds Numeric vector. r2 values evaluated for
#'   threshold recommendation.
#'   Default: \code{c(0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.02)}.
#' @param verbose Logical. If \code{TRUE}, progress messages and a summary
#'   table are printed to the console. Default: \code{TRUE}.
#'
#' @return A named list (returned invisibly) with the following elements:
#' \describe{
#'   \item{\code{fit}}{The \code{nls} model object.}
#'   \item{\code{rho}}{Estimated rho (recombination rate parameter).}
#'   \item{\code{background_r2}}{Mean r2 among distant marker pairs
#'     (distance >= \code{bg_quantile} quantile).}
#'   \item{\code{decay_table}}{A data frame with columns \code{threshold},
#'     \code{decay_kb}, \code{bg_ratio}, and \code{recommendation} for every
#'     value in \code{candidate_thresholds} that lies within the fitted curve
#'     range.}
#'   \item{\code{recommended}}{A single-row data frame from \code{decay_table}
#'     identified as the best threshold.}
#'   \item{\code{gwas_window_kb}}{Suggested GWAS flanking window in kb
#'     (recommended decay distance rounded up to the nearest 50 kb).}
#' }
#'
#' @references
#' Hill, W. G., & Weir, B. S. (1988). Variances and covariances of squared
#' linkage disequilibria in finite populations. \emph{Theoretical Population
#' Biology}, 33(1), 54-78. \doi{10.1016/0040-5809(88)90004-4}
#'
#' @importFrom dplyr select filter mutate summarise pull arrange slice case_when
#' @importFrom rlang sym
#' @importFrom stats nls coef predict quantile
#' @importFrom purrr map_dfr
#'
#' @examples
#' \dontrun{
#' library(readr)
#' library(dplyr)
#'
#' ld <- read_delim("ld_result.txt") |>
#'     select(dist, r2)
#'
#' # Basic usage — distance in bp (default)
#' result <- ld_decay_threshold(ld, n = 805)
#'
#' # Custom column names, distance already in kb
#' result <- ld_decay_threshold(
#'     data      = my_df,
#'     n         = 200,
#'     dist_col  = "distance_bp",
#'     value_col = "R2",
#'     dist_unit = "bp"
#' )
#'
#' # Access results
#' result$recommended      # best threshold row
#' result$gwas_window_kb   # suggested GWAS window
#' result$decay_table      # full threshold table
#' result$background_r2    # estimated background LD
#' }
#'
#' @export
ld_decay_threshold <- function(
        data,
        n,
        dist_col             = "dist",
        value_col            = "r2",
        dist_unit            = c("bp", "kb"),
        max_dist             = NULL,
        bg_quantile          = 0.9,
        start_rho            = c(0.1, 0.01, 1, 0.001, 10),
        candidate_thresholds = c(0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.02),
        verbose              = TRUE
) {

    # ── 0. Input validation ────────────────────────────────────────────────────
    if (!is.data.frame(data)) {
        stop("`data` must be a data frame.", call. = FALSE)
    }
    if (!dist_col %in% names(data)) {
        stop(sprintf("Column '%s' not found in `data`.", dist_col), call. = FALSE)
    }
    if (!value_col %in% names(data)) {
        stop(sprintf("Column '%s' not found in `data`.", value_col), call. = FALSE)
    }
    if (!is.numeric(n) || length(n) != 1L || n <= 0) {
        stop("`n` must be a single positive number.", call. = FALSE)
    }
    if (!is.numeric(bg_quantile) || bg_quantile <= 0 || bg_quantile >= 1) {
        stop("`bg_quantile` must be a numeric value in (0, 1).", call. = FALSE)
    }

    dist_unit <- match.arg(dist_unit)

    # ── 1. Data preparation ────────────────────────────────────────────────────
    df <- data |>
        dplyr::select(
            dist = !!rlang::sym(dist_col),
            r2   = !!rlang::sym(value_col)
        ) |>
        dplyr::filter(!is.na(dist), !is.na(r2), r2 > 0, dist >= 0) |>
        dplyr::mutate(
            dist_kb = if (dist_unit == "bp") dist / 1000 else dist
        )

    if (!is.null(max_dist)) {
        df <- dplyr::filter(df, dist_kb <= max_dist)
    }

    if (nrow(df) < 10L) {
        stop(
            "Fewer than 10 valid observations remain after filtering. ",
            "Check your data or relax `max_dist`.",
            call. = FALSE
        )
    }

    if (verbose) {
        message(sprintf(
            "[1/4] Data prepared: %d observations used for fitting.",
            nrow(df)
        ))
    }

    # ── 2. Hill & Weir nonlinear least-squares fitting ────────────────────────
    hw_formula <- r2 ~ ((10 + rho * dist_kb) /
                            (22 + 13 * rho * dist_kb + rho^2 * dist_kb^2)) *
        (1 + ((3 + rho * dist_kb) /
                  (n * (22 + 13 * rho * dist_kb + rho^2 * dist_kb^2))))

    fit      <- NULL
    rho_used <- NA_real_

    for (rho0 in start_rho) {
        fit <- tryCatch(
            stats::nls(
                formula = hw_formula,
                data    = df,
                start   = list(rho = rho0),
                control = list(maxiter = 1000, tol = 1e-6)
            ),
            error   = function(e) NULL,
            warning = function(w) NULL
        )
        if (!is.null(fit)) {
            rho_used <- rho0
            break
        }
    }

    if (is.null(fit)) {
        stop(
            "Hill & Weir model failed to converge for all starting values: ",
            paste(start_rho, collapse = ", "), ".\n",
            "Consider providing different `start_rho` values or checking data quality.",
            call. = FALSE
        )
    }

    rho_est <- stats::coef(fit)[["rho"]]

    if (verbose) {
        message(sprintf(
            "[2/4] Hill & Weir model converged: rho = %.6f (start = %.4f).",
            rho_est, rho_used
        ))
    }

    # ── 3. Background LD estimation ───────────────────────────────────────────
    dist_bg_cutoff <- stats::quantile(df$dist_kb, probs = bg_quantile)

    background_r2 <- df |>
        dplyr::filter(dist_kb >= dist_bg_cutoff) |>
        dplyr::summarise(bg = mean(r2)) |>
        dplyr::pull(bg)

    if (verbose) {
        message(sprintf(
            "[3/4] Background LD: mean r2 = %.4f (distances >= %.1f kb; top %.0f%% quantile).",
            background_r2, dist_bg_cutoff, bg_quantile * 100
        ))
    }

    # ── 4. Decay distance for each candidate threshold ────────────────────────
    dist_seq  <- seq(0, max(df$dist_kb), length.out = 10000L)
    fit_curve <- data.frame(
        dist_kb = dist_seq,
        r2_fit  = stats::predict(fit, newdata = data.frame(dist_kb = dist_seq))
    )

    decay_table <- purrr::map_dfr(
        sort(candidate_thresholds, decreasing = TRUE),
        function(thr) {
            hit <- dplyr::filter(fit_curve, r2_fit <= thr)
            if (nrow(hit) == 0L) return(NULL)
            data.frame(threshold = thr, decay_kb = round(hit$dist_kb[1L], 1))
        }
    )

    decay_table <- decay_table |>
        dplyr::mutate(
            bg_ratio       = round(threshold / background_r2, 2),
            recommendation = dplyr::case_when(
                threshold < background_r2      ~ "not recommended: below background LD",
                bg_ratio  < 1.5               ~ "caution: close to background LD",
                bg_ratio >= 1.5 & bg_ratio < 4 ~ "recommended",
                bg_ratio >= 4                  ~ "caution: far above background LD",
                TRUE                           ~ NA_character_
            )
        )

    # ── 5. Select the single best threshold ───────────────────────────────────
    recommended_rows <- dplyr::filter(decay_table, recommendation == "recommended")

    if (nrow(recommended_rows) == 0L) {
        warning(
            "No threshold met the recommendation criteria. ",
            "Returning the threshold with bg_ratio closest to 2.",
            call. = FALSE
        )
        recommended_rows <- decay_table
    }

    # Pick the row whose bg_ratio is closest to 2 (signal-to-background sweet spot)
    best_idx       <- which.min(abs(recommended_rows$bg_ratio - 2))
    best_row       <- recommended_rows[best_idx, , drop = FALSE]
    gwas_window_kb <- ceiling(best_row$decay_kb / 50) * 50

    if (verbose) {
        message("[4/4] Threshold recommendation summary:\n")
        print(decay_table, row.names = FALSE)
        message(sprintf(
            "\n  Best threshold : r2 = %.2f  |  decay distance = %.1f kb  |  bg_ratio = %.2f",
            best_row$threshold, best_row$decay_kb, best_row$bg_ratio
        ))
        message(sprintf(
            "  Suggested GWAS window : lead SNP +/- %d kb",
            gwas_window_kb
        ))
    }

    # ── 6. Return ──────────────────────────────────────────────────────────────
    invisible(list(
        fit            = fit,
        rho            = rho_est,
        background_r2  = background_r2,
        decay_table    = decay_table,
        recommended    = best_row,
        gwas_window_kb = gwas_window_kb
    ))
}
