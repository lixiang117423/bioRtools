#' Calculate Alpha Diversity Indices
#'
#' @description
#' Calculate multiple alpha diversity indices from an ASV/OTU abundance table,
#' including Shannon, Simpson, Chao1, ACE, and Observed richness. Returns a
#' long-format data frame with all indices merged with sample metadata.
#'
#' @param data ASV/OTU abundance table. Can be in either orientation:
#'   features × samples (rows = ASVs, first column is feature ID) or
#'   samples × features (rows = samples, columns = features). Auto-detected.
#' @param sample Sample metadata table. Must contain sample names matching
#'   those in data. Default is NULL (no metadata joined).
#' @param index Character vector specifying which indices to calculate.
#'   Options: "shannon", "simpson", "observed", "chao1", "ace".
#'   Default is all: c("shannon", "simpson", "observed", "chao1", "ace").
#' @param verbose Logical indicating whether to print progress information.
#'   Default is TRUE.
#'
#' @return A data frame in long format with columns:
#'   \code{sample}, \code{alpha} (index name), \code{value}, plus all columns
#'   from the sample metadata table (if provided).
#'
#' @details
#' Diversity indices:
#' \itemize{
#'   \item Shannon: \code{vegan::diversity(index = "shannon")}
#'   \item Simpson: \code{vegan::diversity(index = "simpson")}
#'   \item Observed: Number of observed taxa (row sums > 0)
#'   \item Chao1: \code{vegan::estimateR()} — species richness estimator
#'   \item ACE: \code{vegan::estimateR()} — abundance-based coverage estimator
#' }
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' set.seed(42)
#' mat <- matrix(rpois(500, lambda = 5), nrow = 50,
#'   dimnames = list(paste0("ASV", 1:50), paste0("Sample", 1:10)))
#'
#' result <- calc_alpha_diversity(data = mat)
#' head(result)
#'
calc_alpha_diversity <- function(data,
                                 sample = NULL,
                                 index = c("shannon", "simpson",
                                           "observed", "chao1", "ace"),
                                 verbose = TRUE) {

  # ── Input validation ──────────────────────────────────────────────────────
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("'data' must be a data frame or matrix")
  }
  if (!is.null(sample) && !is.data.frame(sample)) {
    stop("'sample' must be a data frame or NULL")
  }

  valid_indices <- c("shannon", "simpson", "observed", "chao1", "ace")
  index <- match.arg(index, valid_indices, several.ok = TRUE)

  # ── Auto-detect orientation and prepare matrix ────────────────────────────
  data <- as.data.frame(data)

  # Determine sample column for matching
  sample_names <- if (!is.null(sample)) as.character(sample[[1]]) else NULL

  if (!is.null(sample_names)) {
    data_colnames <- colnames(data)
    data_rownames <- rownames(data)

    cols_match <- if (ncol(data) > 1)
      sum(data_colnames[-1] %in% sample_names) else 0
    rows_match <- if (!is.null(data_rownames))
      sum(data_rownames %in% sample_names) else 0

    if (cols_match > rows_match && cols_match >= 2) {
      # Features × samples: first column is feature IDs
      if (verbose) message("Detected features x samples format, transposing...")
      feature_names <- as.character(data[[1]])
      data_matrix <- as.matrix(data[, -1, drop = FALSE])
      storage.mode(data_matrix) <- "numeric"
      data_matrix <- t(data_matrix)
      colnames(data_matrix) <- feature_names
    } else if (rows_match >= 2) {
      # Samples × features: rownames are sample IDs
      data_matrix <- as.matrix(data)
      storage.mode(data_matrix) <- "numeric"
    } else {
      stop("Cannot match data to sample names")
    }
  } else {
    # No sample table: assume rows = samples, columns = features
    if (is.character(data[[1]])) {
      # First column might be feature IDs → features × samples
      data_matrix <- as.matrix(data[, -1, drop = FALSE])
      storage.mode(data_matrix) <- "numeric"
      data_matrix <- t(data_matrix)
    } else {
      data_matrix <- as.matrix(data)
      storage.mode(data_matrix) <- "numeric"
    }
  }

  if (verbose) {
    message("Matrix dimensions: ", nrow(data_matrix), " samples x ",
            ncol(data_matrix), " features")
    message("Calculating: ", paste(index, collapse = ", "))
  }

  # ── Calculate indices ─────────────────────────────────────────────────────
  results <- list()
  samples_vec <- rownames(data_matrix)

  # Shannon
  if ("shannon" %in% index) {
    shannon_val <- vegan::diversity(data_matrix, index = "shannon")
    results$shannon <- data.frame(
      sample = samples_vec,
      alpha  = "shannon",
      value  = shannon_val,
      stringsAsFactors = FALSE
    )
  }

  # Simpson
  if ("simpson" %in% index) {
    simpson_val <- vegan::diversity(data_matrix, index = "simpson")
    results$simpson <- data.frame(
      sample = samples_vec,
      alpha  = "simpson",
      value  = simpson_val,
      stringsAsFactors = FALSE
    )
  }

  # Observed richness
  if ("observed" %in% index) {
    observed_val <- rowSums(data_matrix > 0)
    results$observed <- data.frame(
      sample = samples_vec,
      alpha  = "observed",
      value  = observed_val,
      stringsAsFactors = FALSE
    )
  }

  # Chao1 and ACE (from estimateR)
  needs_estimateR <- c("chao1", "ace")
  if (any(needs_estimateR %in% index)) {
    est <- vegan::estimateR(data_matrix)
    est_df <- as.data.frame(est)
    # Row names like "S.obs.ASV1", "S.chao1.ASV1" etc.
    # Actually estimateR returns: S.obs, S.chao1, se.chao1, S.ACE, se.ACE
    # Columns are samples

    if ("chao1" %in% index) {
      chao1_val <- est_df["S.chao1", ]
      results$chao1 <- data.frame(
        sample = samples_vec,
        alpha  = "chao1",
        value  = as.numeric(chao1_val),
        stringsAsFactors = FALSE
      )
    }

    if ("ace" %in% index) {
      ace_val <- est_df["S.ACE", ]
      results$ace <- data.frame(
        sample = samples_vec,
        alpha  = "ace",
        value  = as.numeric(ace_val),
        stringsAsFactors = FALSE
      )
    }
  }

  # ── Combine results ───────────────────────────────────────────────────────
  df_result <- do.call(dplyr::bind_rows, results)

  # ── Join sample metadata ──────────────────────────────────────────────────
  if (!is.null(sample)) {
    df_result <- df_result %>%
      dplyr::left_join(sample, by = "sample")
  }

  if (verbose) {
    message("Done. ", length(unique(df_result$alpha)), " indices for ",
            length(unique(df_result$sample)), " samples")
  }

  df_result
}
