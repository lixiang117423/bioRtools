#' Identify Differentially Abundant Microorganisms Using LEfSe (Pairwise)
#'
#' Performs LEfSe-style analysis on all pairwise group combinations. For each
#' pair, the function first attempts \code{lefser::lefser()}. If that fails
#' (common with small sample sizes or high-dimensional data), it falls back to
#' a manual implementation: Wilcoxon rank-sum test + Fisher's linear
#' discriminant effect size. Results from all pairs are combined into a single
#' data frame with a \code{comparison} column.
#'
#' @param data Feature count matrix (rows = features, cols = samples).
#'   If a data.frame with a non-numeric first column, that column is used as
#'   row names (feature IDs). Accepts raw counts; the function converts to
#'   relative abundances internally.
#' @param sample Data frame of sample metadata. Row names should match
#'   \code{colnames(data)}. If row names are missing, the function will try
#'   to find a column whose values match \code{colnames(data)} and use it as
#'   row names.
#' @param group_col Column name in \code{sample} containing group labels
#'   (default: "group").
#' @param wilcox_threshold P-value threshold for Wilcoxon rank-sum test
#'   (default: 0.05).
#' @param lda_threshold Minimum absolute LDA score (default: 2.0).
#'
#' @return A data frame with columns:
#'   \itemize{
#'     \item \code{features}: Feature identifiers
#'     \item \code{scores}: LDA scores (positive = enriched in first group,
#'       negative = enriched in second group)
#'     \item \code{class}: Group where the feature is enriched
#'     \item \code{pvalue}: Wilcoxon test p-value (NA when using lefser)
#'     \item \code{comparison}: Pair label, e.g. "CK vs Treatment"
#'     \item \code{lda_score}: Absolute LDA score
#'     \item \code{method}: "lefser" or "wilcox_lda"
#'   }
#'
#' @examples
#' set.seed(42)
#' mat <- matrix(floor(runif(60 * 18, 0, 100)), nrow = 60, ncol = 18,
#'   dimnames = list(paste0("ASV_", 1:60), paste0("S", 1:18)))
#' meta <- data.frame(row.names = colnames(mat),
#'   treatment = rep(c("A", "B", "C"), each = 6))
#' res <- find_dams_lefse(mat, meta, group_col = "treatment")
#'
#' @references
#' Segata, N. et al. (2011). Metagenomic biomarker discovery and explanation.
#' \emph{Genome Biology}, 12(6), R60.
#'
#' @seealso \code{\link{find_dams_deseq2}}
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
find_dams_lefse <- function(data, sample, group_col = "group",
                            wilcox_threshold = 0.05,
                            lda_threshold = 2.0) {
  # --- Prepare data ----------------------------------------------------------
  data <- as.data.frame(data)
  if (ncol(data) > 1 && !is.numeric(data[[1]])) {
    feat_names <- data[[1]]
    data <- data[, -1, drop = FALSE]
  } else {
    feat_names <- rownames(data)
  }
  data <- as.matrix(data)
  rownames(data) <- feat_names

  sample <- as.data.frame(sample)

  # Auto-detect sample ID column when row names are missing/uninformative
  rn <- rownames(sample)
  if (is.null(rn) || identical(rn, as.character(seq_len(nrow(sample))))) {
    for (col in names(sample)) {
      if (is.character(sample[[col]]) && all(colnames(data) %in% sample[[col]])) {
        rownames(sample) <- sample[[col]]
        break
      }
    }
  }

  # Align samples
  common <- intersect(colnames(data), rownames(sample))
  if (length(common) == 0) {
    stop("No matching samples between data columns and sample row names")
  }
  data   <- data[, common, drop = FALSE]
  sample <- sample[common, , drop = FALSE]

  if (!group_col %in% names(sample)) {
    stop(sprintf("Column '%s' not found in sample data", group_col))
  }

  groups <- unique(as.character(sample[[group_col]]))
  if (length(groups) < 2) stop("At least 2 groups are required")

  # Convert to relative abundances
  cs <- colSums(data)
  cs[cs == 0] <- 1
  ra <- sweep(data, 2, cs, "/")

  # --- Pairwise analysis -----------------------------------------------------
  pairs <- utils::combn(groups, 2, simplify = FALSE)

  results <- lapply(pairs, function(pair) {
    g1_idx <- sample[[group_col]] == pair[1]
    g2_idx <- sample[[group_col]] == pair[2]
    n1 <- sum(g1_idx)
    n2 <- sum(g2_idx)
    if (n1 < 2 || n2 < 2) return(NULL)

    # Subset
    ra_sub   <- ra[, g1_idx | g2_idx, drop = FALSE]
    samp_sub <- sample[g1_idx | g2_idx, , drop = FALSE]
    samp_sub[[group_col]] <- factor(as.character(samp_sub[[group_col]]))

    # Try lefser first
    res_lefser <- try_lefser(ra_sub, samp_sub, group_col,
                             wilcox_threshold, lda_threshold)
    if (!is.null(res_lefser)) {
      res_lefser$comparison <- paste(pair[1], "vs", pair[2])
      res_lefser$method <- "lefser"
      return(res_lefser)
    }

    # Fallback: Wilcoxon + Fisher's LDA
    res_fb <- fallback_wilcox_lda(ra, g1_idx, g2_idx, pair,
                                  wilcox_threshold, lda_threshold)
    if (!is.null(res_fb)) {
      res_fb$method <- "wilcox_lda"
      return(res_fb)
    }
    NULL
  })

  # --- Combine ---------------------------------------------------------------
  combined <- do.call(rbind, results)

  if (is.null(combined)) {
    message("No significant biomarkers found across any pairwise comparison")
    return(data.frame(
      features   = character(0),
      scores     = numeric(0),
      class      = character(0),
      pvalue     = numeric(0),
      comparison = character(0),
      lda_score  = numeric(0),
      method     = character(0),
      stringsAsFactors = FALSE
    ))
  }

  rownames(combined) <- NULL
  combined[order(combined$comparison, -combined$lda_score), ]
}


# --- Internal: try lefser ----------------------------------------------------
try_lefser <- function(ra, sample_sub, group_col,
                       kruskal_threshold, lda_threshold) {
  if (!requireNamespace("lefser", quietly = TRUE) ||
      !requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    return(NULL)
  }

  # Pre-filter: keep features present in >= 2 samples per group
  grp <- sample_sub[[group_col]]
  levels_grp <- levels(grp)
  keep <- rowSums(ra[, grp == levels_grp[1], drop = FALSE] > 0) >= 2 &
          rowSums(ra[, grp == levels_grp[2], drop = FALSE] > 0) >= 2
  ra_sub <- ra[keep, , drop = FALSE]
  if (nrow(ra_sub) == 0) return(NULL)

  # Remove features constant within any group
  const <- sapply(levels_grp, function(g) {
    apply(ra_sub[, grp == g, drop = FALSE], 1, function(x) var(x) == 0)
  })
  if (is.matrix(const)) {
    keep2 <- !apply(const, 1, any)
  } else {
    keep2 <- !const
  }
  ra_sub <- ra_sub[keep2, , drop = FALSE]
  if (nrow(ra_sub) == 0) return(NULL)

  se <- tryCatch(
    SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = ra_sub),
      colData = sample_sub
    ),
    error = function(e) NULL
  )
  if (is.null(se)) return(NULL)

  res <- suppressMessages(suppressWarnings(tryCatch(
    lefser::lefser(
      se,
      classCol = group_col,
      kruskal.threshold = kruskal_threshold,
      wilcox_threshold  = kruskal_threshold,
      lda_threshold     = lda_threshold
    ),
    error = function(e) NULL
  )))

  if (is.null(res) || nrow(res) == 0) return(NULL)

  res$lda_score <- abs(res$scores)
  res$pvalue <- NA_real_
  res[, c("features", "scores", "class", "pvalue", "lda_score")]
}


# --- Internal: Wilcoxon + Fisher LDA fallback --------------------------------
fallback_wilcox_lda <- function(ra, g1_idx, g2_idx, pair,
                                wilcox_threshold, lda_threshold) {
  n1 <- sum(g1_idx)
  n2 <- sum(g2_idx)

  x1 <- ra[, g1_idx, drop = FALSE]
  x2 <- ra[, g2_idx, drop = FALSE]

  m1 <- rowMeans(x1)
  m2 <- rowMeans(x2)

  v1 <- apply(x1, 1, var)
  v2 <- apply(x2, 1, var)
  pool_sd <- sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2))
  pool_sd[pool_sd == 0] <- NA

  pvals <- apply(ra, 1, function(row) {
    tryCatch(
      stats::wilcox.test(row[g1_idx], row[g2_idx], exact = FALSE)$p.value,
      error = function(e) NA
    )
  })

  lda_raw <- (m1 - m2) / pool_sd * sqrt(n1 * n2 / (n1 + n2))

  sig <- !is.na(pvals) & pvals < wilcox_threshold &
         !is.na(lda_raw) & abs(lda_raw) >= lda_threshold

  if (!any(sig)) return(NULL)

  data.frame(
    features   = rownames(ra)[sig],
    scores     = lda_raw[sig],
    class      = ifelse(m1[sig] > m2[sig], pair[1], pair[2]),
    pvalue     = pvals[sig],
    comparison = paste(pair[1], "vs", pair[2]),
    lda_score  = abs(lda_raw[sig]),
    stringsAsFactors = FALSE
  )
}
