#' Identify Differentially Abundant Microorganisms Using DESeq2 (Pairwise)
#'
#' Performs DESeq2 differential abundance analysis on all pairwise group
#' combinations. For each pair, subsets the count data and runs DESeq2's
#' negative binomial model, then combines results into a single data frame
#' with a \code{comparison} column.
#'
#' @param data Integer count matrix (rows = features, cols = samples).
#'   If a data.frame with a non-numeric first column, that column is used as
#'   row names (feature IDs). Must contain raw counts (non-negative integers).
#' @param sample Data frame of sample metadata. Row names should match
#'   \code{colnames(data)}. If row names are missing, the function will try
#'   to find a column whose values match \code{colnames(data)} and use it as
#'   row names.
#' @param groupCol Column name in \code{sample} containing group labels
#'   (default: "group").
#' @param log2FoldChange Minimum absolute log2 fold change threshold
#'   (default: 1, i.e. 2-fold).
#' @param padj Adjusted p-value (FDR) threshold (default: 0.05).
#' @param min.count Minimum count per sample for feature filtering
#'   (default: 1).
#' @param min.samples Minimum number of samples that must have counts
#'   above \code{min.count} for a feature to be retained (default: 2).
#'
#' @return A data frame with DESeq2 results and a \code{comparison} column,
#'   sorted by comparison then adjusted p-value. Only features classified as
#'   "Enriched" or "Depleted" are returned.
#'
#' @references
#' Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold change
#' and dispersion for RNA-seq data with DESeq2. \emph{Genome Biology}, 15, 550.
#'
#' @seealso \code{\link{find_dams_lefse}}
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' set.seed(123)
#' mat <- matrix(rnbinom(100 * 20, size = 1, mu = 10), nrow = 100, ncol = 20,
#'   dimnames = list(paste0("OTU_", 1:100), paste0("S", 1:20)))
#' meta <- data.frame(row.names = colnames(mat),
#'   treatment = rep(c("A", "B", "C"), times = c(7, 7, 6)))
#'
#' \dontrun{
#' res <- find_dams_deseq2(mat, meta, groupCol = "treatment")
#' head(res)
#' }
find_dams_deseq2 <- function(data, sample, groupCol = "group",
                             log2FoldChange = 1, padj = 0.05,
                             min.count = 1, min.samples = 2) {
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

  # Auto-detect sample ID column
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

  if (!groupCol %in% names(sample)) {
    stop(sprintf("Column '%s' not found in sample data", groupCol))
  }

  groups <- unique(as.character(sample[[groupCol]]))
  if (length(groups) < 2) stop("At least 2 groups are required")

  # Ensure integer counts
  if (!is.integer(data[1, 1])) {
    data <- round(data)
    storage.mode(data) <- "integer"
  }

  # --- Pairwise DESeq2 ------------------------------------------------------
  pairs <- utils::combn(groups, 2, simplify = FALSE)

  results <- lapply(pairs, function(pair) {
    idx <- sample[[groupCol]] %in% pair
    data_sub   <- data[, idx, drop = FALSE]
    sample_sub <- sample[idx, , drop = FALSE]
    sample_sub[[groupCol]] <- factor(as.character(sample_sub[[groupCol]]),
                                     levels = pair)

    if (any(table(sample_sub[[groupCol]]) < 2)) return(NULL)

    # Filter low-count features within this pair
    keep <- rowSums(data_sub >= min.count) >= min.samples
    data_sub <- data_sub[keep, , drop = FALSE]
    if (nrow(data_sub) == 0) return(NULL)

    # Build formula
    formula_str <- paste0("~", groupCol)
    fmla <- as.formula(formula_str)

    dds <- suppressMessages(tryCatch(
      DESeq2::DESeqDataSetFromMatrix(
        countData = data_sub,
        colData   = sample_sub,
        design    = fmla
      ),
      error = function(e) NULL
    ))
    if (is.null(dds)) return(NULL)

    dds <- suppressMessages(tryCatch(
      DESeq2::DESeq(dds, quiet = TRUE),
      error = function(e) NULL
    ))
    if (is.null(dds)) return(NULL)

    res <- tryCatch(
      as.data.frame(DESeq2::results(dds)),
      error = function(e) NULL
    )
    if (is.null(res)) return(NULL)

    # Process
    res <- tibble::rownames_to_column(res, "feature_id")
    res$padj <- ifelse(is.na(res$padj), 1, res$padj)

    res$significance <- ifelse(
      res$log2FoldChange > log2FoldChange & res$padj < padj, "Enriched",
      ifelse(res$log2FoldChange < -log2FoldChange & res$padj < padj, "Depleted", "NS")
    )

    sig <- res$significance != "NS"
    if (!any(sig)) return(NULL)

    res_sig <- res[sig, , drop = FALSE]
    res_sig$comparison <- paste(pair[1], "vs", pair[2])
    res_sig$fold_change <- 2^res_sig$log2FoldChange
    res_sig$abs_log2fc  <- abs(res_sig$log2FoldChange)
    res_sig
  })

  # --- Combine ---------------------------------------------------------------
  combined <- do.call(rbind, results)

  if (is.null(combined)) {
    message("No significant differentially abundant features found")
    empty <- data.frame(
      feature_id    = character(0),
      baseMean      = numeric(0),
      log2FoldChange = numeric(0),
      lfcSE         = numeric(0),
      stat          = numeric(0),
      pvalue        = numeric(0),
      padj          = numeric(0),
      significance  = character(0),
      comparison    = character(0),
      fold_change   = numeric(0),
      abs_log2fc    = numeric(0),
      stringsAsFactors = FALSE
    )
    return(empty)
  }

  rownames(combined) <- NULL
  combined[order(combined$comparison, combined$padj, -combined$abs_log2fc), ]
}
