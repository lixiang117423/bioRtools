#' Identify Differentially Abundant Microorganisms Using DESeq2 (Pairwise)
#'
#' Performs DESeq2 differential abundance analysis on all pairwise group
#' combinations (or against a reference group). For each pair, subsets the
#' count data and runs DESeq2's negative binomial model, then combines results
#' into a single data frame with a \code{comparison} column.
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
#' @param shrink.lfc Logical indicating whether to apply log2 fold change
#'   shrinkage for more accurate estimates (default: TRUE).
#' @param ref_group Optional character string specifying a reference group.
#'   When provided, all comparisons are treatment vs reference instead of
#'   all pairwise combinations.
#'
#' @return A data frame with DESeq2 results for significant features (Enriched
#'   or Depleted), sorted by comparison then adjusted p-value. Columns include:
#'   \itemize{
#'     \item \code{feature_id}: Feature identifiers
#'     \item \code{baseMean}, \code{log2FoldChange}, \code{lfcSE}, \code{stat},
#'       \code{pvalue}, \code{padj}: Standard DESeq2 output
#'     \item \code{significance}: "Enriched" or "Depleted"
#'     \item \code{comparison}: Label in "treatment vs reference" format
#'     \item \code{group}: Treatment group name
#'     \item \code{ref_group}: Reference group name
#'     \item \code{fold_change}: 2^log2FoldChange
#'     \item \code{abs_log2fc}: Absolute log2 fold change
#'     \item \code{group_mean}, \code{group_sd}, \code{group_n}: Treatment group statistics
#'     \item \code{ref_mean}, \code{ref_sd}, \code{ref_n}: Reference group statistics
#'     \item \code{test_method}: "DESeq2-Wald"
#'   }
#'
#' @references
#' Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold change
#' and dispersion for RNA-seq data with DESeq2. \emph{Genome Biology}, 15, 550.
#'
#' @seealso \code{\link{find_dams_lefse}}, \code{\link{find_degs_deseq2}}
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
#' # All pairwise comparisons
#' res <- find_dams_deseq2(mat, meta, groupCol = "treatment")
#'
#' # Against reference group only
#' res_ref <- find_dams_deseq2(mat, meta, groupCol = "treatment", ref_group = "A")
#' }
find_dams_deseq2 <- function(data, sample, groupCol = "group",
                             log2FoldChange = 1, padj = 0.05,
                             min.count = 1, min.samples = 2,
                             shrink.lfc = TRUE,
                             ref_group = NULL) {
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

  # --- Build pairwise pairs --------------------------------------------------
  if (!is.null(ref_group)) {
    if (!ref_group %in% groups) {
      stop("'ref_group' must be one of: ", paste(groups, collapse = ", "))
    }
    others <- setdiff(groups, ref_group)
    pairs <- lapply(others, function(g) c(ref_group, g))
  } else {
    pairs <- utils::combn(groups, 2, simplify = FALSE)
  }

  # --- Pairwise DESeq2 ------------------------------------------------------
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
    fmla <- as.formula(paste0("~", groupCol))

    dds <- tryCatch(
      DESeq2::DESeqDataSetFromMatrix(
        countData = data_sub,
        colData   = sample_sub,
        design    = fmla
      ),
      error = function(e) NULL
    )
    if (is.null(dds)) return(NULL)

    dds <- tryCatch(DESeq2::DESeq(dds, quiet = TRUE), error = function(e) NULL)
    if (is.null(dds)) return(NULL)

    res <- tryCatch({
      if (shrink.lfc) {
        res_raw <- tryCatch(
          DESeq2::lfcShrink(dds,
            coef = DESeq2::resultsNames(dds)[length(DESeq2::resultsNames(dds))],
            type = "apeglm", quiet = TRUE),
          error = function(e) DESeq2::results(dds))
      } else {
        res_raw <- DESeq2::results(dds)
      }
      as.data.frame(res_raw)
    }, error = function(e) NULL)
    if (is.null(res)) return(NULL)

    # Process results
    res <- tibble::rownames_to_column(res, "feature_id")
    res$padj <- ifelse(is.na(res$padj), 1, res$padj)

    res$significance <- ifelse(
      res$log2FoldChange > log2FoldChange & res$padj < padj, "Enriched",
      ifelse(res$log2FoldChange < -log2FoldChange & res$padj < padj, "Depleted", "NS")
    )

    sig <- res$significance != "NS"
    if (!any(sig)) return(NULL)

    res_sig <- res[sig, , drop = FALSE]

    # Comparison label: treatment vs reference
    res_sig$comparison <- paste(pair[2], "vs", pair[1])
    res_sig$fold_change <- 2^res_sig$log2FoldChange
    res_sig$abs_log2fc  <- abs(res_sig$log2FoldChange)
    res_sig$group <- pair[2]
    res_sig$ref_group <- pair[1]

    # Group descriptive statistics
    idx_grp <- sample_sub[[groupCol]] == pair[2]
    idx_ref <- sample_sub[[groupCol]] == pair[1]
    n_grp <- sum(idx_grp)
    n_ref <- sum(idx_ref)

    norm_counts <- DESeq2::counts(dds, normalized = TRUE)
    grp_means <- rowMeans(norm_counts[, idx_grp, drop = FALSE], na.rm = TRUE)
    grp_sds <- apply(norm_counts[, idx_grp, drop = FALSE], 1, sd, na.rm = TRUE)
    ref_means <- rowMeans(norm_counts[, idx_ref, drop = FALSE], na.rm = TRUE)
    ref_sds <- apply(norm_counts[, idx_ref, drop = FALSE], 1, sd, na.rm = TRUE)

    res_sig$group_mean <- round(grp_means[res_sig$feature_id], 2)
    res_sig$group_sd <- round(grp_sds[res_sig$feature_id], 2)
    res_sig$group_n <- n_grp
    res_sig$ref_mean <- round(ref_means[res_sig$feature_id], 2)
    res_sig$ref_sd <- round(ref_sds[res_sig$feature_id], 2)
    res_sig$ref_n <- n_ref
    res_sig$test_method <- "DESeq2-Wald"

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
      group         = character(0),
      ref_group     = character(0),
      group_mean    = numeric(0),
      group_sd      = numeric(0),
      group_n       = integer(0),
      ref_mean      = numeric(0),
      ref_sd        = numeric(0),
      ref_n         = integer(0),
      test_method   = character(0),
      stringsAsFactors = FALSE
    )
    return(empty)
  }

  rownames(combined) <- NULL
  combined[order(combined$comparison, combined$padj, -combined$abs_log2fc), ]
}
