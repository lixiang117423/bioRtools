# Internal helpers for the sample-metadata contract shared by the core
# workflows. They deliberately fail on ambiguous or incomplete IDs instead of
# silently falling back to row position.

has_default_rownames <- function(x) {
  rn <- rownames(x)
  is.null(rn) || identical(rn, as.character(seq_len(nrow(x))))
}

resolve_sample_ids <- function(sample, sample_col = NULL) {
  if (!is.data.frame(sample)) {
    stop("'sample' must be a data frame", call. = FALSE)
  }
  if (nrow(sample) == 0L) {
    stop("'sample' must contain at least one row", call. = FALSE)
  }

  if (!is.null(sample_col)) {
    if (!is.character(sample_col) || length(sample_col) != 1L) {
      stop("'sample_col' must be a single column name", call. = FALSE)
    }
    if (identical(sample_col, "row.names")) {
      ids <- rownames(sample)
    } else if (sample_col %in% names(sample)) {
      ids <- as.character(sample[[sample_col]])
    } else {
      stop("Sample ID column not found in metadata: ", sample_col, call. = FALSE)
    }
  } else {
    preferred <- intersect(
      c("sample", "sample_id", "SampleID", "sample_name", "id", "ID"),
      names(sample)
    )
    if (length(preferred) > 0L) {
      ids <- as.character(sample[[preferred[1]]])
    } else if (!has_default_rownames(sample)) {
      ids <- rownames(sample)
    } else {
      candidates <- names(sample)[vapply(sample, function(x) {
        is.character(x) || is.factor(x)
      }, logical(1))]
      candidates <- candidates[vapply(sample[candidates], function(x) {
        !anyNA(x) && !anyDuplicated(as.character(x))
      }, logical(1))]
      if (length(candidates) == 0L) {
        stop(
          "Cannot identify sample IDs in 'sample'; add a sample ID column or non-default row names",
          call. = FALSE
        )
      }
      ids <- as.character(sample[[candidates[1]]])
    }
  }

  if (length(ids) != nrow(sample) || anyNA(ids) || any(!nzchar(ids))) {
    stop("Sample IDs in metadata must be complete, non-empty, and one per row", call. = FALSE)
  }
  if (anyDuplicated(ids)) {
    dup <- unique(ids[duplicated(ids)])
    stop(
      "Duplicate sample IDs in metadata: ",
      paste(utils::head(dup, 5), collapse = ", "),
      call. = FALSE
    )
  }
  unname(ids)
}

align_sample_metadata <- function(sample, target_ids, sample_col = NULL) {
  target_ids <- as.character(target_ids)
  if (length(target_ids) != nrow(sample)) {
    stop(
      "Number of response rows (", length(target_ids),
      ") must equal number of metadata rows (", nrow(sample), ")",
      call. = FALSE
    )
  }
  if (anyNA(target_ids) || any(!nzchar(target_ids))) {
    stop("Response data must have complete, non-empty sample IDs in row names", call. = FALSE)
  }
  if (anyDuplicated(target_ids)) {
    stop("Response data contains duplicate sample IDs", call. = FALSE)
  }

  metadata_ids <- resolve_sample_ids(sample, sample_col = sample_col)
  missing_ids <- setdiff(target_ids, metadata_ids)
  extra_ids <- setdiff(metadata_ids, target_ids)
  if (length(missing_ids) > 0L || length(extra_ids) > 0L) {
    details <- c(
      if (length(missing_ids)) paste0("missing: ", paste(utils::head(missing_ids, 5), collapse = ", ")),
      if (length(extra_ids)) paste0("extra: ", paste(utils::head(extra_ids, 5), collapse = ", "))
    )
    stop(
      "Sample IDs in response and metadata do not match (",
      paste(details, collapse = "; "),
      ")",
      call. = FALSE
    )
  }

  idx <- match(target_ids, metadata_ids)
  list(
    sample = sample[idx, , drop = FALSE],
    sample_ids = target_ids,
    original_ids = metadata_ids
  )
}
