test_that("pairwise differential analysis rejects multi-factor formulas", {
  set.seed(3)
  counts <- matrix(
    rpois(80 * 12, lambda = 40),
    nrow = 80,
    dimnames = list(paste0("G", 1:80), paste0("S", 1:12))
  )
  meta <- data.frame(
    group = rep(c("A", "B"), each = 6),
    batch = rep(c("b1", "b2"), 6),
    row.names = colnames(counts)
  )

  expect_error(
    find_degs_deseq2(counts, meta, formula = ~ group + batch, pairwise = TRUE),
    "single-factor formula"
  )
  expect_error(
    find_degs_edger(counts, meta, formula = ~ group + batch, pairwise = TRUE),
    "single-factor formula"
  )
})

test_that("DESeq2 wrapper preserves unavailable adjusted p-values", {
  counts <- bioRtools::df.rnaseq.gene
  counts <- counts[rowSums(counts) > 0, , drop = FALSE]
  counts <- counts[seq_len(min(3000, nrow(counts))), , drop = FALSE]
  meta <- bioRtools::df.rnaseq.sample

  res <- suppressMessages(suppressWarnings(
    find_degs_deseq2(
      counts,
      meta,
      formula = ~ group,
      shrink_lfc = FALSE,
      alpha = 0.05,
      pairwise = TRUE
    )
  ))

  expect_true(anyNA(res$padj))
  expect_true(all(res$regulation[is.na(res$padj)] == "NS"))
})

test_that("plot_volcano uses regulation labels when available", {
  df <- data.frame(
    gene = paste0("G", 1:9),
    baseMean = seq(10, 90, length.out = 9),
    log2FoldChange = c(-3, -2, -1, -0.2, 0, 0.2, 1, 2, 3),
    padj = c(1e-5, 1e-4, 0.01, 0.4, 1, 0.4, 0.01, 1e-4, 1e-5),
    regulation = c("Down", "Down", "Down", "NS", "NS", "NS", "Up", "Up", "Up")
  )

  out <- plot_volcano(df, xintercept_abs = 0.5)

  expect_s3_class(out$plot_volcano, "ggplot")
  expect_equal(
    out$data_summary$count[match(c("Down", "NS", "Up"), out$data_summary$regulation)],
    c(3L, 3L, 3L)
  )
  expect_equal(rlang::as_label(out$plot_volcano$mapping$colour), "regulation")
  expect_s3_class(out$plot_volcano$scales$get_scales("colour"), "ScaleDiscrete")
})
