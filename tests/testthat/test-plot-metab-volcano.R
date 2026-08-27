# Migrated from tests/test_plot_metab_volcano.R (bare script).
# plot_metab_volcano(): metabolomics-oriented multi-group volcano. NS points
# are horizontally jittered, Up/Down points are drawn exactly and always on
# top of NS.

# Synthetic opls_analysis differential_analysis-style data
set.seed(1)
n_feat <- 300L
mk <- function(cmp) {
  data.frame(
    comparison = cmp,
    feature_id = paste0("metab_", seq_len(n_feat)),
    log2FoldChange = stats::rnorm(n_feat, 0, 0.8),
    padj = stats::runif(n_feat, 0, 1),
    stringsAsFactors = FALSE
  )
}
df <- rbind(mk("T2 vs T1"), mk("T3 vs T1"))
df$padj[df$log2FoldChange > 1.5] <- 1e-4    # strong Up
df$padj[df$log2FoldChange < -1.5] <- 1e-3   # strong Down
df$padj[seq_len(30)] <- NA                  # untestable features

n_sig_expected <- sum(!is.na(df$padj) & df$padj < 0.05 &
                      abs(df$log2FoldChange) >= 1)

# Regulation-column data: Up/Down/NS labels drive the coloring
df_reg <- df
df_reg$regulation <- ifelse(seq_len(nrow(df_reg)) %% 3 == 0, "Up",
                     ifelse(seq_len(nrow(df_reg)) %% 5 == 0, "Down", "NS"))
n_up <- sum(df_reg$regulation == "Up")
n_down <- sum(df_reg$regulation == "Down")

# Default plot and its built data (deterministic via the default jitter_seed)
p <- plot_metab_volcano(df)
b <- ggplot2::ggplot_build(p)
ns_layer <- b$data[[3]]
sig_layer <- b$data[[4]]
lab_layer <- b$data[[5]]

test_that("returns a ggplot object", {
  expect_true(ggplot2::is_ggplot(p))
})

test_that("both point layers jittered, sig drawn on top of NS", {
  # ggplot2 >= 3.5 implements geom_jitter as GeomPoint + PositionJitter,
  # so check the position: NS jittered, sig jittered too, sig after NS.
  pos_cls <- unname(vapply(p$layers, function(l) class(l$position)[1],
                          character(1)))
  expect_equal(pos_cls[3], "PositionJitter")    # NS layer is jittered
  expect_equal(pos_cls[4], "PositionJitter")    # sig layer is jittered
})

test_that("both layers jitter horizontally but keep exact log2FC", {
  expect_equal(nrow(ns_layer) + nrow(sig_layer), nrow(df))
  # jitter spreads x within each of the 2 groups
  expect_gt(length(unique(round(ns_layer$x, 3))), 2)
  expect_gt(length(unique(round(sig_layer$x, 3))), 2)
  # height = 0: every y (NS and sig) comes straight from the input
  expect_true(all(ns_layer$y %in% df$log2FoldChange))
  expect_true(all(sig_layer$y %in% df$log2FoldChange))
})

test_that("NA padj features are kept visible as NS", {
  n_na_feats <- sum(is.na(df$padj))
  expect_equal(nrow(ns_layer), nrow(df) - n_sig_expected)
  expect_equal(n_na_feats, 30)
})

test_that("significance coloring matches the thresholds", {
  expect_equal(nrow(sig_layer), n_sig_expected)
  expect_gt(nrow(sig_layer), 0)
})

test_that("labels present, one per direction per group", {
  expect_gt(nrow(lab_layer), 0)
  expect_lte(nrow(lab_layer), 2 * 5 * 2)   # top_n=5, 2 directions, 2 groups
})

test_that("custom thresholds, sizes, colors", {
  p2 <- plot_metab_volcano(df, log2fc_thr = 0.5, padj_thr = 0.01,
                           jitter_width = 0.2, ns_color = "grey60",
                           ns_size = 1.2, ns_alpha = 0.8,
                           sig_size = 2.5, sig_alpha = 0.9,
                           up_color = "#E64B35FF", down_color = "#3C5488FF",
                           top_n = 3, label_size = 3)
  b2 <- ggplot2::ggplot_build(p2)
  n_sig2 <- sum(!is.na(df$padj) & df$padj < 0.01 &
                abs(df$log2FoldChange) >= 0.5)
  expect_equal(nrow(b2$data[[4]]), n_sig2)
  expect_lte(nrow(b2$data[[5]]), 2 * 3 * 2)
})

test_that("jitter_width = 0 degenerates to a line (opt-out)", {
  p3 <- plot_metab_volcano(df, jitter_width = 0)
  b3 <- ggplot2::ggplot_build(p3)
  expect_equal(length(unique(round(b3$data[[3]]$x, 3))), 2)
})

test_that("error handling", {
  expect_error(plot_metab_volcano(data.frame(a = 1)),
               "Cannot auto-detect group column")
  expect_error(plot_metab_volcano(df, padj_thr = 2),
               "between 0 and 1")
})

test_that("plot can be saved to a file", {
  out_file <- tempfile(fileext = ".png")
  ggplot2::ggsave(out_file, p, width = 10, height = 6, dpi = 150)
  expect_true(file.exists(out_file))
})

test_that("regulation column (VIP-based) drives coloring", {
  p4 <- plot_metab_volcano(df_reg)
  b4 <- ggplot2::ggplot_build(p4)
  # colored points == Up + Down from the regulation column
  expect_equal(nrow(b4$data[[4]]), n_up + n_down)
  # NS layer == everything else (NA padj rows included)
  expect_equal(nrow(b4$data[[3]]), nrow(df_reg) - n_up - n_down)
  # thresholds must NOT affect coloring when regulation is present
  p5 <- plot_metab_volcano(df_reg, log2fc_thr = 5, padj_thr = 1e-10)
  b5 <- ggplot2::ggplot_build(p5)
  expect_equal(nrow(b5$data[[4]]), n_up + n_down)
  # explicit regulation_col
  p6 <- plot_metab_volcano(df_reg, regulation_col = "regulation")
  expect_equal(nrow(ggplot2::ggplot_build(p6)$data[[4]]), n_up + n_down)
  # missing regulation_col errors
  expect_error(plot_metab_volcano(df_reg, regulation_col = "nope"),
               "not found")
})

test_that("sig_jitter_width = 0 restores the exact group-center line", {
  p7 <- plot_metab_volcano(df_reg, sig_jitter_width = 0)
  b7 <- ggplot2::ggplot_build(p7)
  expect_equal(length(unique(round(b7$data[[4]]$x, 3))), 2)   # 2 group centers
  expect_equal(nrow(b7$data[[4]]), n_up + n_down)             # count unchanged
  # NS layer still jittered in the same call
  expect_gt(length(unique(round(b7$data[[3]]$x, 3))), 2)
  # invalid sig_jitter_width errors
  expect_error(plot_metab_volcano(df_reg, sig_jitter_width = -1),
               "non-negative")
})
