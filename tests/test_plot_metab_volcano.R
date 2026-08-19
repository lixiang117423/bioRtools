# Test script for plot_metab_volcano()
# Metabolomics-oriented multi-group volcano: NS points are horizontally
# jittered, Up/Down points are drawn exactly and always on top of NS.

# Source the function directly (for testing before package installation)
source("R/academic_color_scales.R")
source("R/plot_multi_volcano.R")   # provides internal detect_volcano_cols()
source("R/plot_metab_volcano.R")

library(ggplot2)

# --- Synthetic opls_analysis differential_analysis-style data ---------------
set.seed(1)
n_feat <- 300L
mk <- function(cmp) {
  data.frame(
    comparison = cmp,
    feature_id = paste0("metab_", seq_len(n_feat)),
    log2FoldChange = rnorm(n_feat, 0, 0.8),
    padj = runif(n_feat, 0, 1),
    stringsAsFactors = FALSE
  )
}
df <- rbind(mk("T2 vs T1"), mk("T3 vs T1"))
df$padj[df$log2FoldChange > 1.5] <- 1e-4    # strong Up
df$padj[df$log2FoldChange < -1.5] <- 1e-3   # strong Down
df$padj[seq_len(30)] <- NA                  # untestable features

n_sig_expected <- sum(!is.na(df$padj) & df$padj < 0.05 &
                      abs(df$log2FoldChange) >= 1)

# --- Test 1: returns a ggplot object ----------------------------------------
cat("Test 1: returns ggplot\n")
p <- plot_metab_volcano(df)
stopifnot(inherits(p, "ggplot"))
cat("✓ Test 1 passed\n\n")

# --- Test 2: both point layers jittered, sig drawn on top of NS ------------
cat("Test 2: layer positions and order (Up/Down on top of NS)\n")
# ggplot2 >= 3.5 implements geom_jitter as GeomPoint + PositionJitter,
# so check the position: NS jittered, sig jittered too, sig after NS.
pos_cls <- vapply(p$layers, function(l) class(l$position)[1], character(1))
cat("  layer positions:", paste(pos_cls, collapse = ", "), "\n")
stopifnot(pos_cls[3] == "PositionJitter")    # NS layer is jittered
stopifnot(pos_cls[4] == "PositionJitter")    # sig layer is jittered
cat("✓ Test 2 passed\n\n")

# --- Test 3: both layers jitter horizontally but keep exact log2FC ----------
cat("Test 3: horizontal jitter only, y positions exact\n")
b <- ggplot2::ggplot_build(p)
ns_layer <- b$data[[3]]
sig_layer <- b$data[[4]]
stopifnot(nrow(ns_layer) + nrow(sig_layer) == nrow(df))
# jitter spreads x within each of the 2 groups
stopifnot(length(unique(round(ns_layer$x, 3))) > 2)
stopifnot(length(unique(round(sig_layer$x, 3))) > 2)
# height = 0: every y (NS and sig) comes straight from the input
stopifnot(all(ns_layer$y %in% df$log2FoldChange))
stopifnot(all(sig_layer$y %in% df$log2FoldChange))
cat("✓ Test 3 passed\n\n")

# --- Test 4: NA padj features are kept visible as NS ------------------------
cat("Test 4: NA-padj features plotted as NS\n")
n_na_feats <- sum(is.na(df$padj))
stopifnot(nrow(ns_layer) == nrow(df) - n_sig_expected)
stopifnot(n_na_feats == 30)
cat("✓ Test 4 passed\n\n")

# --- Test 5: significance coloring matches thresholds -----------------------
cat("Test 5: colored points match padj/log2FC thresholds\n")
stopifnot(nrow(sig_layer) == n_sig_expected)
stopifnot(nrow(sig_layer) > 0)
cat("✓ Test 5 passed (", nrow(sig_layer), " colored points)\n\n")

# --- Test 6: labels present, one per direction per group --------------------
cat("Test 6: labels\n")
lab_layer <- b$data[[5]]
stopifnot(nrow(lab_layer) > 0)
stopifnot(nrow(lab_layer) <= 2 * 5 * 2)   # top_n=5, 2 directions, 2 groups
cat("✓ Test 6 passed (", nrow(lab_layer), " labels)\n\n")

# --- Test 7: custom thresholds, sizes, colors -------------------------------
cat("Test 7: custom arguments\n")
p2 <- plot_metab_volcano(df, log2fc_thr = 0.5, padj_thr = 0.01,
                         jitter_width = 0.2, ns_color = "grey60",
                         ns_size = 1.2, ns_alpha = 0.8,
                         sig_size = 2.5, sig_alpha = 0.9,
                         up_color = "#E64B35FF", down_color = "#3C5488FF",
                         top_n = 3, label_size = 3)
b2 <- ggplot2::ggplot_build(p2)
n_sig2 <- sum(!is.na(df$padj) & df$padj < 0.01 &
              abs(df$log2FoldChange) >= 0.5)
stopifnot(nrow(b2$data[[4]]) == n_sig2)
stopifnot(nrow(b2$data[[5]]) <= 2 * 3 * 2)
cat("✓ Test 7 passed\n\n")

# --- Test 8: jitter_width = 0 degenerates to a line (opt-out) ---------------
cat("Test 8: jitter_width = 0\n")
p3 <- plot_metab_volcano(df, jitter_width = 0)
b3 <- ggplot2::ggplot_build(p3)
stopifnot(length(unique(round(b3$data[[3]]$x, 3))) == 2)
cat("✓ Test 8 passed\n\n")

# --- Test 9: error handling -------------------------------------------------
cat("Test 9: error handling\n")
tryCatch({
  plot_metab_volcano(data.frame(a = 1))
  stop("should have thrown an error for missing columns")
}, error = function(e) {
  cat("✓ Test 9a passed: correctly caught -", conditionMessage(e), "\n")
})
tryCatch({
  plot_metab_volcano(df, padj_thr = 2)
  stop("should have thrown an error for invalid padj_thr")
}, error = function(e) {
  cat("✓ Test 9b passed: correctly caught -", conditionMessage(e), "\n")
})
cat("\n")

# --- Test 10: save plot -----------------------------------------------------
cat("Test 10: save plot\n")
ggplot2::ggsave(tempfile(fileext = ".png"), p,
                width = 10, height = 6, dpi = 150)
cat("✓ Test 10 passed\n\n")

# --- Test 11: regulation column drives coloring -----------------------------
cat("Test 11: regulation column (VIP-based) drives coloring\n")
df_reg <- df
df_reg$regulation <- ifelse(seq_len(nrow(df_reg)) %% 3 == 0, "Up",
                     ifelse(seq_len(nrow(df_reg)) %% 5 == 0, "Down", "NS"))
n_up <- sum(df_reg$regulation == "Up")
n_down <- sum(df_reg$regulation == "Down")
p4 <- plot_metab_volcano(df_reg)
b4 <- ggplot2::ggplot_build(p4)
# colored points == Up + Down from the regulation column
stopifnot(nrow(b4$data[[4]]) == n_up + n_down)
# NS layer == everything else (NA padj rows included)
stopifnot(nrow(b4$data[[3]]) == nrow(df_reg) - n_up - n_down)
# thresholds must NOT affect coloring when regulation is present
p5 <- plot_metab_volcano(df_reg, log2fc_thr = 5, padj_thr = 1e-10)
b5 <- ggplot2::ggplot_build(p5)
stopifnot(nrow(b5$data[[4]]) == n_up + n_down)
# explicit regulation_col
p6 <- plot_metab_volcano(df_reg, regulation_col = "regulation")
stopifnot(nrow(ggplot2::ggplot_build(p6)$data[[4]]) == n_up + n_down)
# missing regulation_col errors
tryCatch({
  plot_metab_volcano(df_reg, regulation_col = "nope")
  stop("should have thrown an error for missing regulation_col")
}, error = function(e) {
  cat("✓ Test 11 error case caught -", conditionMessage(e), "\n")
})
cat("✓ Test 11 passed (", n_up, " Up +", n_down, " Down colored)\n\n")

# --- Test 12: sig_jitter_width = 0 restores exact group-center line ---------
cat("Test 12: sig_jitter_width = 0\n")
p7 <- plot_metab_volcano(df_reg, sig_jitter_width = 0)
b7 <- ggplot2::ggplot_build(p7)
stopifnot(length(unique(round(b7$data[[4]]$x, 3))) == 2)   # 2 group centers
stopifnot(nrow(b7$data[[4]]) == n_up + n_down)              # count unchanged
# NS layer still jittered in the same call
stopifnot(length(unique(round(b7$data[[3]]$x, 3))) > 2)
# invalid sig_jitter_width errors
tryCatch({
  plot_metab_volcano(df_reg, sig_jitter_width = -1)
  stop("should have thrown an error for negative sig_jitter_width")
}, error = function(e) {
  cat("✓ Test 12 error case caught -", conditionMessage(e), "\n")
})
cat("✓ Test 12 passed\n\n")

cat("=====================================\n")
cat("All tests completed successfully! ✓\n")
cat("=====================================\n")
