# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.58.1] - 2026-06-30

### Fixed
Crash in the interactive console summary of `enrich_kegg()` / `enrich_go()` when results were non-empty (e.g. `Error in if (nchar(pathway_name) > 50) ... : 参数长度为零`).

- The summary block referenced the **pre-rename** output column names (`Description`, `ID`, `p.adjust`, `gene.count`, `total.genes`, `enrichment.score`) that no longer exist after the columns are renamed to snake_case earlier in the function — so `$Description` was `NULL`, `nchar(NULL)` was length 0, and `if()` crashed.
- Now uses the renamed columns (`description`, `id`, `p_adjust`, `gene_count`, `total_genes`, `enrichment_score`); added an `!is.na(...)` guard on the name-length check.
- Regression test: `tests/test_enrich_summary.R` (exercises the summary loop on real enriched results; verified end-to-end with the interactive path forced on).

## [1.58.0] - 2026-06-30

### Changed (breaking)
`enrich_kegg()` / `enrich_go()` now use **snake_case** for the required annotation-database columns, matching CLAUDE.md §1 and the rest of the package.

- `enrich_kegg()`: `kegg.id` / `kegg.term` / `kegg.category` → `kegg_id` / `kegg_term` / `kegg_category`.
- `enrich_go()`: `go.id` / `go.term` / `go.ontology` → `go_id` / `go_term` / `go_ontology`.
- Backward compatible via a small shim: legacy dot-case column names are still accepted but trigger a one-time deprecation warning (`reconcile_db_columns()`, `@keywords internal`). Existing databases keep working; rename them when convenient.
- Updated the shipped example databases `df.rnaseq.kegg` / `df.rnaseq.go` to snake_case (data unchanged), `data.R` docs, and the `utils::globalVariables()` lists.
- Tests: new `tests/test_enrich_columns.R` (snake_case runs; legacy dot.case reproduces identical results via the shim).

## [1.57.0] - 2026-06-30

### Changed
The ordination / dimension-reduction family now accepts a **features × samples** abundance matrix directly (genes as rows, samples as columns) — the same orientation `find_degs_deseq2` / `find_degs_edger` produce — so an RNA-seq workflow no longer needs a manual `t() %>% as.data.frame()` between DEG and ordination steps.

- Added a `feature_as_row` argument (default `NA` = auto-detect) to `pca_analysis`, `pcoa_analysis`, `opls_analysis`, `rda_analysis`, `tsne_analysis`, and `pairwise_oplsda`. With `sample` metadata supplied, the function matches sample IDs against both the row and column names of `data`; a column match triggers an internal transpose. `TRUE` / `FALSE` force an orientation.
- Extracted the logic into shared internal helpers `detect_feature_as_row()` / `orient_to_sample_row()` (new `R/utils-orientation.R`, `@keywords internal`); refactored `spls_analysis` (which pioneered the pattern in 1.56.0) to use them.
- Backward compatible: default `NA` resolves to samples-as-rows for existing input, so results are unchanged. `opls_analysis` applies orientation only to matrix/data.frame input (SummarizedExperiment / ExpressionSet keep their own path). `NAMESPACE` unchanged.
- Tests: new `tests/test_orientation.R` (helper unit tests + per-function equivalence: feature-row input reproduces the sample-row scores exactly); `tests/test_spls_analysis.R` updated to source the shared helper.

## [1.56.0] - 2026-06-29

### Added
`spls_analysis()` now accepts a feature × sample abundance matrix directly — no more manual `t() %>% as.data.frame()` before the call.

- New `feature_as_row` argument (default `NA` = auto-detect). With `sample` metadata supplied, the function matches sample IDs against both the row and column names of `data`: if they line up with the columns, the matrix is transposed internally. `feature_as_row = TRUE` / `FALSE` force a specific orientation; `NA` is the zero-argument default and resolves to samples-as-rows when `sample` is absent or the signal is ambiguous.
- Added internal helper `detect_feature_as_row()` (`@keywords internal`); `NAMESPACE` unchanged.
- Regression test added: `tests/test_spls_analysis.R` — feature-row input (auto-detected and forced) reproduces the sample-row scores exactly, and sample-row input is not spuriously transposed.

## [1.55.2] - 2026-06-29

### Fixed
`read_data()` no longer crashes on a `.xls`/`.xlsx` file whose content is actually delimited text.

- Symptom: `readxl::read_excel()` → `libxls error: Unable to open file`. Bioinformatics exporters (e.g. RNA-seq FPKM output) routinely emit tab/comma-delimited text with a `.xls` label — `readxl` (which expects a real Excel binary) rejects it.
- Fix: `read_data()` now sniffs the file's magic bytes before routing. Real Excel is detected via ZIP (`PK`, .xlsx) or OLE2 (`D0 CF 11 E0`, .xls) magic and still read by `readxl`; anything else is read as delimited text with the delimiter guessed from the first line. The `raw_names` attribute (original column names) is preserved on the text path too.
- Added internal helpers `is_excel_file()` / `guess_delim()` (`@keywords internal`); no change to the public API or `NAMESPACE`.
- Regression test added: `tests/test_read_data.R`.

## [1.55.1] - 2026-06-26

### Changed
Internal style cleanup per CLAUDE.md §3.3 — no behaviour or API change.

- `net2gephi.R`: two `=` assignments → `<-` (column names `Correlation`/`Direction` kept PascalCase; Gephi import requirement).
- Removed 23 trailing `return(x)` statements (a function body's last expression) across 9 files, replaced with the bare expression. An AST pass distinguished these from legitimate early-exit `return()` inside `if`/`else`/`tryCatch` success blocks, which were kept.

Files touched: `anova_posthoc.R`, `calc_expression_delta_delta_ct.R`, `calc_expression_qpcr_efficiency.R`, `calc_expression_standard_curve.R`, `df2fasta.R`, `extract_tree_hierarchy.R`, `get_lm_stats.R`, `identify_core_microbiome.R`, `net2gephi.R`, `pav_gwas.R`.

## [1.55.0] - 2026-06-25

### Changed (breaking)
Dev-standards naming sweep (per CLAUDE.md §1) — camelCase to snake_case:

- `microbiome_net()`: return list elements renamed. Migration: `adjaMats` → `adja_mats`, `nodeProps` → `node_props`, `globalStats` → `global_stats`, `edgeList` → `edge_list`, `spiecResults` → `spiec_results`. A previously undocumented `edge_list` element is now also described in `@return`.
- `microbiome_net.R` internal variables: `minSampCount` → `min_samp_count`, `X_clr` → `x_clr`, plus the per-group accumulators `nodeProps_list`/`globalStats_list`/`edgeList_list` → snake_case.
- `calc_expression_qpcr_efficiency.R`: `NF_new`/`NF_old`/`mean_M` → `nf_new`/`nf_old`/`mean_m`.
- `calc_alpha_diversity.R`: `needs_estimateR` → `needs_estimate_r`.
- Consumers `plot_microbiome_net.R` / `net2gephi.R` updated to read `net_result$node_props` / `net_result$adja_mats`.

### Preserved (backward compatibility)
- `microbiome_net()`: the previous camelCase element names (`adjaMats`, `nodeProps`, `globalStats`, `edgeList`, `spiecResults`) are retained as **deprecated aliases** pointing to the same objects, so existing scripts (`res$nodeProps`) keep working. Slated for removal in a future major version.
- Single-letter uppercase math variables (`X`, `M`, `N`, `V`, `R`) kept — R statistical convention.
- `pav_gwas()`: `library(stats)` inside `parallel::clusterEvalQ()` kept — standard pattern for loading deps on parallel worker processes, not a side-effect `library()` in ordinary code.
- `Cal*`/`CalExp*` legacy wrapper functions kept unchanged (CLAUDE.md §1.6 backward-compat exception).
- Gephi PascalCase columns in `net2gephi()` kept (format requirement).

## [1.54.0] - 2026-06-24

### Added
- `write_data_by_group()`: writes a grouped data frame to one file per group. Sits after `dplyr::group_by()` in a pipe; the file path is a `glue` template where each `{column}` is filled with that group's value (e.g. `"~/Downloads/result_{CHROM}.xlsx"`). Multi-column groups supported (`"{CHROM}_{TYPE}.csv"`). Delegates each group to `write_data()` unchanged, so all its formats (xlsx/csv/tsv/bed/txt/rds/images) are inherited; extra `...` args forward to `write_data()`. Added `glue` to Imports.

## [1.53.1] - 2026-06-24

### Fixed
Doc/example sweep to catch stragglers from the snake_case renames:

- `get_methylkit_data()`: params `lo.count`/`lo.perc`/`hi.count`/`hi.perc` → snake_case; passed to `methylKit::filterByCoverage` with original external names (`lo.count = lo_count`, etc.).
- `find_dams_lefse()` internal helper `try_lefser`: param `kruskal.threshold` → `kruskal_threshold`; passed to `lefser::lefser` as `krusal.threshold = kruskal_threshold`.
- @examples in `net2gephi.R`, `plot_upset.R`, `plot_microbiome_net.R`: updated calls to `microbiome_net(... group_col = ...)` and `find_dams_deseq2(... group_col = ...)` to use new param names.
- `bioRtools-troubleshooting.R`: updated `log2FoldChange = ` mentions in examples to `log2_fold_change = `; added clarifying comment for the `hist(all_results$log2FoldChange, ...)` line (kept as-is — that's a DESeq2 output column reference).

## [1.53.0] - 2026-06-23

### Changed (breaking)
Snake_case sweep extended to remaining output column names per DEVELOP_GUIDE 1.4.

- `enrich_go()`, `enrich_kegg()`: clusterProfiler-style column names converted to snake_case on output. Migration: `ID` → `id`, `Description` → `description`, `GeneRatio` → `gene_ratio`, `BgRatio` → `bg_ratio`, `p.adjust` → `p_adjust`, `geneID` → `gene_id`, `Count` → `count`, `gene.count` → `gene_count`, `total.genes` → `total_genes`, `enrichment.score` → `enrichment_score`. **Note**: this diverges from raw clusterProfiler output naming.
- `plot_ld_decay()`: internal column standardization now uses lowercase: `Population` → `population`, `Dist` → `dist`, `Mean_r2` → `mean_r2`. Default `pop_col`/`dist_col`/`value_col` also snake_cased (e.g., `"population"` instead of `"Population"`).
- `plot_microbiome_net()`: `Modularity` → `modularity_factor` (intermediate); `modularity2`/`modularity3` kept (already lowercase).
- `calc_expression_qpcr_efficiency()`: internal plotting column `Treatment` → `treatment`.
- `pav_gwas()`: internal plot data column `Variance` → `variance`.

### Preserved
- `net2gephi()`: Gephi requires PascalCase column names (`Id`, `Label`, `Source`, `Target`, `Weight`, `Type`) for CSV import. Renaming would break the function's purpose. Hard exception.

## [1.52.3] - 2026-06-23

### Changed (breaking)
More user-facing output column names snake_cased per DEVELOP_GUIDE 1.4:

- `permanova_test()`: `result_permanova` and `summary_stats` data frames: `Factor` → `factor`, `P_value` → `p_value`, `R_squared` → `r_squared`, `Significance` → `significance`
- `plot_microbiome_net()`: node data frame: `Degree` → `degree`, `Strength` → `strength`

## [1.52.2] - 2026-06-23

### Fixed
- `calc_standard_curve()`: tryCatch scope bug — error handler's `results_list[[i]] <- data.frame(...)` was inside a closure and the assignment never escaped to the parent frame, so failed genes were silently dropped (later filtered as NULL). Restructured: `result_i <- tryCatch({...}, error = function(e) {...}); results_list[[i]] <- result_i` so the assignment happens in the parent. Failed genes now properly recorded as a "Failed" row instead of being silently missing.

### Changed (breaking)
- `calc_standard_curve()`: output column names normalized to snake_case. Migration:
  - `Gene` → `gene`
  - `Formula` → `formula`
  - `Slope` → `slope`, `Intercept` → `intercept`
  - `R2` → `r2`, `P_value` → `p_value`
  - `Max_Cq` → `max_cq`, `Min_Cq` → `min_cq`
  - `Efficiency` → `efficiency`, `Efficiency_Percent` → `efficiency_percent`, `Efficiency_Quality` → `efficiency_quality`
  - `Method_Recommendation` → `method_recommendation`
  - `Date` → `date`

## [1.52.1] - 2026-06-23

### Changed (breaking)
Continued snake_case normalization sweep across remaining exported functions:

- `plot_LDheatmap` → `plot_ld_heatmap` (file renamed too: `R/plot_LDheatmap.R` → `R/plot_ld_heatmap.R`)
- `spls_analysis()`: `keepX` → `keep_x` (passed to mixOmics as `keepX = keep_x`)
- `pca_analysis()`, `tsne_analysis()`: `conf.ellipses` → `conf_ellipses`, `ellipse.level` → `ellipse_level`
- `find_outliers()`: `return.logical` → `return_logical`
- `find_dams_lefse()`: `groupCol` → `group_col`, `wilcox.threshold` → `wilcox_threshold`, `lda.threshold` → `lda_threshold`
- `plot_multi_volcano()`: `groupCol` → `group_col`, `featureCol` → `feature_col`, `log2FCCol` → `log2fc_col`, `padjCol` → `padj_col`, `log2FC_thr` → `log2fc_thr`, `topN` → `top_n`

### Verified clean
Full audit of all exported functions confirms parameters are snake_case (modulo documented exceptions: `na.rm`/`row.names` R base, `CalExp*` backward-compat, S3 methods like `print.core_microbiome`, internal `apply_bioRtools_theme` which references the package name).

## [1.52.0] - 2026-06-23

### Changed (breaking)
Normalize function parameter names to snake_case per DEVELOP_GUIDE 1.1, 1.3. For parameters passed to external packages, the call site preserves the external name (e.g. `pAdjustMethod = p_adjust_method`).

- `microbiome_net()`: `groupCol` → `group_col`, `minSamples` → `min_samples`, `minReads` → `min_reads`, `cor.threshold` → `cor_threshold`, `cor.pvalue` → `cor_pvalue`, `hubQuant` → `hub_quant`, `clustMethod` → `clust_method`, `pulsar.params` → `pulsar_params`
- `pairwise_oplsda()`: `groupCol` → `group_col`, `vip.threshold` → `vip_threshold`, `ortho.components` → `ortho_components`, `cv.folds` → `cv_folds`
- `enrich_go()`, `enrich_kegg()`: `pAdjustMethod` → `p_adjust_method`
- `find_degs_deseq2()`: `log2FoldChange` → `log2_fold_change`, `groupCol` → `group_col`, `shrink.lfc` → `shrink_lfc`, `independent.filtering` → `independent_filtering`
- `find_degs_edger()`: `log2FoldChange` → `log2_fold_change`
- `find_dams_deseq2()`: `log2FoldChange` → `log2_fold_change`, `groupCol` → `group_col`, `shrink.lfc` → `shrink_lfc`, `min.count` → `min_count`, `min.samples` → `min_samples`

### Preserved
- Column-name references to DESeq2 output (e.g. `res$log2FoldChange`, `df$padj`) kept as-is — these are external package output column names per DEVELOP_GUIDE 1.6.
- `CalExp*` family parameters (`cq.table`, `ref.gene`, etc.) — backward-compat functions.
- `na.rm`, `row.names`, `stringsAsFactors` — R base conventions.

## [1.51.1] - 2026-06-22

### Changed (breaking)
- `spls_analysis()`: returned list element names normalized to snake_case per DEVELOP_GUIDE 1.1 (same pattern as the 1.51.0 fix for `rda_analysis`). Migration:
  - `result$result.splsda` → `result$result_splsda`
  - `result$sample.scores` → `result$sample_scores`
  - `result$variance.explained` → `result$variance_explained`
  - `result$variable.loadings` → `result$variable_loadings`
  - `result$classification.performance` → `result$classification_performance`
  - `result$model.parameters` → `result$model_parameters`

## [1.51.0] - 2026-06-22

### Changed (breaking)
- `rda_analysis()`: returned list element names normalized to snake_case per DEVELOP_GUIDE 1.1. Migration:
  - `result$plot.rda` → `result$plot_rda`
  - `result$result.rda` → `result$result_rda`
  - `result$sample.scores` → `result$sample_scores`
  - `result$environmental.scores` → `result$environmental_scores`
  - `result$environmental.fit` → `result$environmental_fit`

### Fixed
- `rda_analysis()`: `environmental_significance` was assigned inside the error handler of `tryCatch` (a closure) and never escaped to the parent frame. If `vegan::envfit()` failed, the final `list()` would error with "object 'environmental_significance' not found". Now pre-defined as an empty data frame before `tryCatch` so the return value is well-formed in all cases.
- `theme_prism()`: `%+replace%` operator was used but not imported from ggplot2. Calling `theme_prism()` without `library(ggplot2)` failed with "没有 %>%+replace% 这个函数". Now imported via `@importFrom ggplot2 %+replace%`.

## [1.50.6] - 2026-06-22

### Changed
- `rda_analysis()`: default plot theme switched from `theme_bio()` to `theme_prism()`.

## [1.50.5] - 2026-06-22

### Fixed
- `pca_analysis()`: `score_plot` failed with `找不到对象'PC1'` because `sample_scores` columns were renamed to `"PC1 (XX.X%)"` format (including variance percent), but the plot's `aes()` referenced plain `PC1`/`PC2`. Column names now kept as `PC1, PC2, ...`. The variance percent is still displayed in plot axis labels and stored in `eigenvalue_data`.

## [1.50.4] - 2026-06-22

### Added
- `add_cor_p()`: per-group correlation and p-value helper. For each combination of grouping columns (specified via `from` + `to` or by existing `dplyr::group_by()`), computes `cor`, `p`, and `method` between two numeric columns and broadcasts to all rows in the group. Replaces the verbose `df %>% group_by(...) %>% mutate(cor = cor(x,y), p = cor.test(x,y)$p.value) %>% ungroup()` pattern. `from`/`to` make explicit which two entities are being correlated. `add_regression = TRUE` adds `lm_r2` and `lm_pvalue` from `lm(y ~ x)` on raw values (for Pearson, `lm_r2 == cor^2` and `lm_pvalue == p` by construction). The `method` column records the analysis method (e.g., `"pearson"` or `"pearson + lm"`) for paper Methods sections.

## [1.50.3] - 2026-06-22

### Changed
- `scale_fill_heatmap()` (and color/colour aliases): middle color default is now `#F0F0F0` (light gray) instead of pure `white`. Pure white was invisible on white plot backgrounds — middle-valued cells disappeared. New `mid_color` parameter lets users override (e.g., `mid_color = "white"` for original pheatmap style, or `mid_color = "yellow"` for RdYlBu style).

## [1.50.2] - 2026-06-22

### Added
- `scale_fill_heatmap()` (with `scale_color_heatmap()` and `scale_colour_heatmap()` aliases): pheatmap-style diverging blue-white-red heatmap fill scale, using the classic `colorRampPalette(c("#0C6291", "white", "#A63446"))` gradient. Parameters: `n` (number of interpolated colors, default 100), `alpha`, `reverse`, plus `...` forwarded to `ggplot2::scale_fill_gradientn()`. Best for diverging data centered at zero; for strictly positive sequential data use `scale_fill_sci_c()` with a sequential palette.

## [1.50.1] - 2026-06-22

### Fixed
- `scale_fill_research_c()` / `scale_color_research_c()` (and the `scale_colour_research_c` alias): calling with no arguments failed with `'arg' must be of length 1` because the full default `palette` vector (length 6) was passed to `scale_*_sci_c()`, whose `match.arg()` requires length 1. Added `palette <- match.arg(palette)` in the wrappers so the default selects the first palette (`research_teal_sequential`) and explicit palette names still work.

## [1.50.0] - 2026-06-22

### Changed (breaking)
- `anova_posthoc()`: output column previously hardcoded as `group` is now the original input name. For `anova_posthoc(iris, Sepal.Length ~ Species)`, output has `Species` instead of `group`. Same for string interface (`group = "treatment"` → output column `treatment`). The internal column was already renamed to `group_anova` for ANOVA fitting; only the final output column was renamed to a fixed `group`. This change preserves the input column name in the result so downstream code can refer to it by its actual name. **Migration**: code doing `result$group` should switch to `result[[<original_name>]]` or `result[[1]]` (the group column is always first).

## [1.49.9] - 2026-06-22

### Changed
- `name_map()`: renamed output column `sanitized` to `name` so it matches `tidyr::pivot_longer()`'s default key column. Enables direct `left_join()` without an intermediate `dplyr::rename()`.

## [1.49.8] - 2026-06-22

### Fixed
- `summarise_stats()`: passing `value` as a string (e.g., `summarise_stats("value")` or `summarise_stats(value = "value")`) silently failed because `{{ value }}` treated the string as a literal rather than a column reference. Switched to `.data[[value_col]]` so both bare names and strings work.

### Changed
- `summarise_stats()`: `value` now defaults to `"value"`, matching the conventional output column name from `tidyr::pivot_longer()`. The user's `pivot_longer() %>% group_by(...) %>% summarise_stats()` pipeline no longer needs to specify `value` explicitly.

## [1.49.7] - 2026-06-22

### Added
- `name_map()`: thin convenience function returning a two-column data frame (sanitized, raw) mapping current column names to originals stored by `read_data()`. Replaces the inline `data.frame(sanitized = names(df), raw = attr(df, "raw_names"))` pattern.

## [1.49.6] - 2026-06-22

### Added
- `read_data()`: for tabular formats (.xlsx/.xls/.csv/.tsv/.txt), the original column names from the file are now preserved as attribute `raw_names`. Use `attr(df, "raw_names")` to retrieve them. Useful when `readxl`/`readr` sanitizes names (e.g., `"Plant height (PH)"` -> `Plant.height..PH..`) — `raw_names` keeps the human-readable originals for plotting, reporting, or restoring before pivot_longer.

  ```r
  df <- read_data("phe.xlsx")
  attr(df, "raw_names")                    # original names
  names(df) <- attr(df, "raw_names")       # restore originals if needed
  ```

## [1.49.5] - 2026-06-22

### Added
- `pivot_longer_from()`: thin pipe-friendly wrapper around `tidyr::pivot_longer()` for the pattern `pivot_longer(cols = start:ncol(.))`. Pass just the starting column index; pivots from there to the last column. Extra args forwarded to `tidyr::pivot_longer()`.

## [1.49.4] - 2026-06-22

### Added
- `cor_analysis()`: new `add_regression` parameter (default FALSE). When TRUE, adds four columns to the result: `lm_slope`, `lm_intercept`, `lm_r2`, `lm_pvalue` from `lm(to ~ from)` on raw values. Computed only after correlation filtering, so cost scales with significant-pair count. For Pearson, `lm_r2 == cor^2` and `lm_pvalue == pvalue` by construction (verified numerically). For Spearman/Kendall, regression is independent of the rank-based correlation.

## [1.49.3] - 2026-06-18

### Changed
- `pairs_by_group()`: default `out_names` is now `c("id1", "id2")` (was `c("Gene1", "Gene2")`). The previous default was too biology-specific; `id1`/`id2` is generic.

## [1.49.2] - 2026-06-18

### Fixed
- `pairs_by_group()`: when `group_col` is passed as a character string (e.g., `group_col = "NLR类型"`), the previous `group_by({{ group_col }})` was treating the string as a literal and creating a constant column instead of grouping by the named column. Switched to `.data[[group_name]]` so both bare names and strings work. Also added clear errors when the group or id column is missing.

## [1.49.1] - 2026-06-18

### Added
- `pairs_by_group()`: for each group in a data frame, generate all C(n, 2) pairwise combinations of an ID column. Useful for KaKs calculation between gene pairs within clades, co-evolution, and within-group network construction. Replaces the verbose `group_by() %>% group_modify(~ combn(...))` pattern.

## [1.49.0] - 2026-06-18

### Changed
- Normalize local variables and output names to snake_case per DEVELOP_GUIDE. **User-facing list element renames (breaking):**
  - `rf_taxa_classification()`: `plot.accuracy` → `plot_accuracy`, `plot.top_features` → `plot_top_features`
  - `permanova_test()`: `result.permanova` → `result_permanova`, `summary.stats` → `summary_stats`, `raw.result` → `raw_result`
  - `pcoa_analysis()`: `result.pcoa` → `result_pcoa`, `plot.pcoa` → `plot_pcoa`, `point.data` → `point_data`, `eigenvalue.pcoa` → `eigenvalue_pcoa`
  - `plot_upset()`: `data.pav` → `data_pav`
- `plot_ld_decay()`: computed columns renamed to snake_case (`Dist_kb` → `dist_kb`, `Dist_bin_kb` → `dist_bin_kb`, `Mean_r2_bin` → `mean_r2_bin`, `Mean_r2_fitted` → `mean_r2_fitted`).
- Many internal local variables renamed across `rf_taxa_classification.R`, `admixture_phylo_analysis.R`, `permanova_test.R`, `plot_pav.R`, `pcoa_analysis.R`, `get_methylkit_data.R`, `anova_posthoc.R` (e.g., `df.admixture` → `df_admixture`, `data.new` → `data_new`). Internal-only, no API impact.

### Preserved (per DEVELOP_GUIDE exceptions)
- `na.rm`, `row.names`, `check.names`, `stringsAsFactors`, `ties.method` — R base conventions
- `branch.length`, `ref.group`, `override.aes`, `rep.num` — external package argument names
- `CalExp*` family (`cq.table`, `design.table`, `ref.gene`) — backward-compat functions
- `legend.position`, `panel.border`, `axis.line` in `theme_bio()`/`theme_prism()` — mirror ggplot2 argument names

## [1.48.4] - 2026-06-16

### Changed
- `colnames_to_df()`: default `name` is now `"sample"` (was `"name"`), matching the dominant use case in omics sample-metadata construction.

## [1.48.3] - 2026-06-16

### Added
- `colnames_to_df()`: extract column names of a data frame or matrix as a single-column data frame. Pipe-friendly replacement for `colnames(df) %>% as.data.frame() %>% magrittr::set_names("sample")`.

## [1.48.2] - 2026-06-15

### Fixed
- `find_degs_deseq2()`: reorder input cleanup so non-numeric columns (e.g., gene_id, description) are stripped to row names *before* the negative-value check. Previously the check ran first and errored with `仅在具有所有类似数值变量的数据框上定义` whenever the input carried a single character column, even though the cleanup logic existed to handle that case.

## [1.48.1] - 2026-06-15

### Added
- `plot_qq()`: new `group_column` parameter. When set, QQ stats are computed per group and the group column is preserved in every layer, enabling `+ facet_wrap(~ group)` directly on the returned plot. Per-group stats are attached as attribute `"group_stats"`.

## [1.48.0] - 2026-06-15

### Added
- `plot_qq()`: generalized QQ plot supporting both p-value mode (`type = "pvalue"`, mirrors `plot_gwas_qq`) and normal mode (`type = "normal"`, works with any numeric data such as transformed phenotypes). Adds `ci_level` parameter for adjustable confidence ribbon.

## [1.47.1] - 2026-06-15

### Changed
- `theme_nature()`, `save_pub()`: mark as internal (not exported)

## [1.47.0] - 2026-06-15

### Added
- `theme_nature()`: minimalist publication theme tuned for Nature-family figures (6.5pt Arial, 0.35pt lines, no grid)
- `save_pub()`: export plot as SVG + PDF + TIFF bundle in one call (default 183x120mm, 600 dpi)

## [1.46.4] - 2026-06-13

### Added
- `plot_pav_phenotype()`: add `show_sample_name` parameter (default FALSE)

## [1.46.3] - 2026-06-13

### Added
- `plot_pav_phenotype()`: PAV heatmap with phenotype annotation, samples clustered by PAV pattern

## [1.46.2] - 2026-06-13

### Changed
- `plot_ld_decay()`: half-decay lines colored by population; default `show_half_decay = FALSE`

## [1.46.1] - 2026-06-13

### Added
- `plot_ld_decay()`: LD decay plot with exponential/Hill-Weir fitting, half-decay distance calculation

## [1.46.0] - 2026-06-13

### Added
- `tsne_analysis()`: t-SNE dimensionality reduction with visualization, following same interface pattern as `pca_analysis()`

## [1.45.5] - 2026-06-12

### Added
- `microbiome_net()`: add `edgeList` to output with correlation, padj, direction, and taxonomy for each edge pair

## [1.45.4] - 2026-06-12

### Changed
- `microbiome_net()`: replace `Hmisc::rcorr()` with `WGCNA::corAndPvalue()` for faster correlation with no minimum sample restriction

## [1.45.3] - 2026-06-12

### Changed
- `microbiome_net()`: lower minimum sample count for correlation from 5 to 3

### Fixed
- `microbiome_net()`: crash when all groups are skipped (returns NULL with warning)
- `microbiome_net()`: pairwise comparison only uses successfully analyzed groups

## [1.45.2] - 2026-06-11

### Added
- `write_data()`: support `.fa` and `.fasta` formats via `df2fasta()`

## [1.45.1] - 2026-06-11

### Changed
- `plot_upset()`, `plot_tax_upset()`: add `engine` parameter (default "ggupset", alternative "ggVennDiagram"); restore ggupset parameters (`fill`, `show_counts`, `order_by`)

## [1.45.0] - 2026-06-11

### Changed
- `plot_upset()`, `plot_tax_upset()`: switch from ggupset to ggVennDiagram for upset plot generation

### Removed
- `plot_upset()`, `plot_tax_upset()`: removed `fill`, `show_counts`, `order_by` parameters (not applicable to ggVennDiagram); added `relative_height`, `relative_width` parameters

## [1.44.6] - 2026-06-11

### Changed
- `calc_beta_nti()`: add `sample_id` parameter (default "sample") to specify sample ID column; compute cophenetic matrix once and reuse across groups; filter zero-count taxa per group to reduce memory usage; add `gc()` between groups

## [1.44.5] - 2026-06-11

### Changed
- `calc_beta_nti()`: when `sample` + `group_col` provided, split by group and calculate βNTI independently within each group (much faster than computing all pairs)

## [1.44.4] - 2026-06-11

### Changed
- `calc_beta_nti()`: add `sample`, `group_col`, `ref_group` parameters for group-annotated βNTI analysis; fix picante compatibility issue by manual tree-data matching instead of `match.phylo.data`

## [1.44.3] - 2026-06-11

### Changed
- `write_data()`: PDF output defaults to `device = cairo_pdf` for better font and Unicode rendering

## [1.44.2] - 2026-06-11

### Changed
- `find_dams_deseq2()`: add `ref_group` parameter for reference-group comparisons; add `shrink.lfc` parameter; fix comparison label to "treatment vs reference" format; add group statistics columns (group_mean, group_sd, ref_mean, ref_sd, etc.)
- `find_degs_deseq2()`: add `groupCol` parameter as convenience shortcut to specify grouping column without writing a formula

## [1.44.1] - 2026-06-11

### Changed
- `write_data()`: support .bed format output (tab-delimited, no header)

## [1.44.0] - 2026-06-11

### Added
- `plot_ternary()`: ternary plot for three-part compositional data, supports discrete and continuous fill with auto-detection
- `plot_dual_axis()`: dual Y-axis plot combining area chart (right axis) and scatter points (left axis), supports phase annotations
- `find_degs_edger()`: differential expression analysis using edgeR

### Changed
- DESCRIPTION: added `edgeR` to Imports, `ggtern` to Suggests

## [1.43.2] - 2026-06-10

### Fixed
- `find_degs_deseq2()`: fix pairwise summary matching wrong regulation labels; fix comparison label order to treatment vs reference

## [1.43.1] - 2026-06-10

### Fixed
- `find_degs_deseq2()`: left_join original expression values (before rounding) instead of rounded integers

## [1.43.0] - 2026-06-10

### Added
- `itol_config()`: unified wrapper for iTOL annotation files via itol.toolkit, supports color strips, symbols, bar charts, and heatmaps

## [1.42.0] - 2026-06-10

### Fixed
- `find_degs_deseq2()`: fix gene rownames lost by `apply()` causing numeric gene column and NA in `group_mean`; fix pairwise stats index after DESeq2 zero-count gene filtering
- `pcoa_analysis()`, `spls_analysis()`, `opls_analysis()`: auto-detect sample ID column when specified column not found

### Changed
- `find_degs_deseq2()`: output now includes original expression values via left_join for downstream analysis

## [1.41.0] - 2026-06-10

### Fixed
- `find_degs_deseq2()`: auto-detect and remove non-numeric columns (e.g., gene_id) in data; auto-detect sample ID column in sample metadata when rownames are missing

## [1.40.0] - 2026-06-10

### Changed
- `find_degs_deseq2()`: add `ref_group` parameter; output now includes `group`, `ref_group`, `group_mean`, `group_sd`, `group_n`, `ref_mean`, `ref_sd`, `ref_n`, `test_method` columns aligned with `opls_analysis()` format

## [1.39.0] - 2026-06-10

### Added
- `plot_pav()`: pangenome presence/absence matrix with category annotation bar

### Fixed
- `microbiome_net()`: use `Hmisc::rcorr()` for fast correlation; skip groups with < 5 samples

## [1.38.0] - 2026-06-10

### Changed
- `microbiome_net()`: new `cor_method` parameter ("spearman", "pearson", "kendall") for correlation-based network construction

## [1.37.0] - 2026-06-10

### Added
- `plot_venn()`: unified Venn/Euler diagram function, supports proportional area (eulerr) and classic equal-size (ggvenn) modes, input is group + feature columns

## [1.36.1] - 2026-06-10

### Fixed
- Add `@examples` to `quick_pfam_plot` documentation
- Update man pages for opls_analysis new parameters

## [1.36.0] - 2026-06-10

### Changed
- `opls_analysis()`: p-value correction now applied per comparison instead of across all; new `p_adjust_method` parameter (default "BH")

## [1.35.0] - 2026-06-10

### Changed
- `rda_analysis()`: add `variance_RDA1` and `variance_RDA2` columns to sample.scores data frame

## [1.34.0] - 2026-06-10

### Changed
- `opls_analysis()`: new `test_method` parameter ("auto", "t-test", "wilcoxon"); "auto" uses shapiro.test() per variable to decide normality

## [1.33.1] - 2026-06-10

### Added
- `opls_analysis()`: add `test_method` column to differential_analysis ("wilcoxon" or "t-test")

## [1.33.0] - 2026-06-10

### Changed
- `opls_analysis()`: add R2X_p1, R2X_o1, R2X_cum, R2Y_cum, Q2_cum as columns in scores data frame (instead of attributes), so all plotting data is in one data frame

## [1.32.0] - 2026-06-10

### Changed
- `opls_analysis()`: attach ropls plot metrics (R2X_p1, R2X_o1, R2X_cum, R2Y_cum, Q2_cum) to scores as `plot_metrics` attribute, for recreating ropls-style score plots; suppress auto-plot output

## [1.31.0] - 2026-06-10

### Changed
- `opls_analysis()`: print ropls diagnostic plots directly (was suppressed); pairwise mode prints each model with "X vs Y" header

## [1.30.0] - 2026-06-10

### Added
- `opls_analysis()`: return `plot` element — ggplot2 score plot with R2Y/Q2Y variance on axes, faceted by comparison in pairwise mode
- `opls_analysis()`: differential analysis now includes `group_mean`, `group_sd`, `group_n`, `ref_mean`, `ref_sd`, `ref_n` columns

## [1.29.4] - 2026-06-09

### Added
- `opls_analysis()`: add `comparison` column to differential_analysis output (e.g. "Cd vs CK")

## [1.29.3] - 2026-06-09

### Fixed
- `opls_analysis()`: restore missing `n_important_vars` computation lost in previous refactor

## [1.29.2] - 2026-06-09

### Fixed
- `opls_analysis()`: pairwise mode differential analysis was overwritten to NULL by shared code block

## [1.29.1] - 2026-06-09

### Fixed
- `opls_analysis()`: remove duplicate `group` column before joining sample_info to avoid `.x`/`.y` suffixes in scores output
- Fix stray `?` in `@export` tags in 8 R files causing roxygen2 lint errors

## [1.29.0] - 2026-06-09

### Added
- `opls_analysis()`: pairwise OPLS-DA mode — when `ref_group` is specified with >2 groups, automatically fits separate binary OPLS-DA models (each group vs ref_group), combines VIP scores and differential analysis across all comparisons
- `opls_analysis()`: new `verbose` parameter to print progress messages during pairwise fitting

## [1.28.0] - 2026-06-09

### Added
- `opls_analysis()`: new `ref_group` parameter for one-vs-reference differential analysis with log2FC, p-value, BH adjustment, and Up/Down regulation labeling
- `spls_analysis()`: improved error message when `group` is accidentally passed as data frame

## [1.27.1] - 2026-06-09

### Fixed
- Fix extra closing parenthesis in 12 files caused by `return()` removal
- Restore accidentally deleted NAMESPACE file
- `read_data()`: add `.name_repair = "universal"` for Excel files to handle empty/duplicate column names

## [1.27.0] - 2026-06-09

### Changed
- `read_data()`: add `delim` parameter for `.txt` files (default `"\t"`), add file existence check, improve `@return` documentation

## [1.26.0] - 2026-06-09

### Changed
- Complete roxygen2 documentation for all exported functions: add missing `@author` (55 files), `@return` (academic_color_scales), `@examples` (plot_manhattan)
- Unify `@author` format to `Xiang LI <lixiang117423@gmail.com>` across all functions

## [1.25.0] - 2026-06-09

### Changed
- Enforce DEVELOP_GUIDE.md compliance across all R source files (45 files touched)
- Remove `library()` calls in function bodies, replace with `pkg::fun()` namespace (pav_gwas, plot_LDheatmap, scale01, row_stat, run_wgcna)
- Fix non-snake_case output column names: `PCoA`→`pcoa`, `mean.sample`→`mean_sample`, `Run/Group/OTU`→`run/group/otu`, `Treatment`→`treatment`
- Fix non-snake_case parameter names: `plot_mantel` (`data.spec`→`data_spec`, `data.env`→`data_env`, `data.sample`→`data_sample`)
- Remove unnecessary `return()` at end of functions across all exported functions

## [1.24.0] - 2026-06-09

### Changed
- Rewrite `DEVELOP_GUIDE.md` as definitive development standard (v2.0), based on tidyverse snake_case conventions


## [1.23.0] - 2026-06-08

### Changed
- Standardize all exported function parameter names to snake_case (affects 10 functions: `calc_beta_nti`, `cor_analysis`, `enrich_kegg`, `enrich_go`, `lm_analysis`, `manhattan_plot`, `pca_analysis`, `opls_analysis`, `anova_posthoc`, `rf_taxa_classification`)
- Standardize output column names to snake_case (`anova.pvalue`→`anova_pvalue`, `Tukey.signif`→`tukey_signif`, `Duncan.signif`→`duncan_signif`, `tax.group`→`tax_group`)


## [1.22.0] - 2026-06-08

### Changed
- Standardize qPCR function output column names to snake_case (`expre`→`relative_expression`, `Expre4Stat`→`expression_value`, `Expression`→`mean_expression`, `SD`→`sd_expression`, `SE`→`se_expression`, `signif`→`significance`, etc.)
- Fix `create_ddct_plot` box/bar plot duplication (box now shows individual data points)
- Replace deprecated `tidyr::spread` with `tidyr::pivot_wider`
- Replace deprecated ggplot2 `size` with `linewidth` for line elements
- Rename data files for naming consistency: `df.call_DAMs_LEfSe.*`→`df.call_dams_lefse.*`, `df.rnaseq.plot_volcano`→`df.rnaseq.plot.volcano`
- Remove redundant bottom-of-file comments from qPCR functions


## [1.21.0] - 2026-06-08

### Added
- `summarise_stats`: grouped descriptive statistics (max, min, mean, median, sd, se, cv, iqr, range, sum, n, n_na)
- `calc_beta_nti`: betaNTI for phylogenetic beta diversity analysis
- `net2gephi`: export network to Gephi CSV format


## [1.20.3] - 2026-06-08

### Added
- `net2gephi`: export microbiome_net result as Gephi-compatible nodes/edges CSV files


## [1.20.2] - 2026-06-08

### Added
- `microbiome_net`: add `method = "cor"` option (Spearman correlation, fast and low memory)
- `microbiome_net`: warn when features/samples ratio > 100:1


## [1.20.1] - 2026-06-08

### Added
- `net2ggnetview`: convert microbiome_net result to ggNetView tbl_graph for visualization


## [1.20.0] - 2026-06-08

### Added
- `microbiome_net`: microbiome network analysis via SpiecEasi (huge + pulsar)
  - CLR transformation + sparse inverse covariance estimation
  - Network properties: centrality, hubs, clusters, global stats
  - Multi-group batch construction and pairwise comparison


## [1.19.4] - 2026-06-08

### Changed
- `write_data`: add col_names parameter (default TRUE), apply to csv/tsv/txt


## [1.19.3] - 2026-06-08

### Fixed
- `write_data`: fix switch syntax error for image formats
- `write_data`: expose width/height/dpi as explicit parameters (default: 8/6/600)


## [1.19.2] - 2026-06-08

### Changed
- `write_data`: add image format support (.pdf/.png/.svg/.tiff/.jpg/.eps) via ggplot2::ggsave


## [1.19.1] - 2026-06-08

### Added
- `replace_na_as`: replace all NA values in a data frame (default: 0)


## [1.19.0] - 2026-06-08

### Added
- `read_data`: read files by extension (.xlsx, .csv, .tsv, .txt, .fasta, .fa, .rds)
- `write_data`: write files by extension (.xlsx, .csv, .tsv, .txt, .sh, .rds)
- `ggsave2`: added dpi=600 default


## [1.18.0] - 2026-06-08

### Added
- `pairwise_oplsda`: pairwise OPLS-DA across all group combinations
- `col2file`: export a data frame column to plain text file (default append=TRUE)
- `ggsave2`: ggsave wrapper with width=8, height=6, dpi=600 defaults
- `find_degs_deseq2`: add pairwise=TRUE for automatic pairwise DESeq2 analysis

### Fixed
- `fasta2df`: strip leading '>' from ID column
- `pairwise_oplsda`: scores output changed from list to data.frame with comparison column
- Removed obsolete `df_to_list`, cleaned up man page duplicates and typos


## [1.17.0] - 2026-06-07

### Added
- `pairwise_oplsda`: pairwise OPLS-DA across all group combinations (later refined: scores as data.frame)

### Fixed
- `fasta2df`: strip leading '>' from ID column
- Removed obsolete `df_to_list`, cleaned up man page duplicates and typos


## [1.16.1] - 2026-05-20

### Changed
- update theme_prism
- Updated files: R/theme_prism.R


## [1.16.0] - 2026-05-20

### Changed
- add theme_prism with visible legend title
- Updated files: R/theme_prism.R


## [1.15.7] - 2026-05-20

### Changed
- fix anova_posthoc Tukey letter order: assign 'a' to highest mean group
- Updated files: R/anova_posthoc.R


## [1.15.6] - 2026-05-20

### Changed
- add 'default' 10-color high-contrast palette as research default
- Updated files: R/academic_color_scales.R,man/pal_sci.Rd,man/scale_color_sci.Rd


## [1.15.5] - 2026-05-20

### Changed
- fix anova_posthoc Tukey letter order: assign 'a' to highest mean group
- Updated files: R/anova_posthoc.R


## [1.15.4] - 2026-05-20

### Changed
- add 'default' 10-color high-contrast palette for research scales
- Updated files: R/academic_color_scales.R


## [1.15.3] - 2026-05-19

### Changed
- update find_outliers with auto quantile type selection (type 2 for small samples, type 7 for large)
- Updated files: R/find_outliers.R


## [1.15.2] - 2026-05-18

### Changed
- remove dev-only files from repo; refactor multiple R modules
- Updated files: _pkgdown.yml,.claude,.DS_Store,.git,.gitignore,.Rbuildignore,.rdevrc,backup_original_files,bioRtools.Rproj,CHANGELOG.md,CODE_REVIEW_REPORT.md,data,DESCRIPTION,deve.R.backup,DEVELOP_GUIDE.md,inst,INSTALL_SYNTENY.md,LICENSE,LICENSE.md,man,missing_author_tags.csv,NAMESPACE,PlotCase,R,README_cn.md,README_en.md,README.md,Rplots.pdf,scripts,SETUP.md,tests,vignettes


## [1.15.1] - 2026-05-18

### Changed
- set coral_teal as default palette for pal_research and update plot_manhattan default colors
- Updated files: R/academic_color_scales.R


## [1.15.0] - 2026-05-18

### Changed
- add pal_research with 7 discrete palettes and 6 gradient scales; update plot_manhattan colors to coral/teal
- Updated files: R/academic_color_scales.R,R/plot_manhattan.R,NAMESPACE,man


## [1.14.1] - 2026-04-27

### Changed
- update plot_manhattan with density plot
- Updated files: R/plot_manhattan.R


## [1.14.0] - 2026-03-27

### Changed
- add get_hap_heatmap module
- Updated files: R/get_hap_from_heatmap.R


## [1.13.0] - 2026-03-25

### Changed
- add ld_decay_threshold module
- Updated files: R/ld_decay_threshold.R


## [1.12.1] - 2026-03-18

### Changed
- upload find_degs_deseq2 module
- Updated files: R/find_degs_deseq2.R


## [1.12.0] - 2026-03-06

### Changed
- add plot_gwas_qq module
- Updated files: R/plot_gwas_qq.R


## [1.11.0] - 2026-02-09

### Changed
- add get_lm_stats_summary  create_heatmap_trees plot_gene_structure


## [1.10.5] - 2026-02-09

### Changed
- update library message


## [1.10.4] - 2026-02-09

### Changed
- move normentR to Suggests


