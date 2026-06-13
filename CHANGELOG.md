# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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


