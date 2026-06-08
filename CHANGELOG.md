# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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


