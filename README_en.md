# bioRtools

[![R](https://img.shields.io/badge/R-%3E%3D2.10-blue)](https://www.r-project.org/)
[![Version](https://img.shields.io/badge/version-1.16.1-green)](https://github.com/lixiang117423/bioRtools)
[![License](https://img.shields.io/badge/license-MIT-yellow)](LICENSE.md)

`bioRtools` is an R package that collects convenience functions for biological data analysis, statistics, and publication-oriented visualization. It covers common workflows for transcriptomics, microbiome analysis, metabolomics, population genetics, gene structure visualization, qPCR analysis, and reusable plotting themes.

## Key Features

- Multi-omics analysis helpers for transcriptomics, microbiomics, metabolomics, and population genetics
- Multivariate workflows: PCA, PCoA, RDA, sPLS-DA, OPLS-DA, PERMANOVA
- Differential analysis: DESeq2-based DEG/DAM detection and LEfSe biomarker analysis
- qPCR workflows: standard curves, delta Ct, delta-delta Ct, and efficiency correction
- Genomics visualization: Manhattan plots, QQ plots, LD heatmaps, synteny, motifs, gene structures, PFAM domains, and pangenome rarefaction
- Publication-ready ggplot themes, academic palettes, and color/fill/colour scales

## Requirements

- R >= 2.10, as declared in `DESCRIPTION`
- R >= 4.0.0 is recommended for a smoother dependency experience

## Installation

Install the development version from GitHub:

```r
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

devtools::install_github("lixiang117423/bioRtools")
```

If dependency installation needs to be handled manually, install the main Bioconductor and CRAN dependencies first:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "clusterProfiler", "DESeq2", "ggtree", "lefser",
  "SummarizedExperiment", "WGCNA"
))

install.packages(c(
  "agricolae", "ape", "broom", "data.table", "dplyr",
  "factoextra", "FactoMineR", "forcats", "ggplot2", "ggprism",
  "ggrepel", "ggsci", "ggtext", "ggthemes", "janitor", "magrittr",
  "mixOmics", "multcomp", "patchwork", "purrr", "rlang", "ropls",
  "rstatix", "scales", "stringr", "tibble", "tidyr", "vegan"
))
```

Optional suggested packages include `chemhelper`, `ggpmisc`, `knitr`, `normentR`, `rmarkdown`, and `testthat`.

## Function Modules

### General Analysis and Statistics

- `pca_analysis()` - Principal component analysis
- `cor_analysis()` - Correlation analysis and heatmap-ready output
- `lm_analysis()`, `get_lm_stats()`, `get_lm_stats_summary()`, `extract_lm_stats()`, `format_lm_stats()` - Linear model analysis helpers
- `anova_posthoc()` - ANOVA with post-hoc tests
- `find_outliers()` - Outlier detection
- `label_signif()`, `label_significance()` - Significance labels
- `reorder_heatmap()`, `gg_heatmap()`, `create_heatmap_trees()` - Heatmap utilities

### Microbiome Analysis

- `pcoa_analysis()` - Principal coordinate analysis
- `rda_analysis()` - Redundancy analysis
- `permanova_test()` - PERMANOVA test
- `find_dams_deseq2()` - DESeq2-based differential abundance analysis
- `find_dams_lefse()` - LEfSe differential abundance analysis
- `rarefy_table()` - Rarefaction
- `top_taxa()` - Top taxa extraction
- `identify_core_microbiome()`, `fit_sloan_neutral_model()` - Core microbiome and neutral model analysis

### Transcriptomics and Enrichment

- `find_degs_deseq2()` - DESeq2 differential expression analysis
- `enrich_go()` - GO enrichment
- `enrich_kegg()` - KEGG enrichment
- `plot_volcano()`, `plot_multi_volcano()` - Volcano plots
- `run_wgcna_analysis()` - WGCNA workflow

### qPCR Analysis

- `calc_standard_curve()` - Standard curve fitting
- `calc_expression_qpcr_efficiency()` - qPCR efficiency calculation
- `calc_expression_delta_ct()` - Delta Ct analysis
- `calc_expression_delta_delta_ct()` - Delta-delta Ct analysis
- `calc_expression_standard_curve()` - Standard-curve-based expression analysis

### Metabolomics and Multivariate Models

- `opls_analysis()` - OPLS-DA
- `spls_analysis()` - sPLS-DA

### Population Genetics and Genomics

- `manhattan_plot()`, `plot_manhattan()` - Manhattan plots
- `plot_gwas_qq()` - GWAS QQ plot
- `plot_LDheatmap()` - LD heatmap
- `ld_decay_threshold()` - LD decay threshold calculation
- `admixture_phylo_analysis()`, `extract_tree_hierarchy()` - Population structure and tree helpers
- `pav_gwas()` - PAV-GWAS helper

### Gene, Motif, Domain, and Synteny Visualization

- `plot_gene_structure()`, `plot_gene_features()`, `plot_gene_features_labeled()` - Gene structure and feature plots
- `plot_motif_location()`, `get_motif_from_meme()` - Motif parsing and visualization
- `plot_pfam()`, `quick_pfam_plot()` - PFAM domain plots
- `plot_synteny()` - Synteny plot
- `plot_pangenome_rarefaction()` - Pangenome rarefaction plot
- `get_hap_from_heatmap()` - Haplotype extraction from heatmap-like data

### Data Conversion and Row Utilities

- `df_to_list()` - Convert data frame to list
- `df2fasta()`, `fasta2df()` - FASTA/data-frame conversion
- `get_methylkit_data()` - Extract methylKit data
- `normalize_int()`, `scale01()`, `scale01_rows()`, `scale01_groups()`, `mutate_scale01()`, `mutate_scale01_named()` - Normalization helpers
- `row_mean()`, `row_sd()`, `row_cv()`, `row_min()`, `row_max()` - Row statistics
- `filter_str()`, `dedup_by_col()` - Regex/string row filtering (supports negation) and duplicate removal by column

### Themes and Palettes

- `theme_bio()`, `theme_prism()` - Plot themes
- `pal_sci()`, `pal_nature()`, `pal_science()`, `pal_cell()`, `pal_jacs()`, `pal_fuel()`, `pal_chem_eng()`, `pal_nat_comm()`, `pal_shinkai()`, `pal_research()` - Academic palettes
- `scale_color_*()`, `scale_colour_*()`, and `scale_fill_*()` variants are available for discrete and continuous palette scales

## Quick Start

### PCA

```r
library(bioRtools)

iris_data <- iris[, 1:4]
sample_info <- data.frame(
  sample_id = paste0("sample_", seq_len(nrow(iris))),
  species = iris$Species
)
rownames(iris_data) <- sample_info$sample_id

pca_result <- pca_analysis(
  data = iris_data,
  sample = sample_info,
  color_by = "species"
)

print(pca_result$plots$score_plot)
print(pca_result$eigenvalues)
```

### Correlation Analysis

```r
cor_result <- cor_analysis(
  data_1 = iris[, 1:2],
  data_2 = iris[, 3:4],
  method = "pearson"
)

print(head(cor_result))
```

### LEfSe Differential Abundance

```r
data(df.call_DAMs_LEfSe.otu)
data(df.call_DAMs_LEfSe.sample)

lefse_result <- find_dams_lefse(
  data = df.call_DAMs_LEfSe.otu,
  sample = df.call_DAMs_LEfSe.sample,
  group_col = "group",
  lda_threshold = 1.0
)

head(lefse_result)
```

### Volcano Plot

```r
data(df.rnaseq.plot_volcano)

volcano_result <- plot_volcano(
  data = df.rnaseq.plot_volcano,
  x = "log2FoldChange",
  y = "padj",
  title = "Differential Expression"
)

print(volcano_result$plot_volcano)
print(volcano_result$data_summary)
```

## RNA-seq Workflow Example

```r
perform_rnaseq_analysis <- function(count_data, sample_data, go_db, kegg_db) {
  degs <- find_degs_deseq2(
    data = count_data,
    sample = sample_data,
    formula = ~group
  )

  sig_genes <- degs$gene[degs$padj < 0.05 & degs$regulation != "NS"]

  go_enrichment <- enrich_go(gene = sig_genes, go_db = go_db)
  kegg_enrichment <- enrich_kegg(gene = sig_genes, kegg_db = kegg_db)

  volcano_result <- plot_volcano(
    data = degs,
    x = "log2FoldChange",
    y = "padj"
  )

  list(
    degs = degs,
    go_enrichment = go_enrichment,
    kegg_enrichment = kegg_enrichment,
    volcano_plot = volcano_result$plot_volcano
  )
}
```

## Data Format Requirements

### Expression, Count, or Abundance Matrix

- Rows are genes, features, taxa, OTUs, ASVs, or markers
- Columns are samples
- Row names and column names should be stable identifiers
- DESeq2 workflows require raw non-negative integer counts

### Sample Metadata

- Rows are samples
- Columns are experimental factors, grouping variables, batches, and covariates
- Sample identifiers must match the input matrix sample names

### Enrichment Annotation Tables

- GO annotation tables should contain `gene`, `go_id`, and `go_term`
- KEGG annotation tables should contain `gene`, `kegg_id`, and `kegg_term`

## Troubleshooting

### Dependency Installation Fails

```r
update.packages(ask = FALSE)
install.packages("BiocManager")
BiocManager::install(ask = FALSE)
```

### Memory Pressure

```r
rm(list = ls())
gc()
```

On Windows, you can also raise the memory limit in older R versions:

```r
memory.limit(size = 16000)
```

### Font Display Issues

```r
theme_bio(base_family = "Arial")
theme_prism(base_family = "Arial")
```

## Version History

See [CHANGELOG.md](CHANGELOG.md) for the full release history. Current package version: `1.16.1`.

## Citation

If you use `bioRtools` in your research, please cite:

```text
Li, X. (2026). bioRtools: Convenience Functions for Biological Data Processing.
R package version 1.16.1. https://github.com/lixiang117423/bioRtools
```

## Contributing

Contributions are welcome:

1. Report bugs in [GitHub Issues](https://github.com/lixiang117423/bioRtools/issues)
2. Suggest features through Issues
3. Fork the repository, create a feature branch, commit your changes, push the branch, and open a pull request

## License

This project is licensed under the MIT License. See [LICENSE.md](LICENSE.md) for details.

## Author

**Xiang LI**

- Maintainer
- Email: lixiang117423@gmail.com
- GitHub: [@lixiang117423](https://github.com/lixiang117423)

## Links

- GitHub Repository: <https://github.com/lixiang117423/bioRtools>
- Documentation: <https://lixiang117423.github.io/bioRtools/>
- Bug Reports: <https://github.com/lixiang117423/bioRtools/issues>
