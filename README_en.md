# bioRtools 🧬

[![R](https://img.shields.io/badge/R-%3E%3D2.10-blue)](https://www.r-project.org/) [![Version](https://img.shields.io/badge/version-1.3.0-green)](https://github.com/lixiang117423/bioRtools) [![License](https://img.shields.io/badge/license-MIT-yellow)](https://claude.ai/chat/LICENSE.md)

---

### 📋 Introduction

`bioRtools` is a comprehensive R package designed for biological data processing, providing a suite of convenience functions for biological data analysis and visualization. The package standardizes complex analytical workflows across omics data types and provides consistent, publication-quality visualization outputs with minimal code requirements.

### ✨ Key Features

- **🧬 Multi-omics Support**: Transcriptomics, metabolomics, microbiomics, population genetics
- **📊 Multivariate Analysis**: PCA, PCoA, RDA, sPLS-DA, OPLS-DA
- **🔬 Statistical Analysis**: ANOVA, correlation analysis, linear regression, PERMANOVA
- **📈 High-Quality Visualization**: Publication-ready plots including volcano plots, Manhattan plots
- **⚡ Easy to Use**: Simplified function interfaces with standardized output formats
- **🎨 Beautiful Themes**: Built-in themes for biological data visualization

### 💻 System Requirements

- **R Version**: R (≥ 2.10)
- **Recommended**: R ≥ 4.0.0 for optimal performance

### 📦 Installation

#### Install Development Version from GitHub (Recommended)

```r
# Install required dependencies
if (!require(devtools)) install.packages("devtools")

# Install bioRtools
devtools::install_github("lixiang117423/bioRtools")
```

#### Manual Dependency Installation

If you encounter dependency issues, install manually:

```r
# Install Bioconductor packages
if (!require(BiocManager)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "clusterProfiler", "lefser", "SummarizedExperiment"))

# Install CRAN packages
install.packages(c("ggplot2", "dplyr", "vegan", "mixOmics", "ropls", 
                   "factoextra", "FactoMineR", "ggsci", "ggprism", "ggrepel",
                   "rstatix", "broom", "scales", "tidyr"))
```

### 🛠️ Main Function Modules

#### 🔧 General Analysis Functions

- `pca_analysis()` - Principal Component Analysis (PCA)
- `cor_analysis()` - Correlation analysis and visualization
- `lm_analysis()` - Linear regression analysis and visualization
- `find_outliers()` - Outlier detection
- `theme_bio()` - Biological data visualization theme
- `reorder_heatmap()` - Heatmap data reordering
- `anova_posthoc()` - ANOVA with post-hoc tests

#### 🧪 Metabolomics Analysis

- `opls_analysis()` - Orthogonal Partial Least Squares Discriminant Analysis (OPLS-DA)
- `spls_analysis()` - Sparse Partial Least Squares Discriminant Analysis (sPLS-DA)

#### 🦠 Microbiome Analysis

- `pcoa_analysis()` - Principal Coordinate Analysis (PCoA)
- `rda_analysis()` - Redundancy Analysis (RDA)
- `find_dams_deseq2()` - DESeq2-based differential abundance analysis
- `find_dams_lefse()` - LEfSe differential analysis
- `permanova_test()` - PERMANOVA permutation test
- `rarefy_table()` - Rarefaction sampling
- `top_taxa()` - Extract top taxa

#### 🧬 Transcriptomics Analysis

- `find_degs_deseq2()` - DESeq2 differential expression analysis
- `enrich_go()` - GO enrichment analysis
- `enrich_kegg()` - KEGG pathway enrichment analysis
- `plot_volcano()` - Volcano plot generation

#### 🧮 Population Genetics Analysis

- `manhattan_plot()` - Manhattan plot generation
- `admixture_phylo_analysis()` - Population structure and phylogenetic analysis
- `plot_LDheatmap()` - LD heatmap plotting

#### 🔧 Utility Functions

- `df_to_list()` - Convert data frame to list
- `plot_manhattan()` - Generic Manhattan plot function
- `scale01()` - Data normalization functions
- `row_mean()`, `row_sd()` - Row-wise statistical functions

### 🚀 Quick Start

#### Principal Component Analysis (PCA)

```r
library(bioRtools)

# Prepare data
data <- iris[,1:4]
sample_info <- data.frame(
  sample = paste0("sample", 1:150),
  species = iris$Species
)

# Perform PCA analysis
pca_result <- pca_analysis(
  data = data, 
  sample = sample_info,
  color.by = "species"
)

# View results
print(pca_result$plots$score_plot)
print(pca_result$eigenvalues)
```

#### Correlation Analysis

```r
# Correlation analysis between two datasets
cor_result <- cor_analysis(
  data.1 = iris[,1:2], 
  data.2 = iris[,3:4],
  method = "pearson"
)

# View correlation heatmap
print(cor_result$plot.cor)
```

#### LEfSe Differential Analysis

```r
# Microbiome differential analysis example
lefse_result <- find_dams_lefse(
  data = abundance_matrix,           # Feature abundance matrix
  sample = sample_metadata,          # Sample information
  groupCol = "treatment",            # Group column name
  lda.threshold = 2.0               # LDA threshold
)

# View significantly different features
print(head(lefse_result))
```

#### Volcano Plot

```r
# Generate volcano plot for differential expression
volcano_plot <- plot_volcano(
  data = deg_results,               # Differential analysis results
  x = "log2FoldChange", 
  y = "padj",
  label = "gene"
)

print(volcano_plot)
```

### 📊 Data Format Requirements

#### Expression/Abundance Matrix

- **Rows**: Genes/Features/OTUs
- **Columns**: Samples
- **Values**: Raw counts, normalized expression, or relative abundance

#### Sample Information Table

- **Rows**: Samples
- **Columns**: Experimental factors and covariates
- **Requirement**: Sample names must match expression matrix column names

#### Annotation Database Format

- **GO Database**: Must contain `gene`, `go.id`, `go.term` columns
- **KEGG Database**: Must contain `gene`, `kegg.id`, `kegg.term` columns

### 🎨 Advanced Usage

#### Custom Visualization Themes

```r
library(ggplot2)

# Using bioRtools theme
p <- ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width)) +
  geom_point(aes(color = Species), size = 2) +
  theme_bio(base_size = 12) +
  labs(title = "Scatter Plot with bioRtools Theme",
       x = "Sepal Length", y = "Sepal Width")

print(p)
```

#### Batch Analysis Workflow

```r
# Complete RNA-seq differential analysis workflow
perform_rnaseq_analysis <- function(count_data, sample_data, go_db, kegg_db) {
  
  # 1. Differential expression analysis
  cat("🔬 Performing differential expression analysis...\n")
  degs <- find_degs_deseq2(
    count.table = count_data,
    sample.table = sample_data,
    design = ~ condition,
    contrast = c("condition", "treatment", "control")
  )
  
  # 2. GO enrichment analysis
  cat("📊 Performing GO enrichment analysis...\n")
  sig_genes <- degs$degs$gene[degs$degs$padj < 0.05]
  go_enrichment <- enrich_go(gene = sig_genes, go.db = go_db)
  
  # 3. KEGG enrichment analysis
  cat("🧬 Performing KEGG enrichment analysis...\n")
  kegg_enrichment <- enrich_kegg(gene = sig_genes, kegg.db = kegg_db)
  
  # 4. Generate volcano plot
  cat("📈 Generating visualization plots...\n")
  volcano_plot <- plot_volcano(
    data = degs$degs,
    x = "log2FoldChange",
    y = "padj"
  )
  
  return(list(
    degs = degs,
    go_enrichment = go_enrichment,
    kegg_enrichment = kegg_enrichment,
    volcano_plot = volcano_plot
  ))
}
```

### ⚠️ Troubleshooting

#### Common Issues and Solutions

**Issue 1**: Dependency installation failure

```r
# Solution: Update R and package managers
update.packages(ask = FALSE)
install.packages("BiocManager")
```

**Issue 2**: Memory insufficient error

```r
# Solution: Increase memory limit (Windows)
memory.limit(size = 16000)
# Or clean workspace
rm(list = ls())
gc()
```

**Issue 3**: Font display issues

```r
# macOS
theme_bio(base_family = "Arial")
# Windows  
theme_bio(base_family = "Arial")
# Ubuntu
theme_bio(base_family = "DejaVu Sans")
```

### 🤝 Contributing

We welcome contributions of all kinds!

1. **🐛 Report Bugs**: Report in [GitHub Issues](https://github.com/lixiang117423/bioRtools/issues)

2. **💡 Feature Requests**: Suggest new features through Issues

3. 💻 Code Contributions

   :

   - Fork the repository
   - Create a feature branch (`git checkout -b feature/amazing-feature`)
   - Commit your changes (`git commit -m 'Add amazing feature'`)
   - Push to branch (`git push origin feature/amazing-feature`)
   - Open a Pull Request

### 📝 Version History

- **v1.3.0**: Current stable version with comprehensive features
- **v0.0.0.5**: Refactored all code, optimized function interfaces
- **v0.0.0.4**: Added population genetics analysis functions
- **v0.0.0.3**: Extended metabolomics analysis tools
- **v0.0.0.2**: Enhanced microbiome analysis
- **v0.0.0.1**: Initial release

### 📖 Citation

If you use bioRtools in your research, please cite:

```
Li, X. (2024). bioRtools: Convenience Functions for Biological Data Processing. 
R package version 1.3.0. https://github.com/lixiang117423/bioRtools
```

### 📄 License

This project is licensed under the MIT License. See [LICENSE.md](https://claude.ai/chat/LICENSE.md) for details.

### 👨‍💻 Author Information

**Xiang LI**

- 🔧 Project Maintainer
- 📧 Email: lixiang117423@gmail.com
- 🐙 GitHub: [@lixiang117423](https://github.com/lixiang117423)

### 🙏 Acknowledgments

Thanks to all developers and users who contributed to this project, and the following excellent R packages:

- **Bioconductor**: DESeq2, clusterProfiler, SummarizedExperiment
- **tidyverse**: ggplot2, dplyr, tidyr, purrr
- **Multivariate Analysis**: vegan, mixOmics, FactoMineR, factoextra
- **Visualization Enhancement**: ggsci, ggprism, ggrepel, ggtext

------

## 🔗 Links

- **📂 GitHub Repository**: https://github.com/lixiang117423/bioRtools
- **📚 Documentation**: https://lixiang117423.github.io/bioRtools/
- **🐛 Bug Reports**: https://github.com/lixiang117423/bioRtools/issues
- **🌐 Package Website**: https://lixiang117423.github.io/bioRtools/

------

<div align="center">

**🌟 Star this repository if you find it useful! 🌟**

Made with ❤️ for the bioinformatics community