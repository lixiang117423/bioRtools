# bioRtools

[English](#english) | [中文](#中文)

---

## 中文

### 简介

`bioRtools` 是一个专为生物数据处理设计的综合性R包，提供便利函数套件用于生物数据分析和可视化。该包涵盖多变量分析、统计检验、差异分析、数据预处理、基因组学可视化、微生物组学分析、代谢组学工作流程以及发表质量的绘图主题。bioRtools 标准化了跨组学数据类型的复杂分析流程，以最少的代码要求提供一致的、发表质量的可视化输出。

### 主要特性

- **🧬 多组学支持**: 转录组学、代谢组学、微生物组学、群体遗传学
- **📊 多变量分析**: PCA、PCoA、RDA、sPLS-DA、OPLS-DA
- **🔬 统计分析**: ANOVA、相关性分析、线性回归、PERMANOVA
- **📈 可视化**: 发表质量的图表，包括火山图、曼哈顿图等
- **⚡ 易于使用**: 简化的函数接口，标准化的输出格式
- **🎨 美观主题**: 内置生物数据可视化主题

### 系统要求

- R (>= 2.10)
- 建议使用 R >= 4.0.0 以获得最佳性能

### 安装

#### 从GitHub安装开发版本（推荐）

```r
# 安装必要的依赖包
if (!require(devtools)) install.packages("devtools")

# 安装bioRtools
devtools::install_github("lixiang117423/bioRtools")
```

#### 安装依赖包

bioRtools依赖多个Bioconductor和CRAN包。如果遇到依赖问题，请手动安装：

```r
# 安装Bioconductor包
if (!require(BiocManager)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "clusterProfiler", "lefser", "SummarizedExperiment"))

# 安装CRAN包
install.packages(c("ggplot2", "dplyr", "vegan", "mixOmics", "ropls", 
                   "factoextra", "FactoMineR", "ggsci", "ggprism"))
```

### 主要功能模块

#### 🔧 通用分析函数
- `pca_analysis()`: 主成分分析 (PCA)
- `cor_analysis()`: 相关性分析和可视化
- `lm_analysis()`: 线性回归分析和可视化
- `find_outliers()`: 异常值检测
- `theme_bio()`: 生物数据可视化主题
- `reorder_heatmap()`: 热图数据重排序
- `anova_posthoc()`: 方差分析及事后检验

#### 🧪 代谢组学分析
- `opls_analysis()`: 正交偏最小二乘判别分析 (OPLS-DA)
- `spls_analysis()`: 稀疏偏最小二乘判别分析 (sPLS-DA)

#### 🦠 微生物组学分析
- `pcoa_analysis()`: 主坐标分析 (PCoA)
- `rda_analysis()`: 冗余分析 (RDA)
- `find_dams_deseq2()`: DESeq2差异丰度分析
- `find_dams_lefse()`: LEfSe差异分析
- `permanova_test()`: PERMANOVA置换检验
- `rarefy_table()`: 稀释抽样
- `top_taxa()`: 获取主要分类群

#### 🧬 转录组学分析
- `find_degs_deseq2()`: DESeq2差异表达分析
- `enrich_go()`: GO功能富集分析
- `enrich_kegg()`: KEGG通路富集分析
- `plot_volcano()`: 火山图绘制

#### 🧮 群体遗传学分析
- `manhattan_plot()`: 曼哈顿图绘制
- `admixture_phylo_analysis()`: 群体结构和系统发育分析
- `plot_LDheatmap`: 绘制LDheatmap

#### 🛠️ 实用工具
- `df_to_list()`: 数据框转换为列表
- `plot_manhattan()`: 通用曼哈顿图绘制函数

### 快速开始

#### 主成分分析 (PCA)

```r
library(bioRtools)

# 准备数据
data <- iris[,1:4]
sample_info <- data.frame(
  sample = paste0("sample", 1:150),
  species = iris$Species
)

# 进行PCA分析
pca_result <- pca_analysis(data = data, sample = sample_info)

# 查看结果
print(pca_result$plot.pca)
print(pca_result$variance.explained)
```

#### 相关性分析

```r
# 两组数据的相关性分析
cor_result <- cor_analysis(
  data.1 = iris[,1:2], 
  data.2 = iris[,3:4],
  method = "pearson"
)

# 查看相关性热图
print(cor_result$plot.cor)
```

#### 线性回归分析

```r
# 线性回归分析并可视化
lm_result <- lm_analysis(
  data = iris,
  x = "Sepal.Length",
  y = "Sepal.Width", 
  group = "Species",
  color = "Species"
)

# 查看回归图
print(lm_result$plot.lm)
print(lm_result$result.lm)
```

#### 差异表达分析

```r
# 使用内置测试数据
data(df.rnaseq.count)
data(df.rnaseq.sample)

# DESeq2差异表达分析
degs_result <- find_degs_deseq2(
  count.table = df.rnaseq.count,
  sample.table = df.rnaseq.sample,
  design = ~ group,
  contrast = c("group", "treatment", "control")
)

# 查看结果
print(head(degs_result$degs))
```

#### GO富集分析

```r
# 使用内置测试数据
data(df.rnaseq.degs)
data(df.rnaseq.go)

# GO富集分析
go_result <- enrich_go(
  gene = df.rnaseq.degs$gene,
  go.db = df.rnaseq.go,
  p.adjust = 0.05
)

# 查看富集结果
print(head(go_result))
```

### 数据格式要求

#### 表达矩阵
- 行：基因/特征
- 列：样本
- 数值：原始计数或标准化表达量

#### 样本信息表
- 必须包含样本分组信息
- 样本名需与表达矩阵列名对应

#### 注释数据库
- GO数据库：包含gene、go.id、go.term列
- KEGG数据库：包含gene、kegg.id、kegg.term列

### 高级用法

#### 自定义可视化主题

```r
library(ggplot2)

# 使用bioRtools主题
p <- ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width)) +
  geom_point(aes(color = Species)) +
  theme_bio() +
  labs(title = "使用bioRtools主题的散点图")

print(p)
```

#### 批量分析流程

```r
# 完整的差异分析流程
perform_differential_analysis <- function(count_data, sample_data, go_db, kegg_db) {
  # 1. 差异表达分析
  degs <- find_degs_deseq2(
    count.table = count_data,
    sample.table = sample_data,
    design = ~ group,
    contrast = c("group", "treatment", "control")
  )
  
  # 2. GO富集分析
  go_enrichment <- enrich_go(
    gene = degs$degs$gene[degs$degs$padj < 0.05],
    go.db = go_db
  )
  
  # 3. KEGG富集分析
  kegg_enrichment <- enrich_kegg(
    gene = degs$degs$gene[degs$degs$padj < 0.05],
    kegg.db = kegg_db
  )
  
  return(list(
    degs = degs,
    go = go_enrichment,
    kegg = kegg_enrichment
  ))
}
```

### 故障排除

#### 常见问题

**问题1**: 安装依赖包失败
```r
# 解决方案：更新R版本或手动安装问题包
update.packages()
```

**问题2**: 内存不足
```r
# 解决方案：增加内存限制
memory.limit(size = 8000)  # Windows
```

**问题3**: 绘图中文显示问题
```r
# 解决方案：设置中文字体
theme_bio(base_family = "STSong")  # macOS
theme_bio(base_family = "SimSun")  # Windows
```

### 贡献指南

我们欢迎各种形式的贡献！

1. **报告Bug**: 在[GitHub Issues](https://github.com/lixiang117423/bioRtools/issues)中报告
2. **功能建议**: 通过Issues提出新功能建议
3. **代码贡献**: 
   - Fork项目
   - 创建功能分支
   - 提交Pull Request

### 版本历史

- **v0.0.0.5**: 重构所有代码，优化函数接口
- **v0.0.0.4**: 添加群体遗传学分析功能
- **v0.0.0.3**: 扩展代谢组学分析工具
- **v0.0.0.2**: 完善微生物组学分析
- **v0.0.0.1**: 初始版本发布

### 引用

如果您在研究中使用了bioRtools，请引用：

```
Li, X. (2024). bioRtools: Convenience Functions for Biological Data Processing. 
R package version 0.0.0.5. https://github.com/lixiang117423/bioRtools
```

### 许可证

本项目采用MIT许可证。详见[LICENSE.md](LICENSE.md)文件。

### 作者信息

**Xiang LI**
- 项目维护者
- Email: lixiang117423@gmail.com
- GitHub: [@lixiang117423](https://github.com/lixiang117423)

### 致谢

感谢所有为此项目做出贡献的开发者和用户，以及以下优秀的R包：
- DESeq2, clusterProfiler (Bioconductor)
- ggplot2, dplyr (tidyverse)
- vegan (生态学分析)
- mixOmics (多变量分析)

---

## English

### Introduction

bioRtools is a comprehensive R package designed for biological data processing, providing a suite of convenience functions for biological data analysis and visualization. It includes tools for multivariate analysis, statistical testing, differential analysis, data preprocessing, genomics visualization, microbiome analysis, metabolomics workflows, and publication-ready plotting themes. The package standardizes complex analytical workflows across omics data types and provides consistent, publication-quality visualization outputs with minimal code requirements.

### Key Features

- **🧬 Multi-omics Support**: Transcriptomics, metabolomics, microbiomics, population genetics
- **📊 Multivariate Analysis**: PCA, PCoA, RDA, sPLS-DA, OPLS-DA
- **🔬 Statistical Analysis**: ANOVA, correlation analysis, linear regression, PERMANOVA
- **📈 Visualization**: Publication-quality plots including volcano plots, Manhattan plots
- **⚡ Easy to Use**: Simplified function interfaces with standardized output formats
- **🎨 Beautiful Themes**: Built-in themes for biological data visualization

### System Requirements

- R (>= 2.10)
- Recommended: R >= 4.0.0 for optimal performance

### Installation

#### Install Development Version from GitHub (Recommended)

```r
# Install required dependencies
if (!require(devtools)) install.packages("devtools")

# Install bioRtools
devtools::install_github("lixiang117423/bioRtools")
```

#### Install Dependencies

bioRtools depends on several Bioconductor and CRAN packages. If you encounter dependency issues, install manually:

```r
# Install Bioconductor packages
if (!require(BiocManager)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "clusterProfiler", "lefser", "SummarizedExperiment"))

# Install CRAN packages
install.packages(c("ggplot2", "dplyr", "vegan", "mixOmics", "ropls", 
                   "factoextra", "FactoMineR", "ggsci", "ggprism"))
```

### Main Function Modules

#### 🔧 General Analysis Functions
- `pca_analysis()`: Principal Component Analysis (PCA)
- `cor_analysis()`: Correlation analysis and visualization
- `lm_analysis()`: Linear regression analysis and visualization
- `find_outliers()`: Outlier detection
- `theme_bio()`: Biological data visualization theme
- `reorder_heatmap()`: Heatmap data reordering
- `anova_posthoc()`: ANOVA with post-hoc tests

#### 🧪 Metabolomics Analysis
- `opls_analysis()`: Orthogonal Partial Least Squares Discriminant Analysis (OPLS-DA)
- `spls_analysis()`: Sparse Partial Least Squares Discriminant Analysis (sPLS-DA)

#### 🦠 Microbiome Analysis
- `pcoa_analysis()`: Principal Coordinate Analysis (PCoA)
- `rda_analysis()`: Redundancy Analysis (RDA)
- `find_dams_deseq2()`: DESeq2-based differential abundance analysis
- `find_dams_lefse()`: LEfSe differential analysis
- `permanova_test()`: PERMANOVA permutation test
- `rarefy_table()`: Rarefaction sampling
- `top_taxa()`: Extract top taxa

#### 🧬 Transcriptomics Analysis
- `find_degs_deseq2()`: DESeq2 differential expression analysis
- `enrich_go()`: GO enrichment analysis
- `enrich_kegg()`: KEGG pathway enrichment analysis
- `plot_volcano()`: Volcano plot generation

#### 🧮 Population Genetics Analysis
- `manhattan_plot()`: Manhattan plot generation
- `admixture_phylo_analysis()`: Population structure and phylogenetic analysis

#### 🛠️ Utility Functions
- `df_to_list()`: Convert data frame to list
- `plot_manhattan()`: Generic Manhattan plot function

### Quick Start

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
pca_result <- pca_analysis(data = data, sample = sample_info)

# View results
print(pca_result$plot.pca)
print(pca_result$variance.explained)
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

#### Linear Regression Analysis

```r
# Linear regression analysis with visualization
lm_result <- lm_analysis(
  data = iris,
  x = "Sepal.Length",
  y = "Sepal.Width", 
  group = "Species",
  color = "Species"
)

# View regression plot
print(lm_result$plot.lm)
print(lm_result$result.lm)
```

#### Differential Expression Analysis

```r
# Using built-in test data
data(df.rnaseq.count)
data(df.rnaseq.sample)

# DESeq2 differential expression analysis
degs_result <- find_degs_deseq2(
  count.table = df.rnaseq.count,
  sample.table = df.rnaseq.sample,
  design = ~ group,
  contrast = c("group", "treatment", "control")
)

# View results
print(head(degs_result$degs))
```

#### GO Enrichment Analysis

```r
# Using built-in test data
data(df.rnaseq.degs)
data(df.rnaseq.go)

# GO enrichment analysis
go_result <- enrich_go(
  gene = df.rnaseq.degs$gene,
  go.db = df.rnaseq.go,
  p.adjust = 0.05
)

# View enrichment results
print(head(go_result))
```

### Data Format Requirements

#### Expression Matrix
- Rows: genes/features
- Columns: samples
- Values: raw counts or normalized expression values

#### Sample Information Table
- Must contain sample grouping information
- Sample names should match expression matrix column names

#### Annotation Databases
- GO database: contains gene, go.id, go.term columns
- KEGG database: contains gene, kegg.id, kegg.term columns

### Advanced Usage

#### Custom Visualization Themes

```r
library(ggplot2)

# Using bioRtools theme
p <- ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width)) +
  geom_point(aes(color = Species)) +
  theme_bio() +
  labs(title = "Scatter plot with bioRtools theme")

print(p)
```

#### Batch Analysis Pipeline

```r
# Complete differential analysis pipeline
perform_differential_analysis <- function(count_data, sample_data, go_db, kegg_db) {
  # 1. Differential expression analysis
  degs <- find_degs_deseq2(
    count.table = count_data,
    sample.table = sample_data,
    design = ~ group,
    contrast = c("group", "treatment", "control")
  )
  
  # 2. GO enrichment analysis
  go_enrichment <- enrich_go(
    gene = degs$degs$gene[degs$degs$padj < 0.05],
    go.db = go_db
  )
  
  # 3. KEGG enrichment analysis
  kegg_enrichment <- enrich_kegg(
    gene = degs$degs$gene[degs$degs$padj < 0.05],
    kegg.db = kegg_db
  )
  
  return(list(
    degs = degs,
    go = go_enrichment,
    kegg = kegg_enrichment
  ))
}
```

### Troubleshooting

#### Common Issues

**Issue 1**: Dependency installation failure
```r
# Solution: Update R version or manually install problematic packages
update.packages()
```

**Issue 2**: Insufficient memory
```r
# Solution: Increase memory limit
memory.limit(size = 8000)  # Windows
```

**Issue 3**: Font display issues in plots
```r
# Solution: Set appropriate font family
theme_bio(base_family = "Arial")  # Universal
```

### Contributing

We welcome contributions of all kinds!

1. **Report Bugs**: Report in [GitHub Issues](https://github.com/lixiang117423/bioRtools/issues)
2. **Feature Requests**: Suggest new features through Issues
3. **Code Contributions**: 
   - Fork the repository
   - Create a feature branch
   - Submit a Pull Request

### Version History

- **v0.0.0.5**: Refactored all code, optimized function interfaces
- **v0.0.0.4**: Added population genetics analysis functions
- **v0.0.0.3**: Extended metabolomics analysis tools
- **v0.0.0.2**: Enhanced microbiome analysis
- **v0.0.0.1**: Initial release

### Citation

If you use bioRtools in your research, please cite:

```
Li, X. (2024). bioRtools: Convenience Functions for Biological Data Processing. 
R package version 0.0.0.5. https://github.com/lixiang117423/bioRtools
```

### License

This project is licensed under the MIT License. See [LICENSE.md](LICENSE.md) for details.

### Author Information

**Xiang LI**
- Project Maintainer
- Email: lixiang117423@gmail.com
- GitHub: [@lixiang117423](https://github.com/lixiang117423)

### Acknowledgments

Thanks to all developers and users who contributed to this project, and the following excellent R packages:
- DESeq2, clusterProfiler (Bioconductor)
- ggplot2, dplyr (tidyverse)
- vegan (ecological analysis)
- mixOmics (multivariate analysis)

---

## Links

- **GitHub Repository**: https://github.com/lixiang117423/bioRtools
- **Documentation**: https://lixiang117423.github.io/bioRtools/
- **Bug Reports**: https://github.com/lixiang117423/bioRtools/issues
- **Package Website**: https://lixiang117423.github.io/bioRtools/