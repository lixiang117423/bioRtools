# bioRtools 🧬

[![R](https://img.shields.io/badge/R-%3E%3D2.10-blue)](https://www.r-project.org/) [![Version](https://img.shields.io/badge/version-1.3.0-green)](https://github.com/lixiang117423/bioRtools) [![License](https://img.shields.io/badge/license-MIT-yellow)](https://claude.ai/chat/LICENSE.md)

[English](./README_en.md) | [中文](./README.md)

------

## 中文 🇨🇳

### 📋 简介

`bioRtools` 是一个专为生物数据处理设计的综合性R包，提供便利函数套件用于生物数据分析和可视化。该包标准化了跨组学数据类型的复杂分析流程，以最少的代码要求提供一致的、发表质量的可视化输出。

### ✨ 主要特性

- **🧬 多组学支持**: 转录组学、代谢组学、微生物组学、群体遗传学
- **📊 多变量分析**: PCA、PCoA、RDA、sPLS-DA、OPLS-DA
- **🔬 统计分析**: ANOVA、相关性分析、线性回归、PERMANOVA
- **📈 高质量可视化**: 发表级别的图表，包括火山图、曼哈顿图等
- **⚡ 易于使用**: 简化的函数接口，标准化的输出格式
- **🎨 美观主题**: 内置生物数据可视化主题

### 💻 系统要求

- **R版本**: R (≥ 2.10)
- **推荐版本**: R ≥ 4.0.0 以获得最佳性能

### 📦 安装

#### 从GitHub安装开发版本（推荐）

```r
# 安装必要的依赖包
if (!require(devtools)) install.packages("devtools")

# 安装bioRtools
devtools::install_github("lixiang117423/bioRtools")
```

#### 手动安装依赖包

如果遇到依赖问题，请手动安装：

```r
# 安装Bioconductor包
if (!require(BiocManager)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "clusterProfiler", "lefser", "SummarizedExperiment"))

# 安装CRAN包
install.packages(c("ggplot2", "dplyr", "vegan", "mixOmics", "ropls", 
                   "factoextra", "FactoMineR", "ggsci", "ggprism", "ggrepel",
                   "rstatix", "broom", "scales", "tidyr"))
```

### 🛠️ 主要功能模块

#### 🔧 通用分析函数

- `pca_analysis()` - 主成分分析 (PCA)
- `cor_analysis()` - 相关性分析和可视化
- `lm_analysis()` - 线性回归分析和可视化
- `find_outliers()` - 异常值检测
- `theme_bio()` - 生物数据可视化主题
- `reorder_heatmap()` - 热图数据重排序
- `anova_posthoc()` - 方差分析及事后检验

#### 🧪 代谢组学分析

- `opls_analysis()` - 正交偏最小二乘判别分析 (OPLS-DA)
- `spls_analysis()` - 稀疏偏最小二乘判别分析 (sPLS-DA)

#### 🦠 微生物组学分析

- `pcoa_analysis()` - 主坐标分析 (PCoA)
- `rda_analysis()` - 冗余分析 (RDA)
- `find_dams_deseq2()` - DESeq2差异丰度分析
- `find_dams_lefse()` - LEfSe差异分析
- `permanova_test()` - PERMANOVA置换检验
- `rarefy_table()` - 稀释抽样
- `top_taxa()` - 获取主要分类群

#### 🧬 转录组学分析

- `find_degs_deseq2()` - DESeq2差异表达分析
- `enrich_go()` - GO功能富集分析
- `enrich_kegg()` - KEGG通路富集分析
- `plot_volcano()` - 火山图绘制

#### 🧮 群体遗传学分析

- `manhattan_plot()` - 曼哈顿图绘制
- `admixture_phylo_analysis()` - 群体结构和系统发育分析
- `plot_LDheatmap()` - 绘制LD热图

#### 🔧 实用工具函数

- `df_to_list()` - 数据框转换为列表
- `plot_manhattan()` - 通用曼哈顿图绘制
- `scale01()` - 数据标准化函数
- `row_mean()`, `row_sd()` - 行统计函数

### 🚀 快速开始

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
pca_result <- pca_analysis(
  data = data, 
  sample = sample_info,
  color.by = "species"
)

# 查看结果
print(pca_result$plots$score_plot)
print(pca_result$eigenvalues)
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

#### LEfSe差异分析

```r
# 微生物组差异分析示例
lefse_result <- find_dams_lefse(
  data = abundance_matrix,           # 特征丰度矩阵
  sample = sample_metadata,          # 样本信息
  groupCol = "treatment",            # 分组列名
  lda.threshold = 2.0               # LDA阈值
)

# 查看显著差异特征
print(head(lefse_result))
```

#### 火山图绘制

```r
# 绘制差异表达基因火山图
volcano_plot <- plot_volcano(
  data = deg_results,               # 差异分析结果
  x = "log2FoldChange", 
  y = "padj",
  label = "gene"
)

print(volcano_plot)
```

### 📊 数据格式要求

#### 表达/丰度矩阵

- **行**: 基因/特征/OTU
- **列**: 样本
- **数值**: 原始计数、标准化表达量或相对丰度

#### 样本信息表

- **行**: 样本
- **列**: 实验因子和协变量
- **要求**: 样本名需与表达矩阵列名对应

#### 注释数据库格式

- **GO数据库**: 包含 `gene`, `go.id`, `go.term` 列
- **KEGG数据库**: 包含 `gene`, `kegg.id`, `kegg.term` 列

### 🎨 高级用法

#### 自定义可视化主题

```r
library(ggplot2)

# 使用bioRtools主题
p <- ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width)) +
  geom_point(aes(color = Species), size = 2) +
  theme_bio(base_size = 12) +
  labs(title = "使用bioRtools主题的散点图",
       x = "萼片长度", y = "萼片宽度")

print(p)
```

#### 批量分析工作流

```r
# 完整的转录组差异分析流程
perform_rnaseq_analysis <- function(count_data, sample_data, go_db, kegg_db) {
  
  # 1. 差异表达分析
  cat("🔬 执行差异表达分析...\n")
  degs <- find_degs_deseq2(
    count.table = count_data,
    sample.table = sample_data,
    design = ~ condition,
    contrast = c("condition", "treatment", "control")
  )
  
  # 2. GO富集分析
  cat("📊 执行GO富集分析...\n")
  sig_genes <- degs$degs$gene[degs$degs$padj < 0.05]
  go_enrichment <- enrich_go(gene = sig_genes, go.db = go_db)
  
  # 3. KEGG富集分析
  cat("🧬 执行KEGG富集分析...\n")
  kegg_enrichment <- enrich_kegg(gene = sig_genes, kegg.db = kegg_db)
  
  # 4. 生成火山图
  cat("📈 生成可视化图表...\n")
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

### ⚠️ 故障排除

#### 常见问题解决方案

**问题1**: 安装依赖包失败

```r
# 解决方案：更新R和包管理器
update.packages(ask = FALSE)
install.packages("BiocManager")
```

**问题2**: 内存不足错误

```r
# 解决方案：增加内存限制（Windows）
memory.limit(size = 16000)
# 或者清理工作环境
rm(list = ls())
gc()
```

**问题3**: 中文字体显示问题

```r
# macOS系统
theme_bio(base_family = "STSong")
# Windows系统  
theme_bio(base_family = "SimSun")
# Ubuntu系统
theme_bio(base_family = "WenQuanYi Micro Hei")
```

### 🤝 贡献指南

我们欢迎各种形式的贡献！

1. **🐛 报告Bug**: 在[GitHub Issues](https://github.com/lixiang117423/bioRtools/issues)中报告

2. **💡 功能建议**: 通过Issues提出新功能建议

3. 💻 代码贡献

   :

   - Fork本项目
   - 创建功能分支 (`git checkout -b feature/amazing-feature`)
   - 提交更改 (`git commit -m 'Add amazing feature'`)
   - 推送到分支 (`git push origin feature/amazing-feature`)
   - 开启Pull Request

### 📝 版本历史

- **v1.3.0**: 当前稳定版本，功能完善
- **v0.0.0.5**: 重构所有代码，优化函数接口
- **v0.0.0.4**: 添加群体遗传学分析功能
- **v0.0.0.3**: 扩展代谢组学分析工具
- **v0.0.0.2**: 完善微生物组学分析
- **v0.0.0.1**: 初始版本发布

### 📖 引用

如果您在研究中使用了bioRtools，请引用：

```
Li, X. (2024). bioRtools: Convenience Functions for Biological Data Processing. 
R package version 1.3.0. https://github.com/lixiang117423/bioRtools
```

### 📄 许可证

本项目采用MIT许可证。详见[LICENSE.md](https://claude.ai/chat/LICENSE.md)文件。

### 👨‍💻 作者信息

**Xiang LI**

- 🔧 项目维护者
- 📧 Email: lixiang117423@gmail.com
- 🐙 GitHub: [@lixiang117423](https://github.com/lixiang117423)

### 🙏 致谢

感谢所有为此项目做出贡献的开发者和用户，以及以下优秀的R包：

- **Bioconductor**: DESeq2, clusterProfiler, SummarizedExperiment
- **tidyverse**: ggplot2, dplyr, tidyr, purrr
- **多变量分析**: vegan, mixOmics, FactoMineR, factoextra
- **可视化增强**: ggsci, ggprism, ggrepel, ggtext

------

</div>