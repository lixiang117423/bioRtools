# bioRtools

[![R](https://img.shields.io/badge/R-%3E%3D2.10-blue)](https://www.r-project.org/)
[![Version](https://img.shields.io/badge/version-1.16.1-green)](https://github.com/lixiang117423/bioRtools)
[![License](https://img.shields.io/badge/license-MIT-yellow)](LICENSE.md)

`bioRtools` 是一个面向生物数据分析、统计和发表级可视化的 R 工具包。它整合了转录组、微生物组、代谢组、群体遗传、基因结构可视化、qPCR 分析以及常用绘图主题等工作流中的便利函数。

## 主要特性

- 支持转录组、微生物组、代谢组和群体遗传等多组学分析场景
- 多变量分析工作流：PCA、PCoA、RDA、sPLS-DA、OPLS-DA、PERMANOVA
- 差异分析：基于 DESeq2 的 DEG/DAM 分析和 LEfSe 标志物分析
- qPCR 分析：标准曲线、delta Ct、delta-delta Ct 和扩增效率校正
- 基因组可视化：曼哈顿图、QQ 图、LD 热图、共线性、motif、基因结构、PFAM 结构域和泛基因组稀释曲线
- 发表级 ggplot 主题、学术配色和 color/fill/colour 标尺

## 系统要求

- `DESCRIPTION` 声明的最低版本为 R >= 2.10
- 推荐使用 R >= 4.0.0，以获得更顺畅的依赖安装体验

## 安装

从 GitHub 安装开发版本：

```r
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

devtools::install_github("lixiang117423/bioRtools")
```

如果依赖安装失败，可先手动安装主要 Bioconductor 和 CRAN 依赖：

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

可选建议包包括 `chemhelper`、`ggpmisc`、`knitr`、`normentR`、`rmarkdown` 和 `testthat`。

## 功能模块

### 通用分析与统计

- `pca_analysis()` - 主成分分析
- `cor_analysis()` - 相关性分析及热图数据输出
- `lm_analysis()`、`get_lm_stats()`、`get_lm_stats_summary()`、`extract_lm_stats()`、`format_lm_stats()` - 线性模型分析辅助函数
- `anova_posthoc()` - 方差分析和事后检验
- `find_outliers()` - 异常值检测
- `label_signif()`、`label_significance()` - 显著性标签
- `reorder_heatmap()`、`gg_heatmap()`、`create_heatmap_trees()` - 热图辅助函数

### 微生物组分析

- `pcoa_analysis()` - 主坐标分析
- `rda_analysis()` - 冗余分析
- `permanova_test()` - PERMANOVA 检验
- `find_dams_deseq2()` - 基于 DESeq2 的差异丰度分析
- `find_dams_lefse()` - LEfSe 差异丰度分析
- `rarefy_table()` - 稀释抽样
- `top_taxa()` - 主要分类群提取
- `identify_core_microbiome()`、`fit_sloan_neutral_model()` - 核心微生物组和中性模型分析

### 转录组与富集分析

- `find_degs_deseq2()` - 基于 DESeq2 的差异表达分析
- `enrich_go()` - GO 富集分析
- `enrich_kegg()` - KEGG 富集分析
- `plot_volcano()`、`plot_multi_volcano()` - 火山图
- `run_wgcna_analysis()` - WGCNA 分析工作流

### qPCR 分析

- `calc_standard_curve()` - 标准曲线拟合
- `calc_expression_qpcr_efficiency()` - qPCR 扩增效率计算
- `calc_expression_delta_ct()` - delta Ct 分析
- `calc_expression_delta_delta_ct()` - delta-delta Ct 分析
- `calc_expression_standard_curve()` - 基于标准曲线的表达量分析

### 代谢组和多变量模型

- `opls_analysis()` - OPLS-DA
- `spls_analysis()` - sPLS-DA

### 群体遗传和基因组分析

- `manhattan_plot()`、`plot_manhattan()` - 曼哈顿图
- `plot_gwas_qq()` - GWAS QQ 图
- `plot_LDheatmap()` - LD 热图
- `ld_decay_threshold()` - LD 衰减阈值计算
- `admixture_phylo_analysis()`、`extract_tree_hierarchy()` - 群体结构和系统树辅助函数
- `pav_gwas()` - PAV-GWAS 辅助函数

### 基因、motif、结构域和共线性可视化

- `plot_gene_structure()`、`plot_gene_features()`、`plot_gene_features_labeled()` - 基因结构和特征图
- `plot_motif_location()`、`get_motif_from_meme()` - motif 解析和可视化
- `plot_pfam()`、`quick_pfam_plot()` - PFAM 结构域图
- `plot_synteny()` - 共线性图
- `plot_pangenome_rarefaction()` - 泛基因组稀释曲线
- `get_hap_from_heatmap()` - 从热图式数据中提取单倍型

### 数据转换和行统计工具

- `df_to_list()` - 数据框转列表
- `df2fasta()`、`fasta2df()` - FASTA 和数据框互转
- `get_methylkit_data()` - 提取 methylKit 数据
- `normalize_int()`、`scale01()`、`scale01_rows()`、`scale01_groups()`、`mutate_scale01()`、`mutate_scale01_named()` - 标准化辅助函数
- `row_mean()`、`row_sd()`、`row_cv()`、`row_min()`、`row_max()` - 行统计函数

### 主题与配色

- `theme_bio()`、`theme_prism()` - 绘图主题
- `pal_sci()`、`pal_nature()`、`pal_science()`、`pal_cell()`、`pal_jacs()`、`pal_fuel()`、`pal_chem_eng()`、`pal_nat_comm()`、`pal_shinkai()`、`pal_research()` - 学术配色
- `scale_color_*()`、`scale_colour_*()` 和 `scale_fill_*()` 系列提供离散和连续配色标尺

## 快速开始

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
  color.by = "species"
)

print(pca_result$plots$score_plot)
print(pca_result$eigenvalues)
```

### 相关性分析

```r
cor_result <- cor_analysis(
  data.1 = iris[, 1:2],
  data.2 = iris[, 3:4],
  method = "pearson"
)

print(cor_result$plot.cor)
```

### LEfSe 差异丰度分析

```r
data(df.call_DAMs_LEfSe.otu)
data(df.call_DAMs_LEfSe.sample)

lefse_result <- find_dams_lefse(
  data = df.call_DAMs_LEfSe.otu,
  sample = df.call_DAMs_LEfSe.sample,
  groupCol = "group",
  lda.threshold = 1.0
)

head(lefse_result)
```

### 火山图

```r
data(df.rnaseq.plot_volcano)

volcano_result <- plot_volcano(
  data = df.rnaseq.plot_volcano,
  x = "log2FoldChange",
  y = "padj",
  title = "Differential Expression"
)

print(volcano_result$plot.volcano)
print(volcano_result$data.summary)
```

## RNA-seq 工作流示例

```r
perform_rnaseq_analysis <- function(count_data, sample_data, go_db, kegg_db) {
  degs <- find_degs_deseq2(
    data = count_data,
    sample = sample_data,
    formula = ~group
  )

  sig_genes <- degs$gene[degs$padj < 0.05 & degs$regulation != "Not significant"]

  go_enrichment <- enrich_go(gene = sig_genes, go.db = go_db)
  kegg_enrichment <- enrich_kegg(gene = sig_genes, kegg.db = kegg_db)

  volcano_result <- plot_volcano(
    data = degs,
    x = "log2FoldChange",
    y = "padj"
  )

  list(
    degs = degs,
    go_enrichment = go_enrichment,
    kegg_enrichment = kegg_enrichment,
    volcano_plot = volcano_result$plot.volcano
  )
}
```

## 数据格式要求

### 表达、计数或丰度矩阵

- 行为基因、特征、分类群、OTU、ASV 或标记
- 列为样本
- 行名和列名应使用稳定的标识符
- DESeq2 工作流要求输入原始非负整数计数

### 样本信息表

- 行为样本
- 列为实验因子、分组变量、批次和协变量
- 样本标识符必须和输入矩阵中的样本名对应

### 富集注释表

- GO 注释表应包含 `gene`、`go_id` 和 `go_term`
- KEGG 注释表应包含 `gene`、`kegg_id` 和 `kegg_term`

## 故障排除

### 依赖安装失败

```r
update.packages(ask = FALSE)
install.packages("BiocManager")
BiocManager::install(ask = FALSE)
```

### 内存压力较大

```r
rm(list = ls())
gc()
```

在旧版 Windows R 中，也可以提高内存限制：

```r
memory.limit(size = 16000)
```

### 字体显示问题

```r
theme_bio(base_family = "Arial")
theme_prism(base_family = "Arial")
```

## 版本历史

完整发布记录见 [CHANGELOG.md](CHANGELOG.md)。当前包版本：`1.16.1`。

## 引用

如果你在研究中使用了 `bioRtools`，请引用：

```text
Li, X. (2026). bioRtools: Convenience Functions for Biological Data Processing.
R package version 1.16.1. https://github.com/lixiang117423/bioRtools
```

## 贡献

欢迎参与贡献：

1. 在 [GitHub Issues](https://github.com/lixiang117423/bioRtools/issues) 中报告问题
2. 通过 Issues 提出功能建议
3. Fork 仓库，创建功能分支，提交修改，推送分支并发起 Pull Request

## 许可证

本项目采用 MIT 许可证。详见 [LICENSE.md](LICENSE.md)。

## 作者

**Xiang LI**

- 项目维护者
- Email: lixiang117423@gmail.com
- GitHub: [@lixiang117423](https://github.com/lixiang117423)

## 链接

- GitHub 仓库：<https://github.com/lixiang117423/bioRtools>
- 文档网站：<https://lixiang117423.github.io/bioRtools/>
- 问题反馈：<https://github.com/lixiang117423/bioRtools/issues>
