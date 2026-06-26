# bioRtools 开发规范

> 版本: 2.0 | 日期: 2026-06-08
> 基于 [tidyverse 风格指南](https://style.tidyverse.org/)，结合项目实际约定

---

## 一、命名规范（核心规则）

### 1.1 通用规则：全部使用 snake_case

函数名、参数名、变量名、输出列名**一律使用下划线分隔的小写字母**。

```r
# Good
calc_beta_nti()
p_threshold
mean_expression
gene_count

# Bad
calcBetaNTI()       # camelCase
p.threshold         # dot.case
MeanExpression      # PascalCase
```

### 1.2 函数命名

- 动词开头，描述功能：`plot_*`, `calc_*`, `find_*`, `enrich_*`, `identify_*`
- 文件名与主函数名一致：`plot_volcano.R` → `plot_volcano()`

```r
# Good
plot_multi_volcano()
enrich_kegg()
find_degs_deseq2()

# Bad
PlotVolcano()       # PascalCase
plot.mv()           # 点号留给 S3 方法
cal_exp()           # 缩写不清晰
```

### 1.3 参数命名

参数名使用 snake_case。传递给外部包的命名参数保留该包的原始名称。

```r
# 我们的参数：snake_case
enrich_kegg(gene, kegg_db, p_adjust = 0.05, min_pathway_size = 3)

# 传递给外部包时，保留其参数名
picante::comdistnt(..., abundance.weighted = abundance_weighted)
clusterProfiler::enricher(..., qvalueCutoff = p_adjust)
```

### 1.4 输出列名

`mutate()`、`summarise()`、`rename()` 等创建的列名使用 snake_case。

```r
# Good
dplyr::mutate(mean_expression = mean(value), sd_expression = sd(value))

# Bad
dplyr::mutate(mean.expre = mean(value), SD = sd(value))
```

### 1.5 数据文件命名

`data/` 目录下的 `.rda` 文件遵循 `df.{analysis_type}.{data_type}.rda` 模式，全部小写。

```
df.pcoa.otu.rda
df.pcoa.sample.rda
df.rnaseq.gene.rda
df.rnaseq.plot.volcano.rda
df.call_dams_lefse.otu.rda
```

### 1.6 例外

以下情况保留原始命名，不做转换：

- **`na.rm`**：R 语言通用惯例，`mean()`, `sd()` 等基础函数均使用 `na.rm`
- **外部包的输出列名**：如 clusterProfiler 的 `p.adjust`、`GeneRatio` 等
- **向后兼容函数**：如 `CalExp2dCt()` 保留旧接口不变
- **S3 方法**：`print.my_class()` 使用点号

---

## 二、函数结构

### 2.1 模板

```r
#' One-line title
#'
#' Multi-paragraph description of what the function does.
#'
#' @param data A data frame containing ...
#' @param threshold Numeric threshold for ... (default: 0.05)
#'
#' @return A list containing:
#'   \item{table}{Data frame with results}
#'   \item{plot}{ggplot object}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- function_name(data = df, threshold = 0.05)
#' result$table
#' }
function_name <- function(data, threshold = 0.05) {
  # 1. Input validation
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }

  # 2. Core logic
  result <- data %>%
    dplyr::filter(.data$p_value < threshold)

  # 3. Return
  return(list(table = result))
}
```

### 2.2 内部辅助函数

使用 `@keywords internal` 标记，不导出。

```r
#' Calculate internal statistics
#' @keywords internal
calculate_stats <- function(x) {
  mean(x, na.rm = TRUE)
}
```

---

## 三、代码风格

### 3.1 管道操作

```r
# Good：每行一步
result <- data %>%
  dplyr::filter(p_value < 0.05) %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(mean_fc = mean(log2_fc))
```

### 3.2 命名空间

始终使用 `pkg::fun()` 显式调用，不依赖 `library()` 的副作用。

```r
# Good
dplyr::mutate()
ggplot2::ggplot()
stats::sd()

# Bad
mutate()   # 依赖 library(dplyr)
```

### 3.3 赋值与 return

- 赋值使用 `<-`，不用 `=`
- 函数末尾的最终结果自动返回，不需要 `return()`
- 仅在提前退出（错误处理、条件分支）时使用 `return()`

### 3.4 格式化

- 2 空格缩进，不用 Tab
- 运算符两侧加空格：`x + y`，`x == y`
- 逗号后加空格：`func(a, b, c)`
- `{` 放在行尾，`}` 独占一行

### 3.5 注释

注释解释 **为什么**（Why），不解释 **是什么**（What）。默认不写注释，代码本身应足够清晰。

```r
# Good
# Use 0.5 as threshold per published protocol
threshold <- 0.5

# Bad
# Set threshold to 0.5
threshold <- 0.5
```

---

## 四、文档规范（roxygen2）

### 4.1 导出函数必需标签

| 标签 | 用途 |
|------|------|
| `@title` | 一行标题 |
| `@description` | 功能描述 |
| `@param` | 每个参数的说明 |
| `@return` | 返回值及列名说明 |
| `@export` | 导出声明 |
| `@examples` | 可运行示例（用 `\dontrun{}` 包裹） |
| `@author` | 作者信息 |

### 4.2 示例规范

- 使用包内自带数据集：`data(df.xxx)`
- 用 `\dontrun{}` 包裹（避免 CRAN 检查问题）
- 展示最常用参数组合

---

## 五、图形函数规范

- 返回 ggplot 对象，**不直接打印**
- 图形参数使用描述性名称：`chr_colors` 而非 `colors`
- 使用 `theme_bio()` 或 `theme_bw()` 作为默认主题

```r
# Good：返回对象
p <- ggplot2::ggplot(data, ggplot2::aes(x, y)) +
  ggplot2::geom_point()
return(p)
```

---

## 六、错误处理

```r
# 参数验证：使用 stop()
if (!is.data.frame(data)) {
  stop("'data' must be a data frame")
}

# 可恢复问题：使用 warning()
if (any(is.na(data$value))) {
  warning("NA values detected, removing them")
}

# 外部操作：使用 tryCatch()
result <- tryCatch({
  risky_operation()
}, error = function(e) {
  warning("Operation failed: ", e$message)
  NULL
})
```

---

## 七、版本控制

- 遵循 [语义化版本](https://semver.org/)：`MAJOR.MINOR.PATCH`
- 每次修改更新 `DESCRIPTION` 版本号和 `CHANGELOG.md`
- Commit message 格式：`version X.Y.Z: brief description`
- 提交前运行 `devtools::document()` 确保文档同步

---

## 八、开发检查清单

### 新函数上线前

- [ ] 函数名 snake_case，动词开头
- [ ] 所有参数名 snake_case
- [ ] 输出列名 snake_case
- [ ] `@param` 每个参数有文档
- [ ] `@return` 说明返回结构
- [ ] `@examples` 可运行
- [ ] 参数验证完整
- [ ] 使用 `pkg::fun()` 显式命名空间
- [ ] `devtools::document()` 通过
- [ ] 包可正常 `pkgload::load_all()`

### 修改现有函数时

- [ ] 参数名改动需同步更新 roxygen `@param` 和 `@examples`
- [ ] 输出列名改动需同步更新下游函数引用
- [ ] 向后兼容函数保持旧接口不变

---

## 参考资源

- [tidyverse 风格指南](https://style.tidyverse.org/) — 本项目主要参考
- [R Packages (Hadley Wickham)](https://r-pkgs.org/) — R 包开发权威指南
- [roxygen2 文档](https://roxygen2.r-lib.org/) — 文档生成工具
