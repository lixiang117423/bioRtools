# bioRtools R代码开发规范文档

## 版本: 1.0
## 日期: 2026-02-09
## 用途: 统一R包的代码结构、命名规范、文档格式

---

## 一、项目结构规范

### 1.1 标准R包目录结构

```
bioRtools/
├── DESCRIPTION          # 包描述文件
├── NAMESPACE            # 命名空间声明
├── LICENSE              # 许可证
├── README.md            # 包说明文档
├── R/                   # R函数源代码目录
│   ├── package-name.R   # 主函数文件
│   └── helper-*.R       # 辅助函数（使用@noRd标记）
├── man/                 # 函数文档目录（自动生成）
├── data/                # 数据文件目录
├── tests/               # 测试文件目录
├── vignettes/           # 教程文档目录
└── inst/                # 额外文件目录
```

### 1.2 文件命名规范

- 使用下划线分隔：`plot_multi_volcano.R`
- 文件名与主函数名保持一致
- 辅助函数使用`@noRd`标记，不导出

---

## 二、函数开发规范

### 2.1 函数定义结构

#### 完整函数模板

```r
#' 函数标题|Function Title
#'
#' 函数描述|Function description
#'
#' @param data 数据框|Data frame containing...
#' @param var1 参数1描述|Parameter 1 description
#' @param var2 参数2描述|Parameter 2 description (default: value)
#'
#' @return 返回值描述|Returns description
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(bioRtools)
#' library(dplyr)
#'
#' # 基本用法示例
#' result <- function_name(
#'   data = df,
#'   var1 = "value1"
#' )
#'
#' print(result)
#' }
function_name <- function(data, var1, var2 = "default") {

  # 输入验证
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }

  if (!var1 %in% names(data)) {
    stop(paste("Column '", var1, "' not found in data"))
  }

  # 核心逻辑
  result <- data %>%
    dplyr::mutate(
      new_col = dplyr::if_else(
        condition,
        true_value,
        false_value
      )
    )

  return(result)
}
```

### 2.2 参数验证规范

```r
# 必需参数验证
if (missing(data)) {
  stop("'data' argument is required")
}

# 数据类型验证
if (!is.data.frame(data)) {
  stop("'data' must be a data frame")
}

# 列存在性验证
if (!all(c("col1", "col2") %in% names(data))) {
  missing_cols <- setdiff(c("col1", "col2"), names(data))
  stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))
}

# 参数范围验证
if (threshold <= 0 || threshold >= 1) {
  stop("'threshold' must be between 0 and 1")
}
```

### 2.3 错误处理规范

```r
# 使用tryCatch处理可能出错的操作
result <- tryCatch({
  risky_operation()
}, error = function(e) {
  warning(paste("Operation failed:", e$message))
  return(NULL)
})

# 提供有意义的错误信息
if (nrow(data) == 0) {
  stop("Data frame is empty. Please provide a non-empty data frame")
}

# 使用warning而不是stop处理可恢复的错误
if (any(is.na(data$column))) {
  warning("Column contains NA values, these will be removed")
  data <- tidyr::drop_na(data, column)
}
```

### 2.4 内部函数规范

```r
#' Internal helper function
#' @noRd
helper_function <- function(x) {
  # 辅助逻辑
  result <- x + 1
  return(result)
}
```

---

## 三、文档规范（roxygen2）

### 3.1 必需文档元素

每个导出的函数必须包含：

- `@title` - 函数标题
- `@description` - 函数描述
- `@param` - 每个参数的详细说明
- `@return` - 返回值说明
- `@export` - 导出声明
- `@examples` - 可运行示例

### 3.2 参数文档格式

```r
#' @param data 输入数据框，必须包含以下列：
#'   \itemize{
#'     \item \code{cluster}: 分组信息
#'     \item \code{gene}: 基因名称
#'     \item \code{avg_log2FC}: log2倍数变化值
#'     \item \code{p_val_adj}: 调整后p值
#'   }
#' @param fc_threshold Fold change阈值，用于筛选差异表达基因（default: 1）
#' @param label_n 每组标注的基因数量（default: 5）
```

### 3.3 返回值文档格式

```r
#' @return A ggplot object with the following layers:
#'   \itemize{
#'     \item Scatter points for each gene
#'     \item Background strips for each cluster
#'     \item Gene labels for significant genes
#'     \item Color-coded up/down regulation
#'   }
```

### 3.4 Examples规范

```r
#' @examples
#' \dontrun{
#' library(bioRtools)
#' library(dplyr)
#'
#' # 基本用法
#' df <- read_tsv("data.tsv")
#' p <- plot_multi_volcano(
#'   df,
#'   fc_column = "avg_log2FC",
#'   cluster_column = "cluster",
#'   gene_column = "gene"
#' )
#'
#' # 保存图形
#' ggplot2::ggsave("volcano.png", p, width = 10, height = 6)
#' }
```

**注意：**
- 使用`\dontrun{}`避免CRAN检查时运行示例
- 示例必须可运行
- 包含必要的library调用
- 展示常用参数组合

---

## 四、命名规范

### 4.1 函数命名

- ✅ 使用动词开头：`plot_*`, `calculate_*`, `format_*`
- ✅ 使用下划线分隔：`plot_multi_volcano`, `anova_posthoc`
- ✅ 名称描述功能：`identify_core_microbiome`
- ❌ 避免使用点号（保留给S3方法）
- ❌ 避免缩写：`plot_mv`（不够清晰）

```r
# Good
plot_multi_volcano()
calculate_qv()
identify_core_microbiome()

# Bad
plotMultiVolcano()  # 驼峰命名
plot_mv()           # 缩写不清晰
plot.multi.volcano() # 点号留给S3
```

### 4.2 变量命名

```r
# Good
gene_names
p_value
log2_fc_threshold
cluster_colors

# Bad
geneNames       # 驼峰命名
p_value         # 不够描述性
x, y, tmp       # 无意义名称（除循环外）
```

### 4.3 常量命名

```r
# Good
DEFAULT_THREADS <- 12
MAX_CLUSTER_NUM <- 20
DEFAULT_PALETTE <- c("#3B9AB2", "#78B7C5")

# 位置：通常放在R/包名-global.R文件中
```

---

## 五、代码风格规范

### 5.1 缩进和空格

**遵循tidyverse风格：**
- 使用2个空格缩进（不使用tab）
- 赋值操作符周围加空格：`x <- 1`
- 逗号后加空格：`func(a, b, c)`
- 中缀运算符周围加空格：`x + y`, `x == y`

```r
# Good
result <- data %>%
  dplyr::filter(column > threshold) %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(mean_val = mean(value))

# Bad
result<-data%>%filter(column>threshold)%>%group_by(cluster)
```

### 5.2 花括号规范

```r
# Good - { 是行尾最后一个字符
if (condition) {
  # 代码缩进2个空格
  result
}

# Good - else 与前一个 } 同行
if (condition) {
  result_a
} else {
  result_b
}

# Bad
if (condition)
{
  result  # 缩进错误
}
```

### 5.3 函数调用规范

**长函数调用换行：**

```r
# Good - 单缩进风格
plot_multi_volcano(
  data = df,
  fc_column = "avg_log2FC",
  cluster_column = "cluster",
  gene_column = "gene",
  fc_threshold = 1
)

# Good - 悬挂缩进风格
plot_multi_volcano(data = df,
                   fc_column = "avg_log2FC",
                   cluster_column = "cluster",
                   gene_column = "gene")

# Bad
plot_multi_volcano(data = df, fc_column = "avg_log2FC",
  cluster_column = "cluster", gene_column = "gene")
```

### 5.4 管道操作规范

```r
# Good - 管道操作每行一个
result <- data %>%
  dplyr::filter(p_val < 0.05) %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(n = dplyr::n())

# Good - 匿名函数使用lambda语法（简短）
data %>%
  dplyr::mutate(new_col = purrr::map_dbl(values, \(x) mean(x^2)))

# Good - 匿名函数使用传统语法（多行）
data %>%
  dplyr::mutate(new_col = purrr::map_dbl(values, function(x) {
    y <- x^2
    mean(y, na.rm = TRUE)
  }))
```

### 5.5 命名空间规范

**始终使用明确的命名空间：**

```r
# Good
dplyr::filter()
ggplot2::aes()
tidyr::drop_na()
purrr::map()

# Bad - 避免使用library()后直接调用
filter()
aes()
```

### 5.6 return()使用规范

```r
# Good - 只在早期返回时使用return()
if (error) {
  return(NULL)
}
x + y  # 最后表达式自动返回

# Bad - 不必要的return()
add <- function(x, y) {
  return(x + y)
}
```

### 5.7 注释规范

```r
# Good - 解释"为什么"
# 使用0.5作为阈值是因为这是文献中常用的cut-off
threshold <- 0.5

# Bad - 解释"什么"
# 设置阈值为0.5
threshold <- 0.5

# Good - 注释符号后加空格
# Default to 12 threads for parallel processing
DEFAULT_THREADS <- 12

# Bad
#Default to 12 threads
DEFAULT_THREADS <- 12
```

---

## 六、图形函数规范

### 6.1 图形函数返回值

```r
# 图形函数应返回ggplot对象，不打印
plot_multi_volcano <- function(data) {
  p <- ggplot2::ggplot(data, ggplot2::aes(x, y)) +
    ggplot2::geom_point()

  return(p)  # 返回对象，让用户决定是否打印
}
```

### 6.2 图形参数命名

```r
# Good - 描述性参数名
cluster_colors  # 而不是 colors
legend_position # 而不是 pos
y_limits        # 而不是 ylim
```

### 6.3 主题参数

```r
# 提供主题参数，而不是硬编码
plot_multi_volcano <- function(
  data,
  theme = ggplot2::theme_bw(),  # 允许自定义主题
  base_size = 12
) {
  # ...
}
```

---

## 七、测试规范

### 7.1 测试文件结构

```r
# tests/testthat/test-plot-multi-volcano.R
test_that("plot_multi_volcano works", {
  # 准备测试数据
  test_data <- tibble::tibble(
    cluster = c("A", "A", "B", "B"),
    gene = c("g1", "g2", "g3", "g4"),
    avg_log2FC = c(1, -1, 2, -2),
    p_val_adj = c(0.01, 0.02, 0.001, 0.03)
  )

  # 测试基本功能
  p <- plot_multi_volcano(
    test_data,
    fc_column = "avg_log2FC",
    cluster_column = "cluster",
    gene_column = "gene"
  )

  # 验证返回值
  expect_s3_class(p, "ggplot")
})

test_that("plot_multi_volcano validates input", {
  # 测试错误处理
  expect_error(
    plot_multi_volcano(NULL),
    "must be a data frame"
  )
})
```

### 7.2 测试覆盖率目标

- 单元测试覆盖率 > 80%
- 所有导出函数必须有测试
- 关键逻辑必须有测试
- 错误处理必须有测试

---

## 八、代码检查清单

### 8.1 函数开发检查

- [ ] 函数命名使用snake_case，动词开头
- [ ] 所有参数有`@param`文档
- [ ] 有清晰的`@return`说明
- [ ] 包含`@export`（导出函数）或`@noRd`（内部函数）
- [ ] 包含可运行的`@examples`
- [ ] 参数验证完整
- [ ] 错误信息清晰有意义
- [ ] 使用明确的命名空间（`dplyr::mutate`）
- [ ] 只在早期返回使用`return()`
- [ ] 注释解释"为什么"而非"什么"

### 8.2 代码风格检查

- [ ] 2个空格缩进
- [ ] 赋值操作符使用`<-`而非`=`
- [ ] 逗号后有空格
- [ ] 花括号位置正确
- [ ] 函数调用格式化合理
- [ ] 变量名使用snake_case
- [ ] 注释符号后有空格

### 8.3 文档检查

- [ ] 标题清晰描述功能
- [ ] 描述说明函数用途
- [ ] 参数说明完整（包含类型、默认值）
- [ ] 返回值说明详细
- [ ] 示例可运行
- [ ] 示例展示常用用法

---

## 九、常见问题

### Q1: 何时使用`::`调用函数？
**A:** 总是使用显式命名空间调用函数（如`dplyr::filter()`），除了：
- 基础R函数：`sum()`, `mean()`, `list()`等
- 包内定义的辅助函数

### Q2: 何时使用`return()`？
**A:**
- ✅ 早期返回（错误处理、条件分支）
- ✅ 函数中间返回
- ❌ 函数最后一行返回（让R自动返回）

### Q3: 如何处理缺失值？
**A:**
```r
# Good - 显式处理
data %>%
  dplyr::filter(!is.na(column))

# Good - 使用drop_na
data %>%
  tidyr::drop_na(column)

# Bad - 隐式处理（可能导致意外结果）
sum(data$column)  # 可能返回NA
```

### Q4: 如何选择默认参数值？
**A:**
- 选择最常用的值作为默认
- 在文档中说明默认值的理由
- 考虑使用NULL表示"不设置"而非某个特定值

### Q5: 如何组织长函数？
**A:**
- 将复杂逻辑提取为内部函数（`@noRd`）
- 使用注释分隔主要步骤
- 每个步骤保持在可理解的长度

---

## 十、版本历史

| 版本 | 日期 | 主要变更 |
|------|------|----------|
| 1.0 | 2026-02-09 | - 初始版本<br>- 基于tidyverse风格指南<br>- 参考bioRtools现有代码 |

---

## 十一、参考资源

- [tidyverse风格指南](https://style.tidyverse.org/)
- [R Packages (Hadley Wickham)](https://r-pkgs.org/)
- [roxygen2文档](https://roxygen2.r-lib.org/)
- [testthat测试](https://testthat.r-lib.org/)

**重要：所有代码修改必须使用git进行版本控制！！！**
