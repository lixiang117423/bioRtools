# bioRtools 开发规范

> 版本: 2.2 | 日期: 2026-08-14
> 基于 [tidyverse 风格指南](https://style.tidyverse.org/)，结合项目实际约定
>
> v2.2 变更：修复模板与规则矛盾；新增测试规范章节；明确依赖声明硬规则；
> 扩展全局状态管理清单（RNG/设备/集群）；补充消息机制与开发中函数约定。

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
- 例外：同族函数（`scale_color_*`、`row_*` 等变体家族）可共享一个文件，
  文件以家族主题命名（如 `academic_color_scales.R`）

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
#' One-line title（首行自由文本即标题，不写 @title 标签）
#'
#' Multi-paragraph description of what the function does.
#' （第二段自由文本即 @description，不写 @description 标签）
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

  # 3. 末尾直接返回，不写 return()
  list(table = result)
}
```

### 2.2 内部辅助函数

使用 `@keywords internal` 标记，不导出。**默认同时加 `@noRd`**，
避免内部实现细节污染参考索引（man/ 目录）。

```r
#' Calculate internal statistics
#' @keywords internal
#' @noRd
calculate_stats <- function(x) {
  mean(x, na.rm = TRUE)
}
```

### 2.3 导出 vs 内部函数判定

- **导出（加 `@export`）**：用户直接调用、属于核心功能、API 稳定不会频繁变动。
- **保持内部（`@keywords internal`，不加 `@export`）**：实现细节、可能重构、仅包内复用。内部函数可假设输入已被上游校验，只检查关键不变量。
- 可选约定：内部辅助函数用 `.` 前缀命名（如 `.compute_metrics()`）与导出函数视觉区分。bioRtools 目前未统一使用，**新代码可选，不强制**。

### 2.4 开发中函数（WIP）

尚未完成、暂不发布的功能遵循以下约定：

- 不写 `@export`，加 `@keywords internal`（可保留 Rd 页面供本地查阅，暂不加 `@noRd`）
- **包内引用禁止使用 `bioRtools::fn()` 自引用**：未导出函数经 `::` 调用在安装版会直接报错（`devtools::load_all()` 下一切正常，极易漏检）。包内一律直接调用 `fn()`
- 已导出的函数不得依赖未导出的 WIP 函数（如需过渡，用直接调用而非 `::`）
- 提交前确认 WIP 函数不在 NAMESPACE 的 export 列表中

> 历史教训：`plot_volcano()` 曾以 `bioRtools::theme_bio()` 引用未导出的
> `theme_bio`，导致安装版自 v1.20.0 起运行时报错。

---

## 三、代码风格

### 3.1 管道操作

统一用 magrittr 的 `%>%`，不混用原生 `|>`（bioRtools 全库一致使用 `%>%`）。

```r
# Good：每行一步
result <- data %>%
  dplyr::filter(p_value < 0.05) %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(mean_fc = mean(log2_fc))
```

### 3.2 命名空间

始终使用 `pkg::fun()` 显式调用，不依赖 `library()` 的副作用。
**例外**：调用本包未导出的内部函数时直接调用（见 2.4）。

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
- 函数末尾的最终结果自动返回，**不写 `return()`**
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

### 3.6 路径与包内数据

严禁硬编码绝对路径（尤其含用户名的路径）。访问包内置数据用 `system.file()`：

```r
# Good — 包内 extdata 数据
path <- system.file("extdata", "example.otu.csv", package = "bioRtools")

# Bad — 硬编码本地路径
read.csv("/Users/xxx/data/example.csv")
```

用户传入的路径参数不要假设固定位置，交给调用方决定。

### 3.7 循环边界

循环一律用 `seq_len(nrow(x))` / `seq_along(x)`，不用 `1:nrow(x)`：
空数据框时 `1:0` 会静默产生 `c(1, 0)`。

```r
# Good
for (i in seq_len(nrow(df))) { ... }

# Bad
for (i in 1:nrow(df)) { ... }   # nrow(df) == 0 时出错
```

---

## 四、文档规范（roxygen2）

### 4.1 导出函数必需标签

| 标签 | 用途 |
|------|------|
| 首行自由文本 | 一行标题（**不写 `@title` 标签**） |
| 第二段自由文本 | 功能描述（**不写 `@description` 标签**） |
| `@param` | 每个参数的说明 |
| `@return` | 返回值及列名说明 |
| `@export` | 导出声明 |
| `@examples` | 可运行示例（用 `\dontrun{}` 包裹） |
| `@author` | 作者信息（可选，推荐） |

### 4.2 示例规范

- 使用包内自带数据集：`data(df.xxx)`
- 用 `\dontrun{}` 包裹（避免 CRAN 检查问题）
- 展示最常用参数组合
- 内部函数示例如需引用自身，用 `bioRtools:::fn()`（`::` 对未导出函数无效）

---

## 五、图形函数规范

- 返回 ggplot 对象，**不直接打印**
- 图形参数使用描述性名称：`chr_colors` 而非 `colors`
- 使用 `theme_bio()` 或 `theme_bw()` 作为默认主题

```r
# Good：返回对象
p <- ggplot2::ggplot(data, ggplot2::aes(x, y)) +
  ggplot2::geom_point()
p
```

- 需要打开/关闭图形设备的函数：只对自己打开的设备执行
  `on.exit(grDevices::dev.off(), add = TRUE)`（见 6.3）

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

### 6.1 错误消息质量

- 消息要**具体、可操作**：点明哪个参数、期望什么类型、实际收到什么。
- 用户面函数抛错时给出修正提示。
- **消息语言统一为英文**（与现有代码一致；文档注释可用中文）。

```r
if (!is.numeric(data[[value_col]])) {
  stop("'value_col' must be numeric, got: ", class(data[[value_col]])[1])
}
```

bioRtools 全库统一用 `stop()`（不引入 cli），保持一致即可。

### 6.2 输入校验分层

- **用户面函数（导出）**：全面校验类型、长度、取值范围、列存在性。
- **内部函数**：假设输入已被上游校验，只断言关键不变量，避免重复开销。

### 6.3 避免全局状态副作用

不要修改全局状态。确需临时修改时，用 `on.exit()` 退出时还原。清单：

| 全局状态 | 规则 |
|----------|------|
| `options()` / `par()` | 保存旧值，`on.exit()` 还原 |
| `set.seed()` | 先捕获 `.Random.seed`，`on.exit()` 还原（见下例） |
| 图形设备 `pdf()/png()/...` | 仅对自己打开的设备 `on.exit(grDevices::dev.off(), add = TRUE)`；错误分支不要重复 `dev.off()` |
| 并行集群 `makeCluster()` | 紧跟 `on.exit(parallel::stopCluster(cl), add = TRUE)` |
| 用户 `.GlobalEnv` | 不写入、不删除用户的变量 |

```r
# par 还原
plot_helper <- function(...) {
  old_par <- graphics::par(mar = c(5, 4, 2, 2))
  on.exit(graphics::par(old_par), add = TRUE)
  # ... 绘图逻辑
}

# RNG 还原（capture_rng_state() / restore_rng_state() 为包内辅助函数）
old_seed <- capture_rng_state()
set.seed(seed)
on.exit(restore_rng_state(old_seed), add = TRUE)

# 集群清理
cl <- parallel::makeCluster(n_cores)
on.exit(parallel::stopCluster(cl), add = TRUE)
```

### 6.4 消息机制

- 进度、摘要、诊断信息用 `message()`（可被 `suppressMessages()` 捕获、便于测试）
- 可恢复问题用 `warning()`
- `cat()` 仅用于 `print.*()` / `summary.*()` 方法的输出
- 导出函数不得通过打印（`print()`/图形）产生唯一副作用，必须返回对象

---

## 七、测试规范

- 统一使用 **testthat 3e**：测试文件位于 `tests/testthat/test-*.R`，
  由 `tests/testthat.R` 入口运行（`usethis::use_testthat(3)` 生成）
- 测试不依赖工作目录：**禁止 `source("R/xxx.R")` 相对路径**（`R CMD check`
  的检查副本中会失败）；测试内调用导出函数，或 `pkgload::load_all()`
- 核心统计/转换函数（`calc_*`、`normalize_*`、`find_*`）必须有测试，
  关键输出可用 `tests/testthat/_snaps/` 或基线 rds 做回归
- 随机算法测试用固定 `set.seed()`，比较输出时注意浮点容差
- `devtools::test()` 必须能发现并执行全部测试；零测试通过视为未达标
- 现有 `tests/*.R` 裸脚本（cat + stopifnot 风格）逐步迁移为 testthat 用例

---

## 八、版本控制

- 遵循 [语义化版本](https://semver.org/)：`MAJOR.MINOR.PATCH`
- 每次修改更新 `DESCRIPTION` 版本号和 `CHANGELOG.md`
- Commit message 格式：`version X.Y.Z: brief description`
- 提交前运行 `devtools::document()` 确保文档同步
- **`Depends: R (>= x.y.z)` 必须真实**：下限不低于代码实际使用的语法/API
  （如使用 `\(x)` 匿名函数则下限至少 R 4.1.0）
- README 及文档中的版本徽章与 `DESCRIPTION` 保持一致
- 开发中功能（2.4）合并前由作者确认，未完成功能不进入发布

---

## 九、开发检查清单

### 新函数上线前

- [ ] 函数名 snake_case，动词开头
- [ ] 所有参数名 snake_case
- [ ] 输出列名 snake_case
- [ ] `@param` 每个参数有文档
- [ ] `@return` 说明返回结构
- [ ] `@examples` 可运行
- [ ] 参数验证完整
- [ ] 使用 `pkg::fun()` 显式命名空间（内部函数直呼，见 2.4）
- [ ] 新增导出前确认无功能近重复的已有函数
- [ ] 内部函数已加 `@keywords internal` + `@noRd`
- [ ] 进度/摘要用 `message()`，未用 `cat()`
- [ ] 已添加 testthat 测试
- [ ] `devtools::document()` 通过
- [ ] 包可正常 `pkgload::load_all()`

### 修改现有函数时

- [ ] 参数名改动需同步更新 roxygen `@param` 和 `@examples`
- [ ] 输出列名改动需同步更新下游函数引用
- [ ] 向后兼容函数保持旧接口不变
- [ ] 新增的 `pkg::` 引用已同步到 DESCRIPTION（见 §10）

### 发布版本前

- [ ] `devtools::document()` 文档已重新生成
- [ ] `devtools::check()` 0 error / 0 warning / 0 note
- [ ] `devtools::test()` 全部测试通过（且测试被真正发现）
- [ ] `urlchecker::url_check()` README 与文档链接有效
- [ ] DESCRIPTION 依赖声明与代码一致（无未声明的 `pkg::` 引用）
- [ ] `DESCRIPTION` 版本号按语义化版本升级（patch / minor / major）
- [ ] README 徽章版本号与 DESCRIPTION 同步
- [ ] `CHANGELOG.md` 已更新
- [ ] commit message 符合 `version X.Y.Z: ...` 格式

---

## 十、依赖管理

bioRtools 已有较多依赖（Imports 近 50 个），新增依赖需谨慎权衡。

### 10.1 依赖声明硬规则

**代码中任何包引用都必须出现在 DESCRIPTION 中**，按使用方式划分：

| 使用方式 | 声明位置 |
|----------|----------|
| `importFrom()`（NAMESPACE） | **Imports**（不可 Suggests，否则加载即失败） |
| `pkg::fun()` 无 `requireNamespace()` 守卫 | **Imports** |
| `pkg::fun()` 有 `requireNamespace("pkg")` 守卫 | **Suggests** |
| 仅文档 `\code{pkg::fun()}` 文字提及 | 无需声明 |
| 仅内部/WIP 函数使用 | 仍须声明，可放 Suggests |

> 历史教训：仓库曾存在约 35 个未声明的 `pkg::` 引用（Rtsne、readr、igraph、
> XML 等），`importFrom` 的 XML/minpack.lm/stats4/tidyselect/treeio 也漏写
> Imports，导致 `R CMD check` 大量 NOTE/WARNING。

### 10.2 何时引入新依赖

引入依赖应满足至少一条：

- 功能收益显著（自己实现成本高、易出错）
- 明显降低维护负担
- 涉及复杂实现（正则、日期解析、web 请求、生信专用算法）

简单工具函数（取长度、判断类型、基础算术）优先用 base R，不为小功能引入新包。

### 10.3 依赖类型划分（DESCRIPTION）

```
# 强依赖：包运行必需，写进 Imports，建议带版本下限
Imports:
    dplyr (>= 1.0.0),
    ggplot2 (>= 3.4.0)

# 可选依赖：仅测试 / vignette / 个别功能用到，写进 Suggests
Suggests:
    testthat (>= 3.0.0),
    knitr,
    rmarkdown
```

- 用 `requireNamespace("pkg", quietly = TRUE)` 检查 Suggests 依赖是否可用，再调用。
- 重型元包（`tidyverse` 整包、`shiny`）只在确有交互需求时引入，避免拖累安装。

---

## 参考资源

- [tidyverse 风格指南](https://style.tidyverse.org/) — 本项目主要参考
- [R Packages (Hadley Wickham)](https://r-pkgs.org/) — R 包开发权威指南
- [roxygen2 文档](https://roxygen2.r-lib.org/) — 文档生成工具
