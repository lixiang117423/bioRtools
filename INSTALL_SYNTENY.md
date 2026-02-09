# 安装和测试 plot_synteny 函数的步骤

## 1. 安装依赖包

```r
# 安装必需的包
install.packages("gggenes")
install.packages("ggforce")

# 或者从GitHub安装开发版本（如果有）
remotes::install_github("wilkox/gggenes")
remotes::install_github("thomasp85/ggforce")
```

## 2. 生成文档和NAMESPACE

```r
# 设置工作目录到包的根目录
setwd("/Users/lixiang/NutstoreFiles/03.编程相关/bioRtools")

# 加载devtools
library(devtools)

# 生成文档（这会创建/更新 NAMESPACE 文件）
document()

# 或者使用函数
devtools::document()
```

## 3. 重新加载包

```r
# 方法1: 使用 devtools
devtools::load_all()

# 方法2: 安装包
devtools::install()

# 然后正常加载
library(bioRtools)
```

## 4. 测试数据导出

```r
# 检查数据是否可以访问
data(df.synteny.gene)
data(df.synteny.link)

# 查看数据
head(df.synteny.gene)
head(df.synteny.link)

# 测试绘图
result <- plot_synteny(
  gene_data = df.synteny.gene,
  syntenic_data = df.synteny.link,
  species_labels = c("A", "B", "C", "D")
)

# 显示图形
result$plot.synteny

# 查看统计信息
result$data.summary
```

## 常见问题解决

### 问题1: 'df.synteny.gene' is not an exported object

**原因**: NAMESPACE文件没有正确导出数据对象

**解决方案**:
```r
# 1. 确保 data.R 文件中有正确的 roxygen2 注释
# 2. 运行 devtools::document() 生成文档
# 3. 检查 NAMESPACE 文件是否包含以下内容：
#    export(df.synteny.gene)
#    export(df.synteny.link)

# 如果没有，手动添加到 NAMESPACE 文件
# 或者重新运行 devtools::document()
```

### 问题2: 找不到 gggenes 包

**解决方案**:
```r
# gggenes 包可能不在 CRAN，需要从 GitHub 安装
remotes::install_github("wilkox/gggenes")

# 或者先安装 remotes
install.packages("remotes")
```

### 问题3: 数据文件加载失败

**解决方案**:
```r
# 确保数据文件在正确的位置
list.files("data", pattern = "synteny")

# 应该看到：
# df.synteny.gene.rda
# df.synteny.link.rda

# 如果文件存在但加载失败，尝试重新创建数据
source("tests/create_synteny_data.R")
```

## 完整测试流程

```r
# 清理环境
rm(list = ls())

# 1. 安装依赖
install.packages(c("gggenes", "ggforce", "devtools"))

# 2. 设置工作目录
setwd("/Users/lixiang/NutstoreFiles/03.编程相关/bioRtools")

# 3. 生成文档
devtools::document()

# 4. 加载包
devtools::load_all()

# 5. 测试数据加载
data(df.synteny.gene)
data(df.synteny.link)

# 6. 测试函数
result <- plot_synteny(
  gene_data = df.synteny.gene,
  syntenic_data = df.synteny.link,
  species_labels = c("A", "B", "C", "D")
)

# 7. 查看结果
print(result$plot.synteny)
print(result$data.summary)

# 8. 保存测试图
ggsave("synteny_test.png", plot = result$plot.synteny,
       width = 10, height = 8, dpi = 300)
```

## 验证安装成功

```r
# 检查函数是否可用
ls("package:bioRtools", pattern = "plot_synteny")

# 检查数据是否可用
ls("package:bioRtools", pattern = "df.synteny")

# 应该看到：
# [1] "plot_synteny"
# [2] "df.synteny.gene"
# [3] "df.synteny.link"
```

## 下一步

测试成功后，你可以：

1. 使用自己的数据
2. 自定义颜色和参数
3. 添加到你的分析流程中
4. 为函数编写更多的测试用例
