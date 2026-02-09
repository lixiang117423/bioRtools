# bioRtools 开发环境配置指南

## 版本: 1.0
## 日期: 2026-02-09

---

## 一、系统信息

### 当前开发环境

- **OS**: macOS (Darwin 25.2.0)
- **R路径**: `/opt/homebrew/bin/R`
- **R版本**: 请运行 `R --version` 查看
- **包路径**: `/Users/lixiang/NutstoreFiles/03.编程相关/bioRtools`

---

## 二、R环境配置

### 2.1 使用Homebrew安装的R

如果你的R是通过Homebrew安装的，路径通常是：
```bash
/opt/homebrew/bin/R  # Apple Silicon Mac
/usr/local/bin/R     # Intel Mac
```

**验证R路径：**
```bash
which R
# 应该显示: /opt/homebrew/bin/R
```

**设置别名（可选）：**
如果需要，可以在 `~/.zshrc` 或 `~/.bash_profile` 中添加：
```bash
alias R='/opt/homebrew/bin/R'
```

### 2.2 在脚本中使用R路径

**方法1：直接指定完整路径**
```bash
/opt/homebrew/bin/Rscript script.R
```

**方法2：使用which R（推荐）**
```bash
#!/bin/bash
R_PATH=$(which R)
$R_PATH --no-save --no-restore -e "source('script.R')"
```

**方法3：假设R在PATH中**
```bash
#!/usr/bin/env Rscript
# script content
```

---

## 三、必需的R包

### 3.1 核心依赖包

```r
# 数据操作
install.packages(c("dplyr", "tidyr", "tibble", "magrittr"))

# 可视化
install.packages(c("ggplot2", "ggprism", "ggrepel", "ggthemes"))

# 开发工具
install.packages(c("devtools", "roxygen2", "testthat", "lintr", "styler"))

# 其他工具
install.packages(c("readr", "rlang", "purrr"))
```

### 3.2 Bioconductor包（如需要）

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("GenomicRanges", "BiocGenerics"))
```

---

## 四、开发工作流

### 4.1 生成文档

```bash
# 方法1：使用完整R路径
/opt/homebrew/bin/R -e "library(roxygen2); roxygenise()"

# 方法2：使用Rscript
/opt/homebrew/bin/Rscript -e "library(roxygen2); roxygenise()"
```

### 4.2 运行测试

```bash
# 运行特定测试脚本
/opt/homebrew/bin/Rscript tests/test_plot_multi_volcano.R

# 使用testthat运行所有测试
/opt/homebrew/bin/R -e "library(testthat); test()"
```

### 4.3 安装开发版本

```bash
# 从当前目录安装
/opt/homebrew/bin/R -e "devtools::install()"

# 或使用Rscript
/opt/homebrew/bin/Rscript -e "devtools::install()"
```

### 4.4 代码检查

```bash
# 使用lintr检查代码风格
/opt/homebrew/bin/Rscript -e "lintr::lint_dir('R/')"

# 使用styler自动格式化代码
/opt/homebrew/bin/Rscript -e "styler::style_dir('R/')"
```

---

## 五、IDE配置

### 5.1 RStudio配置

如果使用RStudio：

1. **设置R路径**：
   - RStudio → Preferences → R → R General
   - R Version：选择 "Custom"
   - 浏览到 `/opt/homebrew/bin/R`

2. **配置Build Tools**：
   - RStudio → Preferences → Tools → Build Tools
   - 确保R路径正确

### 5.2 VS Code配置

如果使用VS Code + R扩展：

1. **安装R扩展**：
   - 搜索 "Yuki Ueda" 的 R扩展

2. **配置R路径**：
   在 `settings.json` 中添加：
   ```json
   {
     "r.rpath.linux": "/opt/homebrew/bin/R",
     "r.rpath.mac": "/opt/homebrew/bin/R"
   }
   ```

---

## 六、Git配置

### 6.1 建议的.gitignore

```
# R specific
.Rproj.user
.Rhistory
.RData
.Ruserdata

# Package build
*.tar.gz
*.Rcheck

# Documentation
docs/

# Test outputs
test_*.png
test_*.pdf

# Temporary files
*~
*.swp
*.swo
```

### 6.2 提交前检查清单

- [ ] 运行 `roxygen2::roxygenise()` 更新文档
- [ ] 运行 `testthat::test()` 确保测试通过
- [ ] 运行 `lintr::lint_dir("R/")` 检查代码风格
- [ ] 更新NAMESPACE（roxygenise会自动完成）
- [ ] 检查DESCRIPTION文件
- [ ] Git commit所有更改

---

## 七、常见问题

### Q1: 如何在命令行中运行R脚本？

**A:** 使用Rscript而不是直接调用R：
```bash
# Good
/opt/homebrew/bin/Rscript script.R

# Bad（会启动交互式会话）
/opt/homebrew/bin/R < script.R
```

### Q2: 如何处理包的依赖？

**A:** 在DESCRIPTION文件中明确声明：
```
Imports:
    dplyr,
    ggplot2,
    rlang
```

### Q3: 如何测试单个函数？

**A:** 创建独立的测试脚本：
```bash
/opt/homebrew/bin/Rscript -e "
source('R/myfunction.R')
# 测试代码
result <- my_function(test_data)
print(result)
"
```

### Q4: 如何在不同环境之间切换R版本？

**A:** 使用符号链接或修改PATH：
```bash
# 临时切换
export PATH="/usr/local/bin/R:$PATH"

# 永久切换（在 ~/.zshrc 中）
export PATH="/opt/homebrew/bin:$PATH"
```

---

## 八、快速开始

### 第一次设置

```bash
# 1. 克隆或进入项目目录
cd /Users/lixiang/NutstoreFiles/03.编程相关/bioRtools

# 2. 安装开发工具
/opt/homebrew/bin/R -e "install.packages(c('devtools', 'roxygen2', 'testthat'))"

# 3. 安装包依赖
/opt/homebrew/bin/R -e "devtools::install_dev_deps()"

# 4. 生成文档
/opt/homebrew/bin/R -e "library(roxygen2); roxygenise()"

# 5. 安装包
/opt/homebrew/bin/R -e "devtools::install()"

# 6. 运行测试
/opt/homebrew/bin/Rscript tests/test_plot_multi_volcano.R
```

### 日常工作流

```bash
# 1. 编辑代码
vim R/new_function.R

# 2. 生成文档
/opt/homebrew/bin/R -e "library(roxygen2); roxygenise()"

# 3. 测试
/opt/homebrew/bin/Rscript tests/test_new_function.R

# 4. 提交
git add .
git commit -m "feat: add new function"
```

---

**最后更新**: 2026-02-09
**维护者**: Xiang LI
