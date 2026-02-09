# BioRtools 开发脚本使用指南

## 版本: 1.0
## 日期: 2026-02-09

---

## 概述

bioRtools 包包含一套完整的开发工具脚本，用于自动化常见的R包开发任务。

---

## 文件结构

```
scripts/
├── setup.R                    # 主脚本（包含所有功能）
├── setup/
│   ├── 01_setup_basic.R      # 基础配置模块
│   └── 02_document.R         # 文档和构建模块
└── README.md                  # 本文档
```

---

## 快速开始

### 第一次设置

```r
# 在R中加载主脚本
source("scripts/setup.R")

# 运行完整设置
setup_package(
  github_user = "your_github_username",
  github_repo = "bioRtools",
  license = "mit"
)
```

这将自动完成：
1. ✅ 配置许可证（MIT）
2. ✅ 设置GitHub链接
3. ✅ 配置 .Rbuildignore 和 .gitignore
4. ✅ 添加常用包的导入声明
5. ✅ 格式化所有代码（styler）
6. ✅ 生成文档（roxygen2）
7. ✅ 构建并安装包

---

## 日常开发工作流

### 场景1: 添加新函数后

```r
source("scripts/setup.R")

# 格式化代码并更新文档
format_and_document()

# 安装更新的包
quick_install()
```

### 场景2: 更新包版本

```r
source("scripts/setup.R")

# 更新为新的次版本（添加新功能）
update_version("minor")

# 或更新为补丁版本（bug修复）
update_version("patch")

# 或更新为主版本（不兼容的API更改）
update_version("major")
```

### 场景3: 准备发布

```r
source("scripts/setup.R")

# 完整设置（包括构建）
setup_package(
  github_user = "yourname",
  github_repo = "bioRtools",
  build = TRUE,
  install = TRUE
)

# 运行测试
run_package_tests()

# 检查包（可选）
devtools::check()
```

---

## 单独功能说明

### 1. setup_basic_config()

配置包的基础设置。

**参数：**
- `github_user` - GitHub用户名
- `github_repo` - GitHub仓库名
- `license` - 许可证类型（"mit", "apache", "gpl"）

**示例：**
```r
setup_basic_config(
  github_user = "lixiang117423",
  github_repo = "bioRtools",
  license = "mit"
)
```

---

### 2. add_package_imports()

添加 `@importFrom` 声明到NAMESPACE。

**参数：**
- `imports` - 导入声明向量（格式：c("dplyr::mutate", "ggplot2::aes")）
           如果为NULL，使用默认集

**示例：**
```r
# 使用默认导入
add_package_imports()

# 添加自定义导入
add_package_imports(c(
  "dplyr::mutate",
  "ggplot2::aes",
  "tidyr::pivot_longer"
))
```

---

### 3. format_package_code()

使用 styler 格式化所有R代码。

**参数：**
- `strict` - 是否使用严格模式（默认：FALSE）

**示例：**
```r
# 格式化所有代码
format_package_code()

# 使用严格模式格式化
format_package_code(strict = TRUE)
```

---

### 4. generate_documentation()

生成 roxygen2 文档。

**参数：**
- `check` - 是否运行 `devtools::check()`（默认：FALSE）

**示例：**
```r
# 仅生成文档
generate_documentation()

# 生成文档并检查
generate_documentation(check = TRUE)
```

---

### 5. build_install_package()

构建和安装包。

**参数：**
- `build` - 是否构建tarball（默认：TRUE）
- `install` - 是否安装包（默认：TRUE）
- `force` - 是否强制重装（默认：TRUE）
- `check` - 是否检查构建的tarball（默认：FALSE）

**示例：**
```r
# 构建并安装
build_install_package()

# 仅安装
build_install_package(build = FALSE)

# 构建、安装并检查
build_install_package(check = TRUE)
```

---

### 6. update_version()

更新包版本号。

**参数：**
- `type` - 版本类型（"major", "minor", "patch", "dev"）

**版本号规则：**
- `major` - 主版本（不兼容的API更改）：1.0.0 → 2.0.0
- `minor` - 次版本（向后兼容的新功能）：1.0.0 → 1.1.0
- `patch` - 补丁版本（bug修复）：1.0.0 → 1.0.1
- `dev` - 开发版本：1.0.0.9000

**示例：**
```r
# 开发版本（默认）
update_version("dev")

# 补丁版本（bug修复）
update_version("patch")

# 次版本（新功能）
update_version("minor")

# 主版本（重大更改）
update_version("major")
```

---

### 7. run_package_tests()

运行包的测试套件。

**参数：**
- `filter` - 测试文件过滤模式
- `reporter` - 报告类型（"summary", "location", "failonly"）

**示例：**
```r
# 运行所有测试
run_package_tests()

# 运行特定测试
run_package_tests(filter = "test_plot_multi_volcano")
```

---

### 8. setup_package()

主设置函数，执行所有步骤。

**参数：**
- `github_user` - GitHub用户名
- `github_repo` - GitHub仓库名
- `license` - 许可证类型
- `build` - 是否构建包
- `install` - 是否安装包
- `format_code` - 是否格式化代码

**示例：**
```r
# 完整设置
setup_package(
  github_user = "yourname",
  github_repo = "bioRtools"
)

# 仅配置，不构建
setup_package(
  github_user = "yourname",
  github_repo = "bioRtools",
  build = FALSE,
  install = FALSE
)
```

---

## 快捷功能

### format_and_document()

快速格式化代码并更新文档。

```r
source("scripts/setup.R")
format_and_document()
```

等同于：
```r
format_package_code()
generate_documentation()
```

---

### quick_install()

快速重新生成文档并安装包。

```r
source("scripts/setup.R")
quick_install()
```

等同于：
```r
generate_documentation(check = FALSE)
devtools::install_local(force = TRUE)
```

---

## 开发工作流示例

### 典型的功能开发流程

```r
# 1. 加载开发工具
source("scripts/setup.R")

# 2. 创建/修改函数
# ... 编辑 R/new_function.R ...

# 3. 格式化和文档化
format_and_document()

# 4. 测试功能
quick_install()

# 5. 提交到git
# git add .
# git commit -m "feat: add new function"
```

### 发布前检查流程

```r
# 1. 加载开发工具
source("scripts/setup.R")

# 2. 更新版本号
update_version("minor")  # 或 "patch"

# 3. 完整设置
setup_package(
  github_user = "yourname",
  github_repo = "bioRtools",
  build = TRUE,
  install = TRUE
)

# 4. 运行所有测试
run_package_tests()

# 5. CRAN检查（如果发布到CRAN）
devtools::check()
```

---

## 故障排除

### 问题：导入声明冲突

**错误：** `namespace import conflict`

**解决：** 删除NAMESPACE文件，重新运行：
```r
file.remove("NAMESPACE")
generate_documentation()
```

---

### 问题：文档生成失败

**错误：** `@return has mismatched braces`

**解决：** 检查roxygen2注释中的括号匹配，然后重新运行：
```r
generate_documentation()
```

---

### 问题：styler格式化问题

**解决：** Styler会自动处理大多数情况。如果有问题，可以：
1. 手动修复特定文件
2. 使用 `strict = TRUE` 进行严格格式化
3. 跳过格式化：`setup_package(format_code = FALSE)`

---

## 最佳实践

1. **定期运行格式化**
   - 每次开发会话后运行 `format_and_document()`
   - 提交前运行一次

2. **版本管理**
   - 开发中使用 `update_version("dev")`
   - 新功能使用 `update_version("minor")`
   - Bug修复使用 `update_version("patch")`

3. **文档优先**
   - 修改函数后立即更新文档注释
   - 运行 `generate_documentation()` 生成 .Rd 文件

4. **测试驱动**
   - 添加功能后立即创建测试
   - 使用 `run_package_tests()` 验证

5. **Git工作流**
   - 在提交前运行 `format_and_document()`
   - 提交信息遵循约定提交格式

---

## 依赖要求

这些脚本需要以下R包：

- `usethis` - 包设置工具
- `styler` - 代码格式化
- `roxygen2` - 文档生成
- `devtools` - 开发工具

如果未安装，脚本会自动安装。

---

## 更新日志

| 版本 | 日期 | 更改内容 |
|------|------|----------|
| 1.0 | 2026-02-09 | 初始版本 - 完全重构开发脚本 |

---

**维护者:** Xiang LI \email{lixiang117423@@foxmail.com}
