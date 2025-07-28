test_that("pca_analysis works correctly", {
  skip_on_cran()
  
  # 准备测试数据
  data <- iris[,1:4]
  sample <- data.frame(
    sample = paste0("sample", 1:150),
    species = iris$Species
  )
  
  # 运行PCA分析
  result <- pca_analysis(data, sample)
  
  # 检查结果结构
  expect_true(is.list(result))
  expect_true(any(grepl("pca", names(result), ignore.case = TRUE)))
})
