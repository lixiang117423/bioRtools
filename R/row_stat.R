#' Scale values to 0-1 range
#'
#' @description
#' This function scales a numeric vector to the range [0, 1] using min-max normalization.
#' The transformation is: (x - min(x)) / (max(x) - min(x))
#'
#' @param x A numeric vector to be scaled. Missing values (NA) are allowed.
#'
#' @return A numeric vector of the same length as input, scaled to [0, 1] range.
#'         If all values in x are identical, returns a vector of zeros.
#'         Missing values in input will result in missing values in output.
#'
#' @details
#' This function performs min-max normalization, which linearly transforms
#' the data to fit within the range [0, 1]. The minimum value becomes 0
#' and the maximum value becomes 1.
#'
#' @examples
#' # Basic usage
#' scale01(c(1, 2, 3, 4, 5))
#'
#' # With negative numbers
#' scale01(c(-2, 0, 2, 4))
#'
#' # With identical values
#' scale01(c(5, 5, 5))
#'
#' # With missing values
#' scale01(c(1, 2, NA, 4, 5))
#'
#' @author Xiang LI \email{lixiang117423@@foxmail.com}
#'
#' @seealso \code{\link{scale}} for standardization (z-score normalization)
#'
#' @export
scale01 <- function(x) {
  if (!is.numeric(x)) {
    stop("Input 'x' must be numeric")
  }

  if (length(x) == 0) {
    return(numeric(0))
  }

  if (all(is.na(x))) {
    return(x)
  }

  min_val <- min(x, na.rm = TRUE)
  max_val <- max(x, na.rm = TRUE)

  if (min_val == max_val) {
    result <- rep(0, length(x))
    result[is.na(x)] <- NA
    return(result)
  }

  (x - min_val) / (max_val - min_val)
}

#' Scale values across rows (rowwise scaling)
#'
#' @description
#' Apply min-max scaling (0-1) to each row independently across specified columns.
#' This is useful when you want to normalize values within each observation/row.
#' If no columns are specified, all numeric columns will be scaled.
#' Row names and column names are preserved.
#'
#' @param data A data frame (can have row names)
#' @param cols Columns to scale. Can be column names, indices, or tidyselect helpers.
#'        If NULL (default), all numeric columns will be scaled.
#' @param suffix Suffix to add to scaled column names (default: "_scaled")
#' @param replace Logical. If TRUE, replace original columns. If FALSE, create new columns (default: TRUE)
#'
#' @return A data frame with scaled values, preserving original row names and column names
#'
#' @examples
#' # Scale all numeric columns for each row (default behavior)
#' df <- data.frame(
#'   var1 = c(1, 5, 2),
#'   var2 = c(3, 2, 8),
#'   var3 = c(2, 7, 4),
#'   row.names = c("gene1", "gene2", "gene3")
#' )
#' scale01_rows(df)
#'
#' # Scale specific columns
#' scale01_rows(df, cols = c("var1", "var2"))
#'
#' @author Xiang LI \email{lixiang117423@@foxmail.com}
#'
#' @export
scale01_rows <- function(data, cols = NULL, suffix = "_scaled", replace = TRUE) {
  # 保存原始的行名和列名
  original_rownames <- rownames(data)
  original_colnames <- colnames(data)

  # 如果没有指定列，自动选择所有数值列
  if (is.null(cols)) {
    # 直接检查哪些列是数值型
    numeric_cols <- sapply(data, is.numeric)
    cols_to_scale <- names(data)[numeric_cols]
  } else {
    # 处理用户指定的列
    if (is.character(cols)) {
      cols_to_scale <- cols[cols %in% colnames(data)]
    } else if (is.numeric(cols)) {
      cols_to_scale <- colnames(data)[cols]
    } else {
      # 对于其他tidyselect语法，使用dplyr
      library(dplyr)
      cols_to_scale <- data %>%
        select({{ cols }}) %>%
        names()
    }
  }

  if (length(cols_to_scale) == 0) {
    warning("No numeric columns found to scale")
    return(data)
  }

  # 提取要标准化的数值数据
  numeric_data <- data[, cols_to_scale, drop = FALSE]

  # 对每行应用scale01
  scaled_matrix <- t(apply(numeric_data, 1, function(row) {
    if (all(is.na(row))) {
      return(row)  # 如果整行都是NA，返回原值
    }
    scale01(row)
  }))

  # 确保scaled_matrix保持矩阵格式并设置正确的维度名
  if (is.vector(scaled_matrix)) {
    # 当只有一行时，t()会返回向量，需要转换回矩阵
    scaled_matrix <- matrix(scaled_matrix, nrow = 1)
  }

  # 设置行名和列名
  rownames(scaled_matrix) <- original_rownames
  colnames(scaled_matrix) <- cols_to_scale

  # 转换为数据框
  scaled_df <- as.data.frame(scaled_matrix)

  # 创建结果数据框
  if (replace) {
    # 替换原列
    result <- data
    result[, cols_to_scale] <- scaled_df
  } else {
    # 添加新列
    result <- data
    new_names <- paste0(cols_to_scale, suffix)
    result[, new_names] <- scaled_df
  }

  # 确保保留原始行名
  rownames(result) <- original_rownames

  return(result)
}

#' Scale values within groups using dplyr
#'
#' @description
#' Apply min-max scaling (0-1) to specified columns within each group.
#' This function is designed to work seamlessly with dplyr's group_by().
#'
#' @param data A grouped or ungrouped data frame
#' @param cols Columns to scale within each group
#' @param suffix Suffix for new scaled columns (default: "_scaled")
#' @param replace Logical. Replace original columns or create new ones (default: FALSE)
#'
#' @return A data frame with within-group scaled values
#'
#' @examples
#' library(dplyr)
#'
#' # Group by category and scale within each group
#' df <- data.frame(
#'   category = rep(c("A", "B"), each = 3),
#'   value1 = c(1, 3, 5, 2, 4, 6),
#'   value2 = c(10, 30, 50, 20, 40, 60)
#' )
#'
#' df %>%
#'   group_by(category) %>%
#'   scale01_groups(cols = c("value1", "value2"))
#'
#' @author Xiang LI \email{lixiang117423@@foxmail.com}
#'
#' @export
scale01_groups <- function(data, cols, suffix = "_scaled", replace = FALSE) {
  library(dplyr)

  result <- data %>%
    mutate(
      across({{ cols }},
        scale01,
        .names = if (replace) "{.col}" else paste0("{.col}", suffix))
    )

  return(result)
}

#' Convenient function for scaling in dplyr pipelines
#'
#' @description
#' A wrapper function that makes it easy to scale multiple columns
#' in dplyr pipelines with automatic "_scaled" suffix.
#'
#' @param .data A data frame (can be grouped)
#' @param ... Column names to scale (unquoted)
#'
#' @return A data frame with additional scaled columns
#'
#' @examples
#' library(dplyr)
#'
#' df <- data.frame(
#'   group = rep(c("A", "B"), each = 3),
#'   var1 = c(1, 3, 5, 2, 4, 6),
#'   var2 = c(10, 30, 50, 20, 40, 60)
#' )
#'
#' # Scale within groups
#' df %>%
#'   group_by(group) %>%
#'   mutate_scale01(var1, var2)
#'
#' @author Xiang LI \email{lixiang117423@@foxmail.com}
#'
#' @export
mutate_scale01 <- function(.data, ...) {
  .data %>% mutate(across(c(...), scale01, .names = "{.col}_scaled"))
}

#' Convenient function for scaling with custom column names
#'
#' @description
#' A wrapper function that allows custom naming when scaling columns
#' in dplyr pipelines.
#'
#' @param .data A data frame (can be grouped)
#' @param ... Named expressions for scaling (e.g., new_name = scale01(old_col))
#'
#' @return A data frame with custom named scaled columns
#'
#' @examples
#' library(dplyr)
#'
#' df <- data.frame(
#'   category = rep(c("A", "B"), each = 3),
#'   revenue = c(100, 300, 500, 200, 400, 600),
#'   cost = c(50, 150, 250, 100, 200, 300)
#' )
#'
#' # Custom names for scaled columns
#' df %>%
#'   group_by(category) %>%
#'   mutate_scale01_named(
#'     revenue_norm = scale01(revenue),
#'     cost_norm = scale01(cost)
#'   )
#'
#' @author Xiang LI \email{lixiang117423@@foxmail.com}
#'
#' @export
mutate_scale01_named <- function(.data, ...) {
  .data %>% mutate(...)
}

#' Calculate row-wise standard deviation
#'
#' @description
#' Calculate the standard deviation for each row across specified columns.
#' This function is designed to work seamlessly with dplyr::mutate().
#'
#' @param data A data frame (can have row names)
#' @param cols Columns to calculate standard deviation across. Can be column names,
#'        indices, or tidyselect helpers. If NULL (default), all numeric columns will be used.
#' @param na.rm Logical. Should missing values be removed? (default: TRUE)
#'
#' @return A numeric vector of standard deviations with names matching the row names.
#'
#' @examples
#' # Use in dplyr::mutate()
#' library(dplyr)
#' df <- data.frame(
#'   sample1 = c(100, 23, 95),
#'   sample2 = c(55, 24, 6),
#'   sample3 = c(93, 20, 59),
#'   row.names = c("Gene1", "Gene2", "Gene3")
#' )
#'
#' # Add standard deviation column using mutate
#' df %>% mutate(sd = row_sd(.))
#'
#' # Calculate SD for specific columns
#' df %>% mutate(sd = row_sd(., cols = c("sample1", "sample2")))
#'
#' @author Xiang LI \email{lixiang117423@@foxmail.com}
#'
#' @export
row_sd <- function(data, cols = NULL, na.rm = TRUE) {
  # 保存原始的行名
  original_rownames <- rownames(data)

  # 如果没有指定列，自动选择所有数值列
  if (is.null(cols)) {
    numeric_cols <- sapply(data, is.numeric)
    cols_to_use <- names(data)[numeric_cols]
  } else {
    # 处理用户指定的列
    if (is.character(cols)) {
      cols_to_use <- cols[cols %in% colnames(data)]
    } else if (is.numeric(cols)) {
      cols_to_use <- colnames(data)[cols]
    } else {
      # 对于其他tidyselect语法，使用dplyr
      library(dplyr)
      cols_to_use <- data %>%
        select({{ cols }}) %>%
        names()
    }
  }

  if (length(cols_to_use) == 0) {
    stop("No numeric columns found to calculate standard deviation")
  }

  if (length(cols_to_use) == 1) {
    warning("Only one column selected. Standard deviation will be 0 for all rows.")
  }

  # 提取数值数据
  numeric_data <- data[, cols_to_use, drop = FALSE]

  # 计算每行的标准差
  row_values <- apply(numeric_data, 1, function(row) {
    if (na.rm) {
      valid_values <- row[!is.na(row)]
      if (length(valid_values) <= 1) {
        return(NA)
      }
      return(sd(valid_values))
    } else {
      if (any(is.na(row))) {
        return(NA)
      }
      return(sd(row))
    }
  })

  # 设置行名
  names(row_values) <- original_rownames

  return(row_values)
}

#' Calculate row-wise coefficient of variation (CV)
#'
#' @description
#' Calculate the coefficient of variation (CV = standard deviation / mean) for each row.
#' This function is designed to work seamlessly with dplyr::mutate().
#'
#' @param data A data frame (can have row names)
#' @param cols Columns to calculate CV across. If NULL (default), all numeric columns will be used.
#' @param na.rm Logical. Should missing values be removed? (default: TRUE)
#'
#' @return A numeric vector of CVs with names matching the row names.
#'
#' @examples
#' library(dplyr)
#' df <- data.frame(
#'   sample1 = c(100, 23, 95),
#'   sample2 = c(55, 24, 6),
#'   sample3 = c(93, 20, 59),
#'   row.names = c("Gene1", "Gene2", "Gene3")
#' )
#'
#' # Add CV column using mutate
#' df %>% mutate(cv = row_cv(.))
#'
#' @author Xiang LI \email{lixiang117423@@foxmail.com}
#'
#' @export
row_cv <- function(data, cols = NULL, na.rm = TRUE) {
  # 保存原始的行名
  original_rownames <- rownames(data)

  # 如果没有指定列，自动选择所有数值列
  if (is.null(cols)) {
    numeric_cols <- sapply(data, is.numeric)
    cols_to_use <- names(data)[numeric_cols]
  } else {
    # 处理用户指定的列
    if (is.character(cols)) {
      cols_to_use <- cols[cols %in% colnames(data)]
    } else if (is.numeric(cols)) {
      cols_to_use <- colnames(data)[cols]
    } else {
      # 对于其他tidyselect语法，使用dplyr
      library(dplyr)
      cols_to_use <- data %>%
        select({{ cols }}) %>%
        names()
    }
  }

  if (length(cols_to_use) == 0) {
    stop("No numeric columns found to calculate coefficient of variation")
  }

  # 提取数值数据
  numeric_data <- data[, cols_to_use, drop = FALSE]

  # 计算每行的CV
  row_values <- apply(numeric_data, 1, function(row) {
    if (na.rm) {
      valid_values <- row[!is.na(row)]
      if (length(valid_values) <= 1) {
        return(NA)
      }
      mean_val <- mean(valid_values)
      if (mean_val == 0) {
        return(Inf)
      }
      return(sd(valid_values) / mean_val)
    } else {
      if (any(is.na(row))) {
        return(NA)
      }
      mean_val <- mean(row)
      if (mean_val == 0) {
        return(Inf)
      }
      return(sd(row) / mean_val)
    }
  })

  # 设置行名
  names(row_values) <- original_rownames

  return(row_values)
}

#' Calculate row-wise mean
#'
#' @description
#' Calculate the mean for each row across specified columns.
#' This function is designed to work seamlessly with dplyr::mutate().
#'
#' @param data A data frame (can have row names)
#' @param cols Columns to calculate mean across. If NULL (default), all numeric columns will be used.
#' @param na.rm Logical. Should missing values be removed? (default: TRUE)
#'
#' @return A numeric vector of means with names matching the row names.
#'
#' @examples
#' library(dplyr)
#' df <- data.frame(
#'   sample1 = c(100, 23, 95),
#'   sample2 = c(55, 24, 6),
#'   sample3 = c(93, 20, 59),
#'   row.names = c("Gene1", "Gene2", "Gene3")
#' )
#'
#' # Add mean column using mutate
#' df %>% mutate(mean_expr = row_mean(.))
#'
#' @author Xiang LI \email{lixiang117423@@foxmail.com}
#'
#' @export
row_mean <- function(data, cols = NULL, na.rm = TRUE) {
  # 保存原始的行名
  original_rownames <- rownames(data)

  # 如果没有指定列，自动选择所有数值列
  if (is.null(cols)) {
    numeric_cols <- sapply(data, is.numeric)
    cols_to_use <- names(data)[numeric_cols]
  } else {
    # 处理用户指定的列
    if (is.character(cols)) {
      cols_to_use <- cols[cols %in% colnames(data)]
    } else if (is.numeric(cols)) {
      cols_to_use <- colnames(data)[cols]
    } else {
      # 对于其他tidyselect语法，使用dplyr
      library(dplyr)
      cols_to_use <- data %>%
        select({{ cols }}) %>%
        names()
    }
  }

  if (length(cols_to_use) == 0) {
    stop("No numeric columns found to calculate mean")
  }

  # 提取数值数据
  numeric_data <- data[, cols_to_use, drop = FALSE]

  # 计算每行的均值
  row_values <- apply(numeric_data, 1, function(row) {
    mean(row, na.rm = na.rm)
  })

  # 设置行名
  names(row_values) <- original_rownames

  return(row_values)
}

#' Calculate row-wise maximum
#'
#' @description
#' Calculate the maximum value for each row across specified columns.
#' This function is designed to work seamlessly with dplyr::mutate().
#'
#' @param data A data frame (can have row names)
#' @param cols Columns to calculate maximum across. If NULL (default), all numeric columns will be used.
#' @param na.rm Logical. Should missing values be removed? (default: TRUE)
#'
#' @return A numeric vector of maximums with names matching the row names.
#'
#' @examples
#' library(dplyr)
#' df <- data.frame(
#'   sample1 = c(100, 23, 95),
#'   sample2 = c(55, 24, 6),
#'   sample3 = c(93, 20, 59),
#'   row.names = c("Gene1", "Gene2", "Gene3")
#' )
#'
#' # Add maximum column using mutate
#' df %>% mutate(max_expr = row_max(.))
#'
#' @author Xiang LI \email{lixiang117423@@foxmail.com}
#'
#' @export
row_max <- function(data, cols = NULL, na.rm = TRUE) {
  # 保存原始的行名
  original_rownames <- rownames(data)

  # 如果没有指定列，自动选择所有数值列
  if (is.null(cols)) {
    numeric_cols <- sapply(data, is.numeric)
    cols_to_use <- names(data)[numeric_cols]
  } else {
    # 处理用户指定的列
    if (is.character(cols)) {
      cols_to_use <- cols[cols %in% colnames(data)]
    } else if (is.numeric(cols)) {
      cols_to_use <- colnames(data)[cols]
    } else {
      # 对于其他tidyselect语法，使用dplyr
      library(dplyr)
      cols_to_use <- data %>%
        select({{ cols }}) %>%
        names()
    }
  }

  if (length(cols_to_use) == 0) {
    stop("No numeric columns found to calculate maximum")
  }

  # 提取数值数据
  numeric_data <- data[, cols_to_use, drop = FALSE]

  # 计算每行的最大值
  row_values <- apply(numeric_data, 1, function(row) {
    max(row, na.rm = na.rm)
  })

  # 设置行名
  names(row_values) <- original_rownames

  return(row_values)
}

#' Calculate row-wise minimum
#'
#' @description
#' Calculate the minimum value for each row across specified columns.
#' This function is designed to work seamlessly with dplyr::mutate().
#'
#' @param data A data frame (can have row names)
#' @param cols Columns to calculate minimum across. If NULL (default), all numeric columns will be used.
#' @param na.rm Logical. Should missing values be removed? (default: TRUE)
#'
#' @return A numeric vector of minimums with names matching the row names.
#'
#' @examples
#' library(dplyr)
#' df <- data.frame(
#'   sample1 = c(100, 23, 95),
#'   sample2 = c(55, 24, 6),
#'   sample3 = c(93, 20, 59),
#'   row.names = c("Gene1", "Gene2", "Gene3")
#' )
#'
#' # Add minimum column using mutate
#' df %>% mutate(min_expr = row_min(.))
#'
#' @author Xiang LI \email{lixiang117423@@foxmail.com}
#'
#' @export
row_min <- function(data, cols = NULL, na.rm = TRUE) {
  # 保存原始的行名
  original_rownames <- rownames(data)

  # 如果没有指定列，自动选择所有数值列
  if (is.null(cols)) {
    numeric_cols <- sapply(data, is.numeric)
    cols_to_use <- names(data)[numeric_cols]
  } else {
    # 处理用户指定的列
    if (is.character(cols)) {
      cols_to_use <- cols[cols %in% colnames(data)]
    } else if (is.numeric(cols)) {
      cols_to_use <- colnames(data)[cols]
    } else {
      # 对于其他tidyselect语法，使用dplyr
      library(dplyr)
      cols_to_use <- data %>%
        select({{ cols }}) %>%
        names()
    }
  }

  if (length(cols_to_use) == 0) {
    stop("No numeric columns found to calculate minimum")
  }

  # 提取数值数据
  numeric_data <- data[, cols_to_use, drop = FALSE]

  # 计算每行的最小值
  row_values <- apply(numeric_data, 1, function(row) {
    min(row, na.rm = na.rm)
  })

  # 设置行名
  names(row_values) <- original_rownames

  return(row_values)
}
