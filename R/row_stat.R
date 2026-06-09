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
      cols_to_use <- data %>%
        dplyr::select({{ cols }}) %>%
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

  row_values
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
      cols_to_use <- data %>%
        dplyr::select({{ cols }}) %>%
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

  row_values
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
      cols_to_use <- data %>%
        dplyr::select({{ cols }}) %>%
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

  row_values
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
      cols_to_use <- data %>%
        dplyr::select({{ cols }}) %>%
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

  row_values
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
      cols_to_use <- data %>%
        dplyr::select({{ cols }}) %>%
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

  row_values
}
