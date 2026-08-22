#' Plot BSA (Bulked Segregant Analysis) QTL Results
#'
#' Reads BSA QTL results and returns a faceted plot of raw values,
#' smoothed fit, and threshold lines, plus a filtered table of
#' significant QTL regions.
#'
#' @param data A file path (character) or data frame from BSA analysis.
#'   Expected columns: Method, Chromosome, Position, Position_Mb,
#'   Raw_Value, Smooth_Value, Threshold, Above_Threshold.
#'
#' @return A named list containing:
#'   \describe{
#'     \item{\code{plot}}{ggplot2 object. Faceted scatter plot of Raw_Value
#'       with Smooth_Value line and Threshold reference line, faceted by
#'       Method and Chromosome.}
#'     \item{\code{table}}{Data frame of significant QTL regions
#'       (Above_Threshold == TRUE) with Chinese column names.}
#'   }
#'
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' \dontrun{
#' res <- bsa_plot_table("path/to/bsa_result.csv")
#' res$plot
#' res$table
#' }
bsa_plot_table <- function(data) {

  if (is.character(data)) {
    df <- read_data(data)
  } else if (is.data.frame(data)) {
    df <- data
  } else {
    stop("'data' must be a file path (character) or data frame")
  }

  required_cols <- c("Method", "Chromosome", "Position", "Raw_Value",
                     "Smooth_Value", "Threshold", "Above_Threshold")
  missing <- setdiff(required_cols, names(df))
  if (length(missing) > 0) {
    stop("Required columns missing: ", paste(missing, collapse = ", "))
  }

  # Natural sort chromosome order
  chr_levels <- stringr::str_sort(unique(df$Chromosome), numeric = TRUE)
  df <- df %>%
    dplyr::mutate(
      Position = .data$Position / 1000000,
      Chromosome = factor(.data$Chromosome, levels = chr_levels)
    )

  p <- df %>%
    ggplot2::ggplot() +
    ggplot2::geom_point(
      ggplot2::aes(.data$Position, .data$Raw_Value),
      color = "grey80"
    ) +
    # Method 已由 facet_grid 分面体现, hline 用固定色即可, 避免 color legend 与 strip 重复
    ggplot2::geom_hline(
      ggplot2::aes(yintercept = .data$Threshold),
      color = "red",
      linetype = "dashed"
    ) +
    ggplot2::geom_line(
      ggplot2::aes(.data$Position, .data$Smooth_Value),
      color = "blue"
    ) +
    ggplot2::facet_grid(Method ~ Chromosome, scales = "free") +
    ggplot2::labs(x = "Position (Mb)", y = "Value") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      strip.text.y = ggplot2::element_text(angle = 0),
      legend.position = "none"
    )

  table <- df %>%
    dplyr::filter(.data$Above_Threshold == TRUE) %>%
    magrittr::set_names(
      c("方法", "染色体", "位置_bp", "位置_Mb",
        "原始值", "拟合值", "方法阈值", "是否超过阈值")
    )

  list(plot = p, table = table)
}
