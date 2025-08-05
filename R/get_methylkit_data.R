# 确保已加载必要的包
# 如果尚未安装，请运行 install.packages("methylKit") 和 install.packages("magrittr")
library(methylKit)
library(magrittr) # 用于 %>% 管道操作符

#' @title 从methylKit对象提取全面的甲基化数据（包含上下文信息）
#' @description 该函数接收一个methylKit对象，进行过滤和标准化，然后生成一个包含
#'              完整信息的宽格式数据框。对于每个样本，它会保留原始的覆盖度(coverage)、
#'              甲基化C(numCs)、未甲基化C(numTs)计数，并计算甲基化百分比。
#'              所有与样本相关的列都将使用'样本名_数据类型'的格式进行命名。
#'
#' @param methyl_obj 一个 `methylRawList` 对象。
#' @param context 字符型。甲基化的序列上下文（如 "CpG", "CHG", 或 "CHH"）。
#'                此信息将被添加为输出数据框的一列。默认为 "CpG"。
#' @param lo.count 数值。最低覆盖度阈值，默认为 10。
#' @param lo.perc 数值。要过滤掉的覆盖度值的较低百分位数。默认为 NULL。
#' @param hi.count 数值。用于过滤的最高覆盖度阈值。默认为 NULL。
#' @param hi.perc 数值。最高覆盖度百分位数，用于去除异常值，默认为 99.9。
#' @param destrand 逻辑值。是否合并正负链数据，默认为 FALSE。
#'
#' @return 一个综合性的数据框，包含：
#'         - chr, start, end, strand: 位点信息。
#'         - context: 甲基化上下文类型。
#'         - 对每个样本，有四列：[样本名]_coverage, [样本名]_numCs, [样本名]_numTs, [样本名]_meth_percent。
#'
#' @examples
#' # 1. 创建一个包含明确样本ID的模拟 methylRawList 对象
#' mock_data1 <- data.frame(
#'   chr = "chr1", start = 1001:1005, end = 1001:1005, strand = "+",
#'   coverage = c(100, 120, 90, 150, 200),
#'   numCs = c(10, 60, 5, 75, 20),
#'   numTs = 0
#' )
#' mock_data1$numTs <- mock_data1$coverage - mock_data1$numCs
#' 
#' mock_data2 <- data.frame(
#'   chr = "chr1", start = 1001:1005, end = 1001:1005, strand = "+",
#'   coverage = c(110, 130, 80, 160, 210),
#'   numCs = c(55, 13, 40, 16, 189),
#'   numTs = 0
#' )
#' mock_data2$numTs <- mock_data2$coverage - mock_data2$numCs
#' 
#' # 假设这个对象是 CpG 数据
#' my_cpg_obj <- methylRawList(
#'    methylRaw(mock_data1, sample.id = "CK1", assembly = "hg19"),
#'    methylRaw(mock_data2, sample.id = "Treat1", assembly = "hg19"),
#'    treatment = c(0,1)
#' )
#'
#' # 2. 使用最终版函数进行处理，指定 context
#' comprehensive_cpg_table <- get_methylkit_data(
#'   methyl_obj = my_cpg_obj,
#'   context = "CpG",
#'   lo.count = 10
#' )
#'
#' # 3. 查看最终结果
#' #    注意新增了 context 列
#' print(comprehensive_cpg_table)
#'
#' # 如果您有 CHG 数据，可以这样调用：
#' # comprehensive_chg_table <- get_methylkit_data(my_chg_obj, context = "CHG")

get_methylkit_data <- function(methyl_obj,
                               context = "CpG",
                               lo.count = 5,
                               lo.perc = NULL,
                               hi.count = NULL,
                               hi.perc = 100,
                               destrand = FALSE) {
  
  # --- 1. 过滤、标准化和合并数据 ---
  message(paste("步骤 1: 正在为", context, "数据按覆盖度过滤、标准化并合并样本..."))
  myobj.final <- methylKit::filterByCoverage(methyl_obj,
                                             lo.count = lo.count,
                                             lo.perc = lo.perc,
                                             hi.count = hi.count,
                                             hi.perc = hi.perc) %>%
    methylKit::normalizeCoverage() %>%
    methylKit::unite(destrand = destrand)
  
  # --- 2. 提取数据和样本信息 ---
  message("步骤 2: 正在提取处理后的数据和样本ID...")
  united_data <- methylKit::getData(myobj.final)
  sample_names <- methylKit::getSampleID(myobj.final)
  num_samples <- length(sample_names)
  
  # --- 3. 创建基础数据框并循环添加所有信息 ---
  message("步骤 3: 正在构建包含上下文、原始计数和百分比的最终表格...")
  
  # 初始化结果数据框，包含固定的位点信息和【新增的上下文信息】
  final_df <- data.frame(
    chr = united_data$chr,
    start = united_data$start,
    end = united_data$end,
    strand = united_data$strand,
    context = context  # 在这里添加上下文列
  )
  
  # 循环遍历每个样本
  for (i in 1:num_samples) {
    current_sample_name <- sample_names[i]
    
    source_cov_col <- paste0("coverage", i)
    source_cs_col <- paste0("numCs", i)
    source_ts_col <- paste0("numTs", i)
    
    dest_cov_col <- paste0(current_sample_name, "_coverage")
    dest_cs_col <- paste0(current_sample_name, "_numCs")
    dest_ts_col <- paste0(current_sample_name, "_numTs")
    dest_pct_col <- paste0(current_sample_name, "_meth_percent")
    
    # 将原始计数数据复制到新列中
    final_df[[dest_cov_col]] <- united_data[[source_cov_col]]
    final_df[[dest_cs_col]] <- united_data[[source_cs_col]]
    final_df[[dest_ts_col]] <- united_data[[source_ts_col]]
    
    # 计算甲基化百分比
    meth_percent <- (united_data[[source_cs_col]] / united_data[[source_cov_col]]) * 100
    
    # 将计算出的百分比数据添加到新列中
    final_df[[dest_pct_col]] <- meth_percent
  }
  
  message("处理完成！最终综合性表格已生成。")
  return(final_df)
}