#' DNA序列编码函数
#' @param dna_strings 一个字符向量，包含DNA序列
#' @return 一个数据框，包含编码后的序列
dna_encoding <- function(dna_strings) {
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("base", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, function(x) factor(x, levels = c("A", "T", "C", "G")))
  return(seq_df)
}

#' 多个样本的m6A预测
#' @param ml_fit 训练好的随机森林模型
#' @param feature_df 包含特征的数据框
#' @param positive_threshold 判断为正的阈值，默认0.5
#' @return 增加预测结果的数据框
#' @import randomForest
#' @export
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5) {
  stopifnot(all(c("gc_content", "evolutionary_conservation", "RNA_type", "RNA_region",
                  "exon_length", "distance_to_junction", "DNA_5mer") %in% colnames(feature_df)))

  feature_df$RNA_type <- factor(feature_df$RNA_type)
  feature_df$RNA_region <- factor(feature_df$RNA_region)

  inferred_prob <- predict(ml_fit, newdata = feature_df, type = "prob")[, 2]
  feature_df$predicted_m6A_prob <- inferred_prob
  feature_df$predicted_m6A_status <- ifelse(inferred_prob > positive_threshold, "Positive", "Negative")

  return(feature_df)
}

#' 单个样本的m6A预测
#' @param ml_fit 训练好的模型
#' @param gc_content GC含量
#' @param RNA_type RNA类型
#' @param RNA_region RNA区域
#' @param exon_length 外显子长度
#' @param distance_to_junction 到连接点的距离
#' @param evolutionary_conservation 保守性得分
#' @param DNA_5mer DNA 5碱基序列
#' @param positive_threshold 阈值，默认0.5
#' @return 一个命名向量，包含预测概率和状态
#' @examples
#' \dontrun{
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
#' prediction_single(ml_fit, gc_content = 0.6, RNA_type = "mRNA", RNA_region = "CDS",
#'                   exon_length = 12, distance_to_junction = 5,
#'                   evolutionary_conservation = 0.8, DNA_5mer = "ATCGAT",
#'                   positive_threshold = 0.5)
#' }
#' @export
#' @import randomForest
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region,
                              exon_length, distance_to_junction,
                              evolutionary_conservation, DNA_5mer,
                              positive_threshold = 0.5) {
  input_df <- data.frame(
    gc_content = gc_content,
    RNA_type = RNA_type,
    RNA_region = RNA_region,
    exon_length = exon_length,
    distance_to_junction = distance_to_junction,
    evolutionary_conservation = evolutionary_conservation,
    DNA_5mer = DNA_5mer
  )
  result <- prediction_multiple(ml_fit, input_df, positive_threshold)
  return(c(predicted_m6A_prob = result$predicted_m6A_prob[1],
           predicted_m6A_status = result$predicted_m6A_status[1]))
}
