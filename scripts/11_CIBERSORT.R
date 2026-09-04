# Dependencies must be installed before execution; see
# environment/requirements-r.txt. This script never installs packages.
required_packages <- c("e1071", "preprocessCore")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages) > 0L) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

# First positional argument: analysis directory (default: current directory).
# Second: mixture matrix; third: reference signature matrix.
cli_args <- commandArgs(trailingOnly = TRUE)
wd_path <- if (length(cli_args) >= 1L) cli_args[[1L]] else getwd()
if (dir.exists(wd_path)) {                     # 如果工作目录存在则
  setwd(wd_path)                              # 设置工作目录
} else {                                       # 否则
  stop("工作目录不存在！")                      # 报错并终止程序
}

# 定义输入文件路径
input_file <- if (length(cli_args) >= 2L) cli_args[[2L]] else "Sample Type Matrix.csv"
reference_file <- if (length(cli_args) >= 3L) cli_args[[3L]] else "refer.txt"
set.seed(20260825)

# 定义核心算法函数，用于对参考矩阵和目标向量进行SVR回归分析（尝试不同的nu值）
coreAlgorithm <- function(refer_matrix, target_vector) {  # 函数接收参考矩阵和目标向量
  num_steps <- 3                               # 定义尝试的步数（对应不同的nu值）
  
  # 定义内部函数，根据索引选择不同的nu参数并训练SVM模型
  runSVM <- function(step_index) {             # 内部函数，接收步长索引
    if (step_index == 1) {                      # 如果索引为1
      nu_param <- 0.25                         # 设置nu参数为0.25
    } else if (step_index == 2) {               # 如果索引为2
      nu_param <- 0.5                          # 设置nu参数为0.5
    } else if (step_index == 3) {               # 如果索引为3
      nu_param <- 0.75                         # 设置nu参数为0.75
    } else {                                    # 否则（冗余判断）
      nu_param <- 0.5                          # 默认设置nu参数为0.5
    }
    # 训练SVM模型，使用线性核且不进行数据缩放
    svm_model <- svm(refer_matrix, target_vector, type="nu-regression", kernel="linear", nu=nu_param, scale=FALSE)
    return(svm_model)                          # 返回训练好的模型
  }  # 结束内部函数
  
  # 根据系统类型选择并行计算方式（Windows系统不支持多核并行）
  if (Sys.info()['sysname'] == 'Windows') {     # 如果操作系统为Windows
    models_list <- mclapply(1:num_steps, runSVM, mc.cores=1)  # 使用单核运行
  } else {                                      # 否则（Linux或Mac）
    models_list <- mclapply(1:num_steps, runSVM, mc.cores=num_steps)  # 使用多个核心并行计算
  }
  
  # 初始化存储每个模型的均方根误差（RMSE）和相关系数的向量
  rmse_vec <- rep(0, num_steps)                # 初始化RMSE向量
  corr_vec <- rep(0, num_steps)                # 初始化相关系数向量
  
  # 对每个训练好的模型计算加权系数及评价指标
  for (i in 1:num_steps) {                     # 循环遍历每个模型
    current_model <- models_list[[i]]          # 获取当前模型
    # 计算模型权重（系数的转置乘以支持向量）
    weight_vals <- t(current_model$coefs) %*% current_model$SV  
    weight_vals[which(weight_vals < 0)] <- 0    # 将负权重置为0
    norm_weights <- weight_vals / sum(weight_vals)  # 归一化权重
    # 对参考矩阵每列按对应权重进行加权
    weighted_refer <- sweep(refer_matrix, MARGIN=2, norm_weights, '*')
    # 按行求和，得到估计的目标向量
    est_target <- apply(weighted_refer, 1, sum)
    # 计算当前模型的RMSE
    rmse_vec[i] <- sqrt(mean((est_target - target_vector)^2))
    # 计算当前模型的皮尔逊相关系数
    corr_vec[i] <- cor(est_target, target_vector)
  }
  
  # 选择RMSE最小的模型作为最佳模型
  best_idx <- which.min(rmse_vec)              # 获取RMSE最小模型的索引
  best_model <- models_list[[best_idx]]          # 选择最佳模型
  
  # 重新计算最佳模型的归一化系数
  best_weights <- t(best_model$coefs) %*% best_model$SV  
  best_weights[which(best_weights < 0)] <- 0     # 将负权重置为0
  final_weights <- best_weights / sum(best_weights)  # 归一化最终权重
  
  # 存储最佳模型的RMSE和相关系数
  best_rmse <- rmse_vec[best_idx]              # 最佳RMSE
  best_corr <- corr_vec[best_idx]              # 最佳相关系数
  
  # 将最终结果存入列表并返回
  result_list <- list("final_weights" = final_weights, "best_rmse" = best_rmse, "best_corr" = best_corr)
  return(result_list)                         # 返回结果列表
}  # 结束coreAlgorithm函数

# 定义置换检验函数，用于生成经验零分布
doPermutation <- function(num_perm, refer_matrix, mix_matrix) {  # 函数接收置换次数、参考矩阵和混合矩阵
  perm_counter <- 1                          # 初始化置换计数器
  target_list <- as.list(data.matrix(mix_matrix))  # 将混合矩阵转换为列表，方便采样
  corr_distribution <- numeric()             # 初始化存储相关系数的向量
  
  # 初始化置换进度条
  perm_progress <- txtProgressBar(min = 0, max = num_perm, style = 3)
  
  # 进行置换循环
  while (perm_counter <= num_perm) {          # 当置换计数器小于等于置换次数时
    # 随机采样混合数据（允许重复采样），生成新的目标向量
    permuted_target <- as.numeric(target_list[sample(length(target_list), nrow(refer_matrix), replace = TRUE)])
    # 判断标准差是否为0，避免除零错误
    if (sd(permuted_target) == 0) {           # 如果标准差为0
      std_target <- permuted_target          # 保持原样
    } else {                                  
      std_target <- (permuted_target - mean(permuted_target)) / sd(permuted_target)  # 标准化目标向量
    }
    # 调用核心算法函数处理置换数据
    perm_result <- coreAlgorithm(refer_matrix, std_target)
    current_corr <- perm_result$best_corr     # 提取当前相关系数
    # 将当前相关系数追加到分布向量中
    corr_distribution <- c(corr_distribution, current_corr)
    # 更新进度条显示
    setTxtProgressBar(perm_progress, perm_counter)
    perm_counter <- perm_counter + 1         # 增加计数器
  }
  
  close(perm_progress)                        # 关闭置换进度条
  sorted_corr <- sort(corr_distribution)       # 将相关系数分布进行排序
  return(list("null_distribution" = sorted_corr))  # 返回经验零分布列表
}  # 结束doPermutation函数

# 定义主函数，用于执行CIBERSORT分析
CIBERSORTModified <- function(refer_file, mix_file, num_perm = 0, do_QN = TRUE) {  # 主函数接收参考文件、混合文件、置换次数和是否进行量化归一化的标志
  library(e1071)                              # 加载e1071包
  library(parallel)                           # 加载parallel包，用于并行计算
  library(preprocessCore)                     # 加载preprocessCore包，用于量化归一化
  
  # 检查参考文件是否存在
  if (!file.exists(refer_file)) {              # 如果参考文件不存在
    stop("参考文件未找到！")                    # 报错并停止程序
  }
  # 检查混合文件是否存在
  if (!file.exists(mix_file)) {              # 如果混合文件不存在
    stop("混合文件未找到！")                    # 报错并停止程序
  }
  
  # 读取参考矩阵数据（假设为制表符分隔）
  refer_data <- read.table(refer_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  # 读取混合矩阵数据（假设为逗号分隔）
  mix_data <- read.table(mix_file, header = TRUE, sep = ",", row.names = 1, check.names = FALSE)
  
  # 将数据框转换为矩阵
  refer_data <- data.matrix(refer_data)
  mix_data <- data.matrix(mix_data)
  
  # 对参考数据和混合数据按照行名进行排序
  refer_data <- refer_data[order(rownames(refer_data)), ]
  mix_data <- mix_data[order(rownames(mix_data)), ]
  
  # 保存置换次数到变量中
  perm_count <- num_perm
  
  # 如果混合数据最大值小于50，则进行反对数转换
  if (max(mix_data) < 50) {
    mix_data <- 2^mix_data
  }
  
  # 如果启用量化归一化，则对混合数据进行归一化处理
  if (do_QN == TRUE) {
    orig_colnames <- colnames(mix_data)      # 保存原始列名
    orig_rownames <- rownames(mix_data)      # 保存原始行名
    mix_data <- normalize.quantiles(mix_data)  # 进行量化归一化
    colnames(mix_data) <- orig_colnames       # 恢复列名
    rownames(mix_data) <- orig_rownames       # 恢复行名
  }
  
  # 找出参考数据和混合数据中共同存在的基因
  refer_genes <- rownames(refer_data)
  mix_genes <- rownames(mix_data)
  common_genes <- intersect(refer_genes, mix_genes)  # 取交集
  if (length(common_genes) == 0) {             # 如果没有共同基因，则终止程序
    stop("参考数据与混合数据之间没有共同基因！")
  }
  
  # 子集化参考矩阵和混合矩阵，只保留共同基因
  refer_data <- refer_data[common_genes, , drop = FALSE]
  mix_data <- mix_data[common_genes, , drop = FALSE]
  
  # 对参考矩阵进行标准化（z-score标准化）
  refer_data <- (refer_data - mean(refer_data)) / sd(as.vector(refer_data))
  
  # 如果置换次数大于0，则计算经验零分布
  null_dist <- NULL                         # 初始化空的零分布
  if (perm_count > 0) {
    perm_results <- doPermutation(perm_count, refer_data, mix_data)
    null_dist <- perm_results$null_distribution  # 获取置换得到的零分布
  }
  
  # 定义输出表头
  output_header <- c("Mixture", colnames(refer_data), "P-value", "Correlation", "RMSE")
  output_matrix <- NULL                     # 初始化输出矩阵
  
  # 初始化混合样本处理进度条
  total_samples <- ncol(mix_data)           # 混合样本总数
  sample_progress <- txtProgressBar(min = 0, max = total_samples, style = 3)
  
  # 循环处理每个混合样本
  for (sample_idx in 1:total_samples) {
    curr_sample <- mix_data[, sample_idx]   # 提取当前混合样本
    # 判断当前样本标准差是否为0，避免除零错误
    if (sd(curr_sample) == 0) {
      std_sample <- curr_sample             # 如果标准差为0，则保持原样
    } else {
      std_sample <- (curr_sample - mean(curr_sample)) / sd(curr_sample)  # 否则进行标准化
    }
    # 运行核心算法处理当前样本
    alg_result <- coreAlgorithm(refer_data, std_sample)
    weight_vec <- alg_result$final_weights  # 提取最终权重
    corr_val <- alg_result$best_corr          # 提取相关系数
    rmse_val <- alg_result$best_rmse          # 提取RMSE值
    
    # 计算经验p值（当置换结果可用时）
    p_val <- NA                             # 初始化p值为NA
    if (!is.null(null_dist)) {              # 如果零分布存在
      p_val <- 1 - (which.min(abs(null_dist - corr_val)) / length(null_dist))  # 计算经验p值
    }
    
    # 整理当前样本结果，组合样本名、权重、p值、相关系数和RMSE
    curr_result <- c(colnames(mix_data)[sample_idx], weight_vec, p_val, corr_val, rmse_val)
    if (is.null(output_matrix)) {           # 如果输出矩阵为空（第一次循环）
      output_matrix <- curr_result          # 则初始化输出矩阵
    } else {                                # 否则
      output_matrix <- rbind(output_matrix, curr_result)  # 追加当前结果到输出矩阵
    }
    
    # 更新混合样本进度条
    setTxtProgressBar(sample_progress, sample_idx)
  }
  
  close(sample_progress)                    # 关闭样本处理进度条
  
  # 将输出表头与结果合并后写入CSV文件
  write.table(rbind(output_header, output_matrix), file = "CIBERSORT-Results.csv", sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # 将结果转化为数值矩阵（去掉样本名称和表头），便于后续筛选
  result_numeric <- rbind(output_header, output_matrix)
  result_numeric <- result_numeric[, -1]    # 移除第一列（样本名称）
  result_numeric <- result_numeric[-1, ]     # 移除表头行
  result_numeric <- matrix(as.numeric(unlist(result_numeric)), nrow = nrow(output_matrix))
  rownames(result_numeric) <- colnames(mix_data)  # 设置行名为混合样本名
  colnames(result_numeric) <- c(colnames(refer_data), "P-value", "Correlation", "RMSE")  # 设置列名
  
  return(result_numeric)                    # 返回数值结果矩阵
}  # 结束CIBERSORTModified函数

# 运行修改后的CIBERSORT分析，并将结果存储到result_table中
result_table <- CIBERSORTModified(reference_file, input_file, num_perm = 1000, do_QN = TRUE)

# 筛选出p-value小于0.05的结果
filtered_table <- result_table[result_table[, "P-value"] < 0.05, ]
# 提取权重信息（排除最后三列：P-value、Correlation、RMSE）
weight_only <- as.matrix(filtered_table[, 1:(ncol(filtered_table) - 3)])
# 将权重矩阵的列名作为首行添加到结果中
final_output <- rbind(id = colnames(weight_only), weight_only)
# 将最终结果写入CSV文件，注意不输出引号和列名
write.table(final_output, file = "CIBERSORT-Results.csv", sep = ",", quote = FALSE, col.names = FALSE)
