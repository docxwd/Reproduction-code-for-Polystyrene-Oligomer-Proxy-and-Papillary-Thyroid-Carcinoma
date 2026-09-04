# -------------------------------
# 步骤 0：设置总步骤数及进度条
# -------------------------------
# -------------------------------
# The first positional argument is the analysis directory (default: current dir).
cli_args <- commandArgs(trailingOnly = TRUE)
work_dir <- if (length(cli_args) >= 1L) cli_args[[1L]] else getwd()
# 如果目录不存在，则创建目录（包括子目录）


total_steps <- 9                                     # 定义总共9个主要步骤
progress_bar <- txtProgressBar(min = 0, max = total_steps, style = 3)  # 创建进度条对象
current_step <- 0                                    # 初始化当前步骤

# -------------------------------
# 步骤 1：设置随机种子，确保结果可重复
# -------------------------------
set.seed(12345)                                     # 固定随机种子
current_step <- current_step + 1                     # 更新进度条步骤
setTxtProgressBar(progress_bar, current_step)       # 更新进度条显示

# -------------------------------
# 步骤 2：检查并加载必要的 R 包（glmnet）
# -------------------------------
# Dependencies must be installed before execution; the script never changes
# the user's R library at runtime.
if (!requireNamespace("glmnet", quietly = TRUE)) {
  stop("Required package 'glmnet' is not installed; see environment/requirements-r.txt.")
}
suppressPackageStartupMessages(library(glmnet))
current_step <- current_step + 1                     # 更新进度条步骤
setTxtProgressBar(progress_bar, current_step)       # 更新进度条显示

# -------------------------------


if (!dir.exists(work_dir)) {                        # 检查工作目录是否存在
  dir.create(work_dir, recursive = TRUE)            # 创建目录（递归创建）
}
setwd(work_dir)                                     # 设置当前工作目录
current_step <- current_step + 1                     # 更新进度条步骤
setTxtProgressBar(progress_bar, current_step)       # 更新进度条显示

# -------------------------------
# 步骤 4：加载输入数据文件，并进行必要的判断
# -------------------------------
input_file <- if (length(cli_args) >= 2L) cli_args[[2L]] else "geneexp.csv"
# 判断数据文件是否存在
if (!file.exists(input_file)) {                    # 如果文件不存在则报错
  stop(paste("错误：输入文件", input_file, "不存在，请检查路径。"))
}
data_raw <- read.table(input_file,                   # 读取CSV文件
                       header = TRUE,                # 第一行为表头
                       sep = ",",                    # 分隔符为逗号
                       check.names = FALSE,          # 保持列名原样
                       row.names = 1)                # 第一列作为行名
current_step <- current_step + 1                     # 更新进度条步骤
setTxtProgressBar(progress_bar, current_step)       # 更新进度条显示

# -------------------------------
# 步骤 5：数据预处理 —— 转置数据矩阵
# -------------------------------
data_transposed <- t(data_raw)                      # 转置数据，使行代表样本，列代表基因
# 冗余：检查转置后的前6个基因名称
print(head(colnames(data_transposed)))             # 输出基因名称，便于调试
current_step <- current_step + 1                     # 更新进度条步骤
setTxtProgressBar(progress_bar, current_step)       # 更新进度条显示

# -------------------------------
# 步骤 6：准备Lasso回归模型所需的数据
# -------------------------------
X_matrix <- as.matrix(data_transposed)              # 将转置数据转为矩阵格式（输入变量）
# 利用正则表达式提取样本标签（假设样本名格式为 "XXX_分组信息"）
Y_labels <- gsub("(.*)\\_(.*)", "\\2", row.names(data_transposed))  # 输出分类标签
# 冗余：打印部分目标标签以供检查
print(head(Y_labels))                              # 查看前几个样本的分组信息
current_step <- current_step + 1                     # 更新进度条步骤
setTxtProgressBar(progress_bar, current_step)       # 更新进度条显示

# -------------------------------
# 步骤 7：构建Lasso回归模型及交叉验证模型
# -------------------------------
# 拟合Lasso回归模型（采用二分类，alpha=1 表示Lasso）
lasso_model <- glmnet(X_matrix, Y_labels, family = "binomial", alpha = 1)
# 手动设置模型中系数矩阵的行名，确保基因名称正确显示
dimnames(lasso_model$beta)[[1]] <- colnames(X_matrix)
# 进行10折交叉验证以选取最佳惩罚参数lambda
cv_model <- cv.glmnet(X_matrix, Y_labels, family = "binomial", alpha = 1,
                      type.measure = 'deviance', nfolds = 10)
current_step <- current_step + 1                     # 更新进度条步骤
setTxtProgressBar(progress_bar, current_step)       # 更新进度条显示

# -------------------------------
# 步骤 8：绘图——生成Lasso路径图并添加基因名称标签
# -------------------------------
# 打开PDF设备，准备保存Lasso路径图
jpeg(filename = "lasso_modified.jpeg", width = 6, height = 5.5, units = "in", res = 600, quality = 95, bg = "white")
# 绘制Lasso回归路径图，横轴为lambda值（对数刻度），自动标注部分曲线
plot(lasso_model, xvar = "lambda", label = TRUE)
# 冗余：再次提取系数矩阵以确保准确性
coef_matrix <- as.matrix(lasso_model$beta)
# 获取最后一个lambda值（图中右侧位置），并计算其对数值（因为横轴为log(lambda)）
lambda_last <- tail(lasso_model$lambda, 1)
log_lambda_last <- log(lambda_last)
# 遍历每个基因，如果在最后一个lambda处系数不为0，则手动添加标签
for (i in 1:nrow(coef_matrix)) {                   # 循环每个基因
  # 判断当前基因在最后一个lambda处的系数是否非零
  if (coef_matrix[i, length(lasso_model$lambda)] != 0) {  
    text(x = log_lambda_last,                     # 设置标签的x坐标（图右侧）
         y = coef_matrix[i, length(lasso_model$lambda)],  # 对应的系数值作为y坐标
         labels = rownames(coef_matrix)[i],       # 使用基因名称作为标签
         pos = 4,                                 # 标签显示在点的右侧
         cex = 0.7,                               # 设置标签字体大小
         col = "red")                             # 设置标签颜色
  }
}
dev.off()                                          # 关闭PDF设备，保存文件
current_step <- current_step + 1                     # 更新进度条步骤
setTxtProgressBar(progress_bar, current_step)       # 更新进度条显示

# -------------------------------
# 步骤 9：绘图——生成交叉验证误差图及输出结果表格
# -------------------------------
# 绘制交叉验证误差图，保存为PDF
jpeg(filename = "cvfit_modified.jpeg", width = 6, height = 5.5, units = "in", res = 600, quality = 95, bg = "white")
plot(cv_model)                                     # 绘制交叉验证曲线
dev.off()                                          # 关闭PDF设备

# 提取最优lambda下的所有基因系数
optimal_coefs <- coef(lasso_model, s = cv_model$lambda.min)
# 判断是否存在除截距外的非零系数
if (sum(optimal_coefs != 0) <= 1) {                 # 如果只有截距项为非零
  warning("警告：未检测到任何非零的基因系数，可能需要检查数据或模型参数。")
} else {
  # 获取所有非零系数对应的基因名称（排除截距项）
  nonzero_indices <- which(optimal_coefs != 0)      # 获取非零系数的索引
  gene_names <- rownames(optimal_coefs)[nonzero_indices]  # 提取基因名称
  gene_names <- gene_names[gene_names != "(Intercept)"]      # 排除截距项
  # 输出关键基因列表到文本文件
  write.table(gene_names,
              file = "LASSO_genes_modified.txt",
              sep = "\t",
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)
}

# 输出最优lambda下的完整基因回归系数表
full_coef_matrix <- as.matrix(coef(lasso_model, s = cv_model$lambda.min))
coef_dataframe <- data.frame(Gene = rownames(full_coef_matrix), 
                             Coefficient = full_coef_matrix[,1])
write.table(coef_dataframe,
            file = "LASSO_coefficients_modified.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)

# 输出交叉验证详细结果，包含lambda值及对应的均方误差等信息
cv_results_table <- data.frame(lambda = cv_model$lambda,
                               cvm = cv_model$cvm,
                               cvup = cv_model$cvup,
                               cvlo = cv_model$cvlo)
write.table(cv_results_table,
            file = "CV_results_modified.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)
current_step <- current_step + 1                     # 更新进度条步骤
setTxtProgressBar(progress_bar, current_step)       # 更新进度条显示

# 获取最优lambda下的回归系数
optimal_coefs <- coef(lasso_model, s = cv_model$lambda.min)

# 移除截距项并提取非零系数的基因名称
nonzero_indices <- which(optimal_coefs != 0)  # 获取非零系数的索引
gene_names <- rownames(optimal_coefs)[nonzero_indices]  # 提取基因名称
gene_names <- gene_names[gene_names != "(Intercept)"]  # 排除截距项

# 输出非零系数基因的个数
cat("筛选出的非零系数的基因个数：", length(gene_names), "\n")

# -------------------------------
# 步骤 10：清理工作与结束
# -------------------------------
close(progress_bar)                                # 关闭进度条
cat("所有步骤已完成，结果文件已保存。\n")         # 输出完成提示信息
