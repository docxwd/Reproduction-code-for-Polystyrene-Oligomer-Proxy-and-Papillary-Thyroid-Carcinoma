## =========================================================
## Clean pipeline (same outputs as your original script)
## - Input: geneexp.csv (rows=genes, cols=samples, first col=ID or rownames)
## - Input: IntersectionGenes.csv (one gene per line)
## - Output:
##   gene_analysis_results.csv
##   boxplot.<gene>.jpeg
##   ROC.<gene>.jpeg
##   density_plot.<gene>.jpeg
##   combined_<gene>.jpeg   (boxplot + ROC)
##   combined_boxplot_all_genes.jpeg
##   All_Genes_Combined_ROC.jpeg
## =========================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
  library(pROC)
  library(gridExtra)
})

# Bootstrap is intentionally unseeded; CI endpoints may vary between runs.
# Plot-only jitter seeds below do not fix the statistical bootstrap RNG.

## -------- 0) set working dir --------
cli_args <- commandArgs(trailingOnly = TRUE)
work_dir <- if (length(cli_args) >= 1L) cli_args[[1L]] else getwd()
if (!dir.exists(work_dir)) stop("Analysis directory not found: ", work_dir)
setwd(work_dir)

output_folder <- getwd()

## -------- 1) read expression matrix --------
expFile <- if (length(cli_args) >= 2L) cli_args[[2L]] else "geneexp.csv"
rt0 <- read.csv(expFile, header = TRUE, check.names = FALSE)

# 兼容两种格式：第一列叫ID / 或者已经是row.names
if ("ID" %in% colnames(rt0)) {
  rownames(rt0) <- rt0$ID
  rt <- rt0[, setdiff(colnames(rt0), "ID"), drop = FALSE]
} else {
  # 如果你本来就是row.names=1读入的，也能兼容
  rt <- rt0
}

# 清理列名空格
colnames(rt) <- trimws(colnames(rt))

## -------- 2) build group (Type) SAFELY from EACH sample name --------
# 从列名提取后缀（最后一个下划线后的部分）
suffix <- sub("^.*_", "", colnames(rt))

# 统一标签：con / tre / treat 都归一到 con / treat
Type <- ifelse(tolower(suffix) == "con", "con",
               ifelse(tolower(suffix) %in% c("tre", "treat"), "treat", NA))

if (any(is.na(Type))) {
  bad <- colnames(rt)[is.na(Type)]
  stop(
    "有些样本列名后缀不是 con/tre/treat，无法分组：\n",
    paste(bad, collapse = ", "),
    "\n请确保列名类似 xxx_con 或 xxx_tre(或xxx_treat)"
  )
}

Type <- factor(Type, levels = c("con", "treat"))

# y：给ROC用（0=con, 1=treat）
y <- ifelse(Type == "con", 0, 1)

# 最关键的防坑检查：必须只有两组
if (nlevels(droplevels(Type)) != 2) {
  stop("Type 不是两组！实际水平：", paste(levels(droplevels(Type)), collapse = ", "))
}

## -------- 3) read gene list & subset --------
geneFile <- if (length(cli_args) >= 3L) cli_args[[3L]] else "IntersectionGenes.csv"
geneRT <- read.csv(geneFile, header = FALSE, check.names = FALSE, stringsAsFactors = FALSE)

selectedGenes <- unique(trimws(as.character(geneRT[[1]])))
selectedGenes <- selectedGenes[selectedGenes != ""]

# 只保留表达矩阵里存在的基因
selectedGenes_in <- intersect(selectedGenes, rownames(rt))
if (length(selectedGenes_in) == 0) {
  stop("IntersectionGenes.csv 里的基因在 geneexp.csv 中一个都没找到。请检查基因名是否一致。")
}
rt_filtered <- rt[selectedGenes_in, , drop = FALSE]

cat("Samples:", ncol(rt_filtered), " | Genes:", nrow(rt_filtered), "\n")
print(table(Type))

## -------- 4) results table --------
result_df <- data.frame(
  Gene = character(),
  P_Value = numeric(),
  Mean_Con = numeric(),
  Mean_Treat = numeric(),
  AUC = numeric(),
  AUC_Lower_CI = numeric(),
  AUC_Upper_CI = numeric(),
  stringsAsFactors = FALSE
)

## -------- 5) progress bar --------
total_genes <- nrow(rt_filtered)
cat("开始基因差异分析...\n")
pb <- txtProgressBar(min = 0, max = total_genes, style = 3)

## -------- 6) per-gene loop: t-test + plots + ROC --------
roc_list <- list()
auc_leg <- character()
my_cols <- c("red", "deepskyblue", "forestgreen", "orange",
             "purple", "gray40", "black", "magenta", "gold", "brown")
if (length(my_cols) < total_genes) my_cols <- rep(my_cols, length.out = total_genes)

for (g_idx in seq_len(total_genes)) {
  gene <- rownames(rt_filtered)[g_idx]
  expr <- as.numeric(rt_filtered[gene, ])
  
  rt1 <- data.frame(expression = expr, Type = Type)
  
  ## ---- t-test (robust) ----
  p_value <- NA_real_
  mean_con <- mean(rt1$expression[rt1$Type == "con"], na.rm = TRUE)
  mean_treat <- mean(rt1$expression[rt1$Type == "treat"], na.rm = TRUE)
  
  # 只有当两组都有足够样本时才做t检验
  if (sum(rt1$Type == "con") >= 2 && sum(rt1$Type == "treat") >= 2) {
    tt <- tryCatch(t.test(expression ~ Type, data = rt1), error = function(e) NULL)
    if (!is.null(tt)) p_value <- tt$p.value
  }
  
  ## ---- ROC + CI (robust) ----
  auc_val <- NA_real_
  ci_low <- NA_real_
  ci_high <- NA_real_
  
  roc1 <- tryCatch(
    pROC::roc(response = Type, predictor = expr, levels = c("con", "treat"), quiet = TRUE),
    error = function(e) NULL
  )
  
  if (!is.null(roc1)) {
    auc_val <- as.numeric(pROC::auc(roc1))
    # AUC方向自动修正（保证 >=0.5）
    if (!is.na(auc_val) && auc_val < 0.5) auc_val <- 1 - auc_val
    
    ci1 <- tryCatch(pROC::ci.auc(roc1, method = "bootstrap", boot.n = 1000), error = function(e) NULL)
    if (!is.null(ci1)) {
      ciVec <- as.numeric(ci1)
      ci_low <- ciVec[1]; ci_high <- ciVec[3]
    }
  }
  
  ## ---- save stats ----
  result_df <- rbind(result_df, data.frame(
    Gene = gene,
    P_Value = p_value,
    Mean_Con = mean_con,
    Mean_Treat = mean_treat,
    AUC = auc_val,
    AUC_Lower_CI = ci_low,
    AUC_Upper_CI = ci_high,
    stringsAsFactors = FALSE
  ))
  
  ## ---- boxplot (with stat_compare_means) ----
  colors <- c("con" = "#FF6347", "treat" = "#4682B4")
  group_labels <- c("con" = "Con", "treat" = "Tumor")
  
  bp <- ggplot(rt1, aes(x = Type, y = expression, fill = Type)) +
    geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.75) +
    geom_jitter(
      color = "black",
      size = 1.8,
      position = position_jitter(width = 0.18, height = 0, seed = 20260127)
    ) +
    scale_x_discrete(labels = group_labels) +
    scale_fill_manual(values = colors, labels = group_labels) +
    theme_minimal(base_size = 14) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14),
          legend.position = "none") +
    labs(y = paste(gene, "expression")) +
    stat_compare_means(method = "t.test")
  
  jpeg(filename = paste0("boxplot.", gene, ".jpeg"), width = 3.6, height = 4.6, units = "in", res = 600, quality = 95, bg = "white")
  print(bp)
  dev.off()
  
  ## ---- density plot ----
  dp <- ggplot(rt1, aes(x = expression, fill = Type)) +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = colors, labels = group_labels) +
    theme_minimal(base_size = 14) +
    labs(title = paste(gene, "Density Plot"), x = "Expression", y = "Density")
  
  jpeg(filename = paste0("density_plot.", gene, ".jpeg"), width = 6, height = 6, units = "in", res = 600, quality = 95, bg = "white")
  print(dp)
  dev.off()
  
  ## ---- ROC plot ----
  if (!is.null(roc1)) {
    jpeg(filename = paste0("ROC.", gene, ".jpeg"), width = 5, height = 5, units = "in", res = 600, quality = 95, bg = "white")
    plot(roc1, print.auc = TRUE, col = "orange", legacy.axes = TRUE, main = gene, lwd = 2)
    abline(a = 0, b = 1, lty = 2, col = "gray70")
    if (is.finite(ci_low) && is.finite(ci_high)) {
      text(0.60, 0.15,
           paste0("95% CI: ", sprintf("%.3f", ci_low), " - ", sprintf("%.3f", ci_high)),
           col = "orange")
    }
    dev.off()
  }
  
  ## ---- combined (boxplot + ROC) ----
  if (!is.null(roc1)) {
    roc_plot <- pROC::ggroc(roc1, legacy.axes = TRUE) +
      ggtitle(paste0(gene, "  AUC=", sprintf("%.3f", auc_val))) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "none")
    
    comb <- gridExtra::grid.arrange(bp, roc_plot, ncol = 2, top = paste("Gene:", gene))
    ggsave(filename = file.path(output_folder, paste0("combined_", gene, ".jpeg")),
           plot = comb, width = 10, height = 6, dpi = 600, device = "jpeg", quality = 95, bg = "white")
  }
  
  ## ---- store for combined ROC ----
  if (!is.null(roc1)) {
    roc_list[[gene]] <- roc1
    auc_leg <- c(auc_leg, paste0(gene, "  AUC=", sprintf("%.3f", auc_val)))
  }
  
  setTxtProgressBar(pb, g_idx)
}

close(pb)

## -------- 7) save summary table --------
write.csv(result_df, file = "gene_analysis_results.csv", row.names = FALSE)
cat("\n差异分析结果已保存：gene_analysis_results.csv\n")

## -------- 8) combined boxplot for all genes --------
combined_data <- do.call(rbind, lapply(rownames(rt_filtered), function(gene) {
  data.frame(
    expression = as.numeric(rt_filtered[gene, ]),
    Type = Type,
    Gene = gene,
    stringsAsFactors = FALSE
  )
}))

cb <- ggplot(combined_data, aes(x = Gene, y = expression, fill = Type)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.75) +
  geom_point(
    color = "black",
    size = 1.0,
    alpha = 0.75,
    position = position_jitterdodge(
      jitter.width = 0.12,
      jitter.height = 0,
      dodge.width = 0.65,
      seed = 20260127
    ),
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c("con" = "#FF6347", "treat" = "#4682B4"),
                    labels = c("Con", "Tumor")) +
  theme_minimal(base_size = 14) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(y = "Expression", title = "Differential Expression (Selected Genes)") +
  stat_compare_means(method = "t.test", label = "p.signif") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave("combined_boxplot_all_genes.jpeg", plot = cb, width = 9, height = 6, dpi = 600, device = "jpeg", quality = 95, bg = "white")
cat("合并箱线图已保存：combined_boxplot_all_genes.jpeg\n")

## -------- 9) All genes combined ROC (overlay) --------
if (length(roc_list) > 0) {
  jpeg(filename = "All_Genes_Combined_ROC.jpeg", width = 6, height = 6, units = "in", res = 600, quality = 95, bg = "white")
  par(mar = c(5, 6, 4, 2) + 0.1, cex = 1.1)
  
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
       xlab = "1 - Specificity", ylab = "Sensitivity",
       main = "All Genes Combined ROC")
  abline(0, 1, lty = 2, col = "gray70", lwd = 2)
  
  genes_to_plot <- names(roc_list)
  for (i in seq_along(genes_to_plot)) {
    gene <- genes_to_plot[i]
    roc1 <- roc_list[[gene]]
    lines(1 - roc1$specificities, roc1$sensitivities, col = my_cols[i], lwd = 2.5)
  }
  
  legend("bottomright",
         legend = auc_leg,
         col = my_cols[seq_along(genes_to_plot)],
         lwd = 2.5, bty = "n", cex = 0.9)
  
  dev.off()
  cat("所有基因合并ROC图已保存：All_Genes_Combined_ROC.jpeg\n")
} else {
  cat("提示：没有成功生成任何ROC对象，跳过 All_Genes_Combined_ROC.jpeg\n")
}

cat("\n全部完成 ✅\n")
