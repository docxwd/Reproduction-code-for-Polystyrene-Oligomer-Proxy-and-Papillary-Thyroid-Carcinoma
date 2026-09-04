#---------------------#
# 1. 包加载函数 #
# 分析包（Bioconductor）
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ComplexHeatmap)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ggpubr)
#------------------------------#
# 2. 初始化参数/环境与进度条    #
#------------------------------#
cli_args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1L]]) else "scripts/04_GO_enrichment.R"
package_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
geneFile <- normalizePath(
  if (length(cli_args) >= 2L) cli_args[[2L]] else file.path(package_root, "data", "processed", "shared_targets_75.csv"),
  winslash = "/",
  mustWork = TRUE
)
wdPath <- if (length(cli_args) >= 1L) cli_args[[1L]] else file.path(package_root, "reproduced_outputs", "GO_submitted_75")
wdPath <- normalizePath(wdPath, winslash = "/", mustWork = FALSE)
if (!dir.exists(wdPath)) dir.create(wdPath, recursive = TRUE)
setwd(wdPath)

updateProgress <- function(pb, value, message) {
  setTxtProgressBar(pb, value)
  cat(message, "\n")
}
totalSteps <- 9
pb <- txtProgressBar(min = 0, max = totalSteps, style = 3)

pvalThreshold <- 0.05
padjThreshold <- 0.05
colorParameter <- if (padjThreshold > 0.05) "pvalue" else "p.adjust"
updateProgress(pb, 1, "第1步：设置参数和工作环境完成")

#------------------------------#
# 3. 读取基因数据               #
#------------------------------#
if (!file.exists(geneFile)) stop("错误：基因文件 不存在")
geneData <- read.csv(geneFile, header = TRUE, check.names = FALSE)
updateProgress(pb, 2, "第2步：读取基因数据完成")

geneSymbols <- unique(as.vector(if ("Gene" %in% names(geneData)) geneData$Gene else geneData[, 1]))
if (length(geneSymbols) == 0) stop("错误：未找到有效的基因符号")

#------------------------------#
# 4. 基因符号转Entrez ID       #
#------------------------------#
entrezMapping <- mget(geneSymbols, org.Hs.egSYMBOL2EG, ifnotfound = NA)
entrezIDs <- as.character(entrezMapping)
validGenes <- entrezIDs[!is.na(entrezIDs) & entrezIDs!="NA"]
if (length(validGenes) == 0) stop("错误：无有效的 Entrez 基因 ID")
updateProgress(pb, 3, "第3步：基因符号转换为 EntrezID完成")

#------------------------------#
# 5. GO富集与过滤              #
#------------------------------#
goAnalysis <- enrichGO(
  gene = validGenes,
  OrgDb = org.Hs.eg.db,
  pvalueCutoff = 1, qvalueCutoff = 1,
  ont = "all",
  readable = TRUE
)
goResult <- as.data.frame(goAnalysis)
if (nrow(goResult)==0) stop("警告：未检测到任何富集结果")
filteredGO <- goResult[goResult$pvalue < pvalThreshold & goResult$p.adjust < padjThreshold, ]
updateProgress(pb, 4, "第4步：GO 富集分析及结果过滤完成")

#------------------------------#
# 6. 结果输出                  #
#------------------------------#
outputFile <- "GO_results.txt"
write.table(filteredGO, file = outputFile, sep = "\t", quote = FALSE, row.names = FALSE)
updateProgress(pb, 5, "第5步：富集结果写入文件完成")

#------------------------------#
# 7. 柱状图/气泡图 PDF绘制     #
#------------------------------#
# 柱状图
jpeg("GO_barplot.jpeg", width=8, height=10, units="in", res=600, quality=95, bg="white")
barPlot <- barplot(
  goAnalysis,
  drop=TRUE, showCategory=10, label_format=50,
  split="ONTOLOGY", color=colorParameter
) +
  facet_grid(ONTOLOGY ~ ., scale = 'free') +
  scale_fill_gradientn(colors = c("#FF6666", "#FFB266", "#FFFF99", "#99FF99", "#6666FF", "#7F52A0", "#B266FF"))
print(barPlot)
dev.off()

# 气泡图
jpeg("GO_bubble.jpeg", width=8, height=10, units="in", res=600, quality=95, bg="white")
bubblePlot <- dotplot(
  goAnalysis,
  showCategory = 10,
  orderBy = "GeneRatio",
  label_format = 50,
  split = "ONTOLOGY",
  color = colorParameter
) +
  facet_grid(ONTOLOGY ~ ., scale = 'free') +
  scale_color_gradientn(colors = c("#FFB266", "#FFFF99", "#99FF99", "#6666FF", "#7F52A0", "#B266FF"))
print(bubblePlot)
dev.off()
updateProgress(pb, 6, "第6步：柱状图和气泡图生成完成")

#------------------------------#
# 8. 分组条形图 ggbarplot       #
#------------------------------#
topGO <- filteredGO %>% group_by(ONTOLOGY) %>% slice_head(n = 10)
jpeg("GO_grouped_barplot.jpeg", width=11, height=8, units="in", res=600, quality=95, bg="white")
groupBarPlot <- ggbarplot(
  topGO,
  x="Description", y="Count", fill="ONTOLOGY", color="white",
  xlab="", palette="aaas",
  legend="right", sort.val="desc", sort.by.groups=TRUE,
  position=position_dodge(0.9)
) +
  rotate_x_text(75) +
  theme(panel.background = element_blank(),
        axis.text.x = element_text(size=10, color="black")) +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_discrete(expand=c(0,0)) +
  geom_text(
    aes(label=Count),
    position=position_dodge(0.9), vjust=-0.3, size=3
  )
print(groupBarPlot)
dev.off()
updateProgress(pb, 7, "第7步：分组条形图生成完成")

#------------------------------#
# 9. chord弦图 PDF绘制         #
#------------------------------#
jpeg("GO_chord_diagram.jpeg", width=12, height=12, units="in", res=600, quality=95, bg="white")
go <- read.delim("GO_results.txt", header=TRUE, stringsAsFactors=FALSE)
top_terms <- go %>%
  group_by(ONTOLOGY) %>%
  slice_min(order_by = p.adjust, n = 10) %>%
  ungroup()
insert_linebreak <- function(text, line_length=35) {
  if(nchar(text) <= line_length) return(text)
  paste(strwrap(text, width=line_length), collapse = "\n")
}
top_terms$Description_new <- sapply(top_terms$Description, insert_linebreak)
mat <- table(
  factor(top_terms$ONTOLOGY, levels=c("BP","CC","MF")),
  factor(top_terms$Description_new, levels=unique(top_terms$Description_new))
)
n_ont <- 3
n_term <- ncol(mat)
gap.deg <- c(rep(1, n_ont-1), 10, rep(1, n_term-1), 10)
grid.col <- c(BP="#E69F00", CC="#56B4E9", MF="#009E73",
              setNames(rep("#BBBBBB", n_term), colnames(mat)))
circos.clear()
circos.par(gap.degree = gap.deg, start.degree = 90)
chordDiagram(
  mat,
  grid.col = grid.col,
  transparency = 0.4,
  annotationTrack = c("", "grid"),
  preAllocateTracks = list(track.height = 0.2)
)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  sector.name <- get.cell.meta.data("sector.index")
  xlim <- get.cell.meta.data("xlim")
  ylim <- get.cell.meta.data("ylim")
  circos.text(mean(xlim), ylim[1] + 0.1, sector.name,
              facing="clockwise", niceFacing=TRUE,
              adj=c(0,0.5), cex=0.60)
}, bg.border=NA)
title("GO Ontology")
circos.clear()
dev.off()
updateProgress(pb, 8, "第8步：chord弦图生成完成")

# 载入必要的包
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# 读取GO富集分析的结果（你的filteredGO）
filteredGO <- read.table("GO_results.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)

# TRUE则全部三类（BP, CC, MF），如果只要其中某几类，改为只留需要的类别名称
filteredGO_sub <- filteredGO %>% filter(ONTOLOGY %in% c("BP", "CC", "MF"))

# 分别对每个大类取Top10最显著的GO term
topGO <- filteredGO_sub %>%
  group_by(ONTOLOGY) %>%
  slice_min(order_by = p.adjust, n = 10) %>%
  ungroup() %>%
  arrange(ONTOLOGY, p.adjust)

# 高级：保证每一类内，Description是从最显著到最不显著排序，显示时不乱
topGO <- topGO %>% 
  group_by(ONTOLOGY) %>%
  mutate(Description = factor(Description, levels = rev(unique(Description)))) %>%
  ungroup()

# GeneRatio转数值（主要用于气泡图横轴）
topGO$GeneRatio_num <- sapply(topGO$GeneRatio, function(x) {
  sp <- unlist(strsplit(as.character(x), "/"))
  as.numeric(sp[1]) / as.numeric(sp[2])
})

# GO富集条形图（分面，和KEGG一样）
bar_colors <- brewer.pal(7, "YlOrRd")
p_bar <- ggplot(topGO, aes(x = Description, y = -log10(p.adjust), fill = -log10(p.adjust))) +
  geom_bar(stat = "identity", width = 0.8) +
  coord_flip() +
  facet_wrap(~ ONTOLOGY, scales = "free_y", ncol = 1) +
  scale_fill_gradientn(colors = bar_colors) +
  labs(x = "GO Term", y = "-log10(Adjusted p-value)", title = "GO Enrichment Analysis (BP/CC/MF)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", color = "darkblue"),
    strip.text = element_text(size = 15, face = "bold", color = "black"),
    axis.text.y = element_text(size = 11)
  )
ggsave("GO_barplot_custom.jpeg", p_bar, width = 12, height = 12, dpi = 600, device = "jpeg", quality = 95, bg = "white")

# GO富集气泡图（分面，和KEGG一样）
dot_colors <- brewer.pal(7, "Spectral")
p_dot <- ggplot(topGO, aes(x = GeneRatio_num, y = Description, size = Count, color = p.adjust)) +
  geom_point(alpha = 0.8) +
  facet_wrap(~ ONTOLOGY, scales = "free_y", ncol = 1) +
  scale_color_gradientn(colors = dot_colors) +
  labs(x = "Gene Ratio", y = "GO Term", title = "GO Dotplot (BP/CC/MF)") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", color = "darkred"),
    strip.text = element_text(size = 15, face = "bold", color = "black"),
    axis.text.y = element_text(size = 11)
  )
ggsave("GO_dotplot_custom.jpeg", p_dot, width = 12, height = 12, dpi = 600, device = "jpeg", quality = 95, bg = "white")

# ---------------------------
# GO富集气泡图/条形图 分面 分别显示BP/CC/MF 通路名字加粗
# ---------------------------


# 2. 读取数据（用你自己的GO富集结果路径）
filteredGO <- read.table("GO_results.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)

# 3. 只保留BP,CC,MF三类（如想只看部分，改filter内容即可）
filteredGO_sub <- filteredGO %>% filter(ONTOLOGY %in% c("BP", "CC", "MF"))

# 4. 每类top10最显著GO term
topGO <- filteredGO_sub %>%
  group_by(ONTOLOGY) %>%
  slice_min(order_by = p.adjust, n = 10) %>%
  ungroup() %>%
  arrange(ONTOLOGY, p.adjust)

# 5. 保证每个类别内部按显著性排序
topGO <- topGO %>% 
  group_by(ONTOLOGY) %>%
  mutate(Description = factor(Description, levels = rev(unique(Description)))) %>%
  ungroup()

# 6. GeneRatio数值化
topGO$GeneRatio_num <- sapply(topGO$GeneRatio, function(x) {
  sp <- unlist(strsplit(as.character(x), "/"))
  as.numeric(sp[1]) / as.numeric(sp[2])
})

# 7. 条形图 (KEGG style, 通路名字加粗)
bar_colors <- brewer.pal(7, "YlOrRd")
p_bar <- ggplot(topGO, aes(x = Description, y = -log10(p.adjust), fill = -log10(p.adjust))) +
  geom_bar(stat = "identity", width = 0.8) +
  coord_flip() +
  facet_wrap(~ ONTOLOGY, scales = "free_y", ncol = 1) +
  scale_fill_gradientn(colors = bar_colors) +
  labs(x = "GO Term", y = "-log10(adjusted p-value)", title = "GO Enrichment Analysis (BP/CC/MF)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", color = "darkblue"),
    strip.text = element_text(size = 15, face = "bold", color = "black"),
    axis.text.y = element_text(size = 11, face = "bold") # 通路名字加粗
  )
ggsave("GO_barplot_custom.jpeg", p_bar, width = 12, height = 12, dpi = 600, device = "jpeg", quality = 95, bg = "white")

# 8. 气泡图 (KEGG style, 通路名字加粗)
dot_colors <- brewer.pal(7, "Spectral")
p_dot <- ggplot(topGO, aes(x = GeneRatio_num, y = Description, size = Count, color = p.adjust)) +
  geom_point(alpha = 0.8) +
  facet_wrap(~ ONTOLOGY, scales = "free_y", ncol = 1) +
  scale_color_gradientn(colors = dot_colors) +
  labs(x = "Gene Ratio", y = "GO Term", title = "GO Dotplot (BP/CC/MF)") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", color = "darkred"),
    strip.text = element_text(size = 15, face = "bold", color = "black"),
    axis.text.y = element_text(size = 11, face = "bold") # 通路名字加粗
  )
ggsave("GO_dotplot_custom.jpeg", p_dot, width = 12, height = 12, dpi = 600, device = "jpeg", quality = 95, bg = "white")


# 读取GO富集分析结果
filteredGO <- read.table("GO_results.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)

# 只保留BP、CC和MF类别
filteredGO_sub <- filteredGO %>% filter(ONTOLOGY %in% c("BP", "CC", "MF"))

# 对每个类别筛选Top最显著的GO term
topGO <- filteredGO_sub %>%
  group_by(ONTOLOGY) %>%
  slice_min(order_by = p.adjust, n = 20) %>%
  ungroup() %>%
  arrange(ONTOLOGY, p.adjust)

# 确保每个类别内按p.adjust排序
topGO <- topGO %>% 
  group_by(ONTOLOGY) %>%
  mutate(Description = factor(Description, levels = rev(unique(Description)))) %>%
  ungroup()

# GeneRatio转为数值格式（主要用于气泡图横轴）
topGO$GeneRatio_num <- sapply(topGO$GeneRatio, function(x) {
  sp <- unlist(strsplit(as.character(x), "/"))
  as.numeric(sp[1]) / as.numeric(sp[2])
})

# 1. 绘制CC类（Cellular Component）的条形图
ccGO <- topGO %>% filter(ONTOLOGY == "CC")

jpeg("GO_barplot_CC.jpeg", width=12, height=8, units="in", res=600, quality=95, bg="white")
ccBarPlot <- ggbarplot(
  ccGO,
  x = "Description", y = "Count", fill = "ONTOLOGY", color = "white",
  xlab = "", palette = "Set1", legend = "none", sort.val = "desc", sort.by.groups = TRUE,
  position = position_dodge(0.9)
) + 
  rotate_x_text(75) +
  ggtitle("GO Enrichment Analysis - Cellular Component (CC)") +  # 添加标题
  theme(
    panel.background = element_blank(),
    axis.text.x = element_text(size = 10, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "darkblue") # 设置标题样式
  )
print(ccBarPlot)
dev.off()

# 2. 绘制BP类（Biological Process）的条形图
bpGO <- topGO %>% filter(ONTOLOGY == "BP")

jpeg("GO_barplot_BP.jpeg", width=12, height=8, units="in", res=600, quality=95, bg="white")
bpBarPlot <- ggbarplot(
  bpGO,
  x = "Description", y = "Count", fill = "ONTOLOGY", color = "white",
  xlab = "", palette = "Set2", legend = "none", sort.val = "desc", sort.by.groups = TRUE,
  position = position_dodge(0.9)
) + 
  rotate_x_text(75) +
  ggtitle("GO Enrichment Analysis - Biological Process (BP)") +  # 添加标题
  theme(
    panel.background = element_blank(),
    axis.text.x = element_text(size = 10, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "darkblue") # 设置标题样式
  )
print(bpBarPlot)
dev.off()

# 3. 绘制MF类（Molecular Function）的条形图
mfGO <- topGO %>% filter(ONTOLOGY == "MF")

jpeg("GO_barplot_MF.jpeg", width=12, height=15, units="in", res=600, quality=95, bg="white")
mfBarPlot <- ggbarplot(
  mfGO,
  x = "Description", y = "Count", fill = "ONTOLOGY", color = "white",
  xlab = "", palette = "Set3", legend = "none", sort.val = "desc", sort.by.groups = TRUE,
  position = position_dodge(0.9)
) + 
  rotate_x_text(75) +
  ggtitle("GO Enrichment Analysis - Molecular Function (MF)") +  # 添加标题
  theme(
    panel.background = element_blank(),
    axis.text.x = element_text(size = 10, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "darkblue") # 设置标题样式
  )
print(mfBarPlot)
dev.off()

cat('CC, BP, 和 MF分类的条形图已分别保存为三个PDF文件，并添加了标题。\n')



#------------------------------#
# 10. 关闭进度条               #
#------------------------------#
close(pb)
cat('全部分析及绘图完成！\n')
