

# ---------------------------------------------------
# Load necessary packages
# ---------------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(RColorBrewer)
library(circlize)

# ---------------------------------------------------
# 3. Parameter settings
# ---------------------------------------------------
pvalueFilter <- 0.05
adjPvalFilter <- 0.05
colorSel <- "p.adjust"
if (adjPvalFilter > 1) {
  colorSel <- "pvalue"
}

# ---------------------------------------------------
# 4. Set working directory from the first positional argument (or current dir)
# ---------------------------------------------------
cli_args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1L]]) else "scripts/05_KEGG_enrichment.R"
package_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
input_file <- normalizePath(
  if (length(cli_args) >= 2L) cli_args[[2L]] else file.path(package_root, "data", "processed", "shared_targets_75.csv"),
  winslash = "/",
  mustWork = TRUE
)
work_dir <- if (length(cli_args) >= 1L) cli_args[[1L]] else file.path(package_root, "reproduced_outputs", "KEGG_submitted_75")
if (!dir.exists(work_dir)) dir.create(work_dir, recursive = TRUE)
setwd(work_dir)

# ---------------------------------------------------
# 5. Read the input file (assume the first column contains gene symbols)
# ---------------------------------------------------
if (!file.exists(input_file)) stop("Input gene file not found: ", input_file)
rt <- read.csv(input_file, header = TRUE, check.names = FALSE)
if ("Gene" %in% names(rt)) names(rt)[names(rt) == "Gene"] <- "genes"
if (!("genes" %in% names(rt))) names(rt)[1] <- "genes"

# ---------------------------------------------------
# 6. Extract gene names and convert to Entrez IDs
# ---------------------------------------------------
genes <- unique(as.vector(rt$genes))
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound = NA)
entrezIDs <- as.character(entrezIDs)
rt$entrezIDs <- entrezIDs
rt <- rt[rt$entrezIDs != "NA", ]
gene <- rt$entrezIDs  # Filtered Entrez IDs

# ---------------------------------------------------
# 7. Perform KEGG enrichment analysis
# ---------------------------------------------------
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff = 1, qvalueCutoff = 1)

# ---------------------------------------------------
# 8. Convert enrichment results to a readable format (replace Entrez IDs with gene symbols)
# ---------------------------------------------------
kkx <- setReadable(kk, 'org.Hs.eg.db', 'ENTREZID')
KEGG <- as.data.frame(kkx)

# Optionally, if manual conversion is needed, use:
# KEGG$geneID <- sapply(KEGG$geneID, function(x) {
#   id_vector <- unlist(strsplit(x, "/"))
#   idx <- match(id_vector, rt$entrezIDs)
#   paste(rt$genes[idx], collapse = "/")
# })

# Filter significant enrichment results using the manuscript's BH-FDR criterion.
KEGG <- KEGG[(KEGG$pvalue < pvalueFilter & KEGG$p.adjust < adjPvalFilter), ]

# ---------------------------------------------------
# 9. Output the significant enrichment results (with gene symbols)
# ---------------------------------------------------
write.table(KEGG, file = "KEGG.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Define the number of pathways to display
showNum <- 30
if (nrow(KEGG) < showNum) {
  showNum <- nrow(KEGG)
}

# ---------------------------------------------------
# 10. Create a bar plot (enhanced version with English labels)
# ---------------------------------------------------
topKEGG <- KEGG[order(KEGG$p.adjust), ][1:showNum, ]
topKEGG$Description <- factor(topKEGG$Description, levels = rev(topKEGG$Description))

p_bar <- ggplot(topKEGG, aes(x = Description, y = -log10(p.adjust), fill = -log10(p.adjust))) +
  geom_bar(stat = "identity", width = 0.8) +
  coord_flip() +
  scale_fill_gradientn(colors = brewer.pal(7, "YlOrRd")) +
  theme_minimal(base_size = 14) +
  labs(x = "KEGG Pathway", y = "-log10(Adjusted p-value)", title = "KEGG Pathway Enrichment Analysis") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", color = "darkblue"))
ggsave("barplot_custom.jpeg", p_bar, width = 8, height = 7, dpi = 600, device = "jpeg", quality = 95, bg = "white")

# ---------------------------------------------------
# 11. Create a dot plot (enhanced version with English labels)
# ---------------------------------------------------
if(all(grepl("/", topKEGG$GeneRatio))){
  topKEGG$GeneRatio_num <- sapply(topKEGG$GeneRatio, function(x) {
    parts <- unlist(strsplit(x, "/"))
    as.numeric(parts[1]) / as.numeric(parts[2])
  })
} else {
  topKEGG$GeneRatio_num <- as.numeric(topKEGG$GeneRatio)
}

p_dot <- ggplot(topKEGG, aes(x = GeneRatio_num, y = Description, size = Count, color = p.adjust)) +
  geom_point(alpha = 0.8) +
  scale_color_gradientn(colors = brewer.pal(7, "Spectral")) +
  theme_classic(base_size = 14) +
  labs(x = "Gene Ratio", y = "KEGG Pathway", title = "KEGG Dotplot") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", color = "darkred"))
ggsave("dotplot_custom.jpeg", p_dot, width = 8, height = 7, dpi = 600, device = "jpeg", quality = 95, bg = "white")

# ---------------------------------------------------
# 12. Create a Gene-Pathway Network plot (cnetplot, enhanced version with English labels)
# ---------------------------------------------------
jpeg(filename = "cnetplot_custom.jpeg", width = 9, height = 5.25, units = "in", res = 600, quality = 95, bg = "white")
kkx <- setReadable(kk, 'org.Hs.eg.db', 'ENTREZID')
p_cnet <- cnetplot(kkx, circular = TRUE, showCategory = 5, colorEdge = TRUE)
print(p_cnet + 
        ggtitle("Gene-Pathway Network") +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", color = "purple")))
dev.off()

# ---------------------------------------------------
# 13. Create an enhanced Chord Diagram with English labels
# ---------------------------------------------------
# Select the top 5 pathways (adjust as needed)
topKEGG_chord <- KEGG[1:min(10, nrow(KEGG)), ]
gene_path_df <- data.frame(Pathway = character(), Gene = character(), stringsAsFactors = FALSE)
for (i in 1:nrow(topKEGG_chord)) {
  pathway <- topKEGG_chord$Description[i]
  genes <- unlist(strsplit(topKEGG_chord$geneID[i], "/"))
  temp_df <- data.frame(Pathway = rep(pathway, length(genes)), Gene = genes, stringsAsFactors = FALSE)
  gene_path_df <- rbind(gene_path_df, temp_df)
}
print(head(gene_path_df))

chord_matrix <- as.matrix(table(gene_path_df$Pathway, gene_path_df$Gene))
rownames(chord_matrix) <- as.character(rownames(chord_matrix))
colnames(chord_matrix) <- as.character(colnames(chord_matrix))
print(dim(chord_matrix))
print(chord_matrix)

allSectors <- union(rownames(chord_matrix), colnames(chord_matrix))
print(allSectors)

nPathways <- length(rownames(chord_matrix))
nGenes <- length(colnames(chord_matrix))
pathway_colors <- brewer.pal(n = min(nPathways, 8), name = "Set2")
if(nPathways > length(pathway_colors)){
  pathway_colors <- rep(pathway_colors, length.out = nPathways)
}
gene_colors <- brewer.pal(n = min(nGenes, 8), name = "Pastel1")
if(nGenes > length(gene_colors)){
  gene_colors <- rep(gene_colors, length.out = nGenes)
}

grid.col <- c(setNames(pathway_colors, rownames(chord_matrix)),
              setNames(gene_colors, colnames(chord_matrix)))
grid.col <- grid.col[allSectors]
print(grid.col)

jpeg(filename = "chordDiagram_custom.jpeg", width = 8, height = 8, units = "in", res = 600, quality = 95, bg = "white")
chordDiagram(chord_matrix, grid.col = grid.col, transparency = 0.25, 
             annotationTrack = "grid", preAllocateTracks = list(track.height = 0.05))
title("Gene-Pathway Chord Diagram")
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  sector.name <- get.cell.meta.data("sector.index")
  circos.text(CELL_META$xcenter, CELL_META$ylim[1] - mm_y(5),
              sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.8)
}, bg.border = NA)
dev.off()

# ---------------------------------------------------
# End of analysis
# ---------------------------------------------------
