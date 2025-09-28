# set working directory
setwd('/Users/HP/Documents/HACKBIO/PlantProject/')

# Load packages to work with
library(DESeq2)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(KEGGREST)
library(DOSE)
library(org.At.tair.db)
library(enrichplot)

# Load data

a_t_count <- read.delim('counts.txt', header = T)
a_t_meta <- read.delim('metadata.tsv', header = T, stringsAsFactors = TRUE)

#preview
head(a_t_count)
head(a_t_meta)

raw_counts <- a_t_count[ , c(1, 7:12) ]

# Add gene IDs as rownames
rownames(raw_counts) <- a_t_count$Geneid

# Remove the Gene id column
raw_counts <- raw_counts[ , -1 ]

# preview raw_counts
head(raw_counts)

a_t_meta$condition <- factor(a_t_meta$condition, levels = c("control", "treatment")) 

dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = a_t_meta,
                              design = ~ condition)

#Run DESeq

dds <- DESeq(dds)
final_res <- results(dds)    

# Inspect results
head(final_res)
summary(final_res)



# Subset upregulated and downregulated genes
upregulated   <- subset(final_res, padj < 0.05 & log2FoldChange > 1)
downregulated <- subset(final_res, padj < 0.05 & log2FoldChange < -1)

# Basic volcano plot with all genes in grey
plot(final_res$log2FoldChange,  
     -log10(final_res$padj), 
     pch = 19,
     cex = 0.4,
     col = "grey",
     xlab = "Log2 Fold Change (Treatment vs Control)", 
     ylab = "-log10 Adjusted P-value",
     main = "Volcano Plot of Differential Expression")

# Add threshold lines
abline(v = c(-1, 1), col = "darkblue", lwd = 0.7, lty = 2)
abline(h = -log10(0.05), col = "darkred", lwd = 0.7, lty = 2)

# Highlight significant genes
points(upregulated$log2FoldChange,   -log10(upregulated$padj),   col = "salmon",  pch = 19, cex = 0.6)
points(downregulated$log2FoldChange, -log10(downregulated$padj), col = "skyblue", pch = 19, cex = 0.6)

legend("topright", 
       legend = c("Upregulated", "Downregulated"),
       col = c("salmon", "skyblue"), 
       pch = 19, 
       cex = 0.8)

# Heat Map
# VST-normalized counts
vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)

# Top 50 DEGs by padj
sig_degs <- subset(final_res, padj < 0.05)
top50_genes <- head(sig_degs[order(sig_degs$padj), ], 50)
vsd_top50 <- vsd_mat[rownames(top50_genes), ]

# Remove "gene:" prefix from gene names (handles "gene:AT1..." and "gene: AT1..." variants)
clean_names <- gsub("^gene:?\\s*", "", rownames(vsd_top50))
rownames(vsd_top50) <- clean_names

# Sample annotation (no visible title)
sample_group <- data.frame(Condition = colData(dds)$condition)
rownames(sample_group) <- colnames(vsd_top50)
colnames(sample_group) <- " "   # set to a single space so no visible title appears

# Annotation colors must match the (blank) column name
ann_colors <- list(
  " " = c("control" = "skyblue", "treatment" = "tomato")
)

# Plot
pheatmap(vsd_top50,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         annotation_col = sample_group,
         annotation_colors = ann_colors,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Top 50 DEGs Heatmap (UVC Treatment vs Control)")

sig_degs <- final_res %>%
  as.data.frame() %>%
  filter(!is.na(padj)) %>%
  filter(abs(log2FoldChange) > 2.5, padj < 0.05)

top100_degs <- sig_degs %>%
  arrange(padj) %>%
  head(100)

head(top100_degs)
# Enrichement analysis
# Extract gene IDs from top100 DEGs
gene_list <- gsub("^gene:", "", rownames(top100_degs))

gene_entrez <- mapIds(org.At.tair.db,
                      keys = gene_list,
                      column = "ENTREZID",
                      keytype = "TAIR",
                      multiVals = "first")
entrez_ids <- na.omit(gene_entrez)

ego_BP <- enrichGO(
  gene          = entrez_ids,
  OrgDb         = org.At.tair.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

ego_MF <- enrichGO(
  gene          = entrez_ids,
  OrgDb         = org.At.tair.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

ego_CC <- enrichGO(
  gene          = entrez_ids,
  OrgDb         = org.At.tair.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

BP_on <- as.data.frame(ego_BP); BP_on$Ontology <- "BP"
MF_on <- as.data.frame(ego_MF); MF_on$Ontology <- "MF"
CC_on <- as.data.frame(ego_CC); CC_on$Ontology <- "CC"
GO_merged <- rbind(BP_on, MF_on, CC_on)


# Dotplot for GO Biological Process
dotplot(ego_BP, showCategory = 15) +
  ggtitle("GO Enrichment: Biological Process (Top 15)") +
  theme(
    plot.title   = element_text(size = 18, face = "bold"),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 8),
    legend.text  = element_text(size = 12),
    legend.title = element_text(size = 14)
  )

# Dotplot for GO Molecular Function
dotplot(ego_MF, showCategory = 10) +
  ggtitle("GO Enrichment: Molecular Function (Top 10)") +
  theme(
    plot.title   = element_text(size = 18, face = "bold"),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 7),
    legend.text  = element_text(size = 12),
    legend.title = element_text(size = 14)
  )

# Dotplot for GO Cellular Component
dotplot(ego_CC, showCategory = 20) +
  ggtitle("GO Enrichment: Cellular Component (Top 20)") +
  theme(
    plot.title   = element_text(size = 18, face = "bold"),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 10),
    legend.text  = element_text(size = 12),
    legend.title = element_text(size = 14)
  )

# KEGG enrichment using TAIR IDs
tair_ids <- gsub("^gene:?\\s*", "", rownames(final_res))  # remove "gene:" prefix
de_tair  <- tair_ids[final_res$padj < 0.05 & abs(final_res$log2FoldChange) > 1]


ekegg <- enrichKEGG(
  gene         = de_tair,
  organism     = "ath",
  keyType      = "kegg", 
  pAdjustMethod= "BH",
  qvalueCutoff = 0.05
)

# Show top 5 enriched pathways
top5_kegg <- as.data.frame(ekegg) %>%
  arrange(p.adjust) %>%
  head(5)
top5_kegg

# Dotplot for visualization
dotplot(ekegg, showCategory = 30) + 
  ggtitle("KEGG Enrichment: Top 100 DEGs") +
  theme(
    plot.title   = element_text(size = 18, face = "bold"),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 10),
    legend.text  = element_text(size = 12),
    legend.title = element_text(size = 14)
  )

# Outputs
write.csv(final_res,         "Final DESeq2 Results.csv")
write.csv(upregulated,          "Upregulated_Genes.csv")
write.csv(downregulated,        "Downregulated_Genes.csv")
write.csv(raw_counts,  "Raw_counts.csv")
write.csv(top100_degs, file = "Top100_DEGs.csv", row.names = TRUE)
write.csv(as.data.frame(ego_BP), file="GO_BP_results.csv", row.names=FALSE)
write.csv(as.data.frame(ego_MF), file="GO_MF_results.csv", row.names=FALSE)
write.csv(as.data.frame(ego_CC), file="GO_CC_results.csv", row.names=FALSE)
write.csv(GO_merged, file="GO_Enrichment_Merged.csv", row.names=FALSE)
write.csv(as.data.frame(ekegg), file="KEGG_Top100_DEGs.csv", row.names=FALSE)
write.csv(top5_kegg, "Top 5 KEGG pathways.csv", row.names = FALSE)