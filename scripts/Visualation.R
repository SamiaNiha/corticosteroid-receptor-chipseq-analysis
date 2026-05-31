# ==== Load required libraries ====
library(DiffBind)
library(ChIPseeker)
library(TxDb.Rnorvegicus.UCSC.rn6.refGene)
library(org.Rn.eg.db)
library(AnnotationDbi)
library(GenomeInfoDb)
library(dplyr)
library(ggplot2)

# ==== Load GR object ====
gr <- readRDS("/Volumes/BAM Files/gr_counted.rds")
gr <- dba.contrast(gr, categories = DBA_CONDITION, reorderMeta = list(Condition = "FS"))
gr <- dba.analyze(gr)

# ==== Report Peaks and Annotate ====
gr_peaks <- dba.report(gr, contrast = 3)
seqlevelsStyle(gr_peaks) <- "UCSC"
txdb <- TxDb.Rnorvegicus.UCSC.rn6.refGene
gr_anno <- annotatePeak(gr_peaks, TxDb = txdb, tssRegion = c(-3000, 3000), annoDb = "org.Rn.eg.db")
anno_df <- as.data.frame(gr_anno)
gr_report <- as.data.frame(gr_peaks)
gr_report$SYMBOL <- anno_df$SYMBOL
symbol2ensembl <- mapIds(org.Rn.eg.db, keys = gr_report$SYMBOL,
                         column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
gr_report$ENSEMBL <- symbol2ensembl[gr_report$SYMBOL]

# ==== Prepare data for IPA ====
ipa_df <- gr_report %>%
  select(ENSEMBL, SYMBOL, Fold, p.value, FDR) %>%
  rename(`Ensembl ID` = ENSEMBL,
         `Gene Symbol` = SYMBOL,
         `log2 Fold Change` = Fold,
         `p-value` = p.value,
         `FDR` = FDR) %>%
  filter(!is.na(`Gene Symbol`)) %>%
  mutate(absFC = abs(`log2 Fold Change`)) %>%
  group_by(`Gene Symbol`) %>%
  slice_max(absFC, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(-absFC)

write.csv(ipa_df, "/Volumes/BAM Files/GR_FS_vs_BLPM_IPA_Ready.csv", row.names = FALSE)
write.table(unique(ipa_df$`Gene Symbol`),
            "/Volumes/BAM Files/GR_FS_vs_BLPM_Gene_List.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# ==== Save visualisation plots as PNG ====

# PCA plot
png("/Volumes/BAM Files/GR_PCA_plot.png", width = 800, height = 600)
dba.plotPCA(gr, contrast = 3, label = DBA_ID)
dev.off()

# Sample correlation heatmap
png("/Volumes/BAM Files/GR_Heatmap_All.png", width = 800, height = 600)
dba.plotHeatmap(gr)
dev.off()

# Differential peaks heatmap
png("/Volumes/BAM Files/GR_Heatmap_DB_Only.png", width = 800, height = 600)
dba.plotHeatmap(gr, contrast = 3, correlations = FALSE, scale = "row")
dev.off()

# MA plot
png("/Volumes/BAM Files/GR_MA_plot.png", width = 800, height = 600)
dba.plotMA(gr, contrast = 3)
dev.off()

# Volcano plot
png("/Volumes/BAM Files/GR_Volcano_plot.png", width = 800, height = 600)
dba.plotVolcano(gr, contrast = 3)
dev.off()

# Boxplot of read distributions
png("/Volumes/BAM Files/GR_Boxplot_Reads.png", width = 800, height = 600)
dba.plotBox(gr)
dev.off()

# Annotation pie chart
png("/Volumes/BAM Files/GR_Annotation_PieChart.png", width = 800, height = 600)
plotAnnoPie(gr_anno)
dev.off()

# TSS distance barplot
png("/Volumes/BAM Files/GR_TSS_Barplot.png", width = 800, height = 600)
plotDistToTSS(gr_anno, title = "GR FS vs BLPM: Distance to TSS")
dev.off()










# ==== Load required libraries ====
library(DiffBind)
library(ChIPseeker)
library(TxDb.Rnorvegicus.UCSC.rn6.refGene)
library(org.Rn.eg.db)
library(AnnotationDbi)
library(GenomeInfoDb)
library(dplyr)
library(ggplot2)

# ==== Load MR object ====
mr <- readRDS("/Volumes/BAM Files/mr_counted.rds")
mr <- dba.contrast(mr, categories = DBA_CONDITION, reorderMeta = list(Condition = "FS"))
mr <- dba.analyze(mr)

# ==== Report Peaks and Annotate ====
mr_peaks <- dba.report(mr, contrast = 3)
seqlevelsStyle(mr_peaks) <- "UCSC"
txdb <- TxDb.Rnorvegicus.UCSC.rn6.refGene
mr_anno <- annotatePeak(mr_peaks, TxDb = txdb, tssRegion = c(-3000, 3000), annoDb = "org.Rn.eg.db")
anno_df <- as.data.frame(mr_anno)
mr_report <- as.data.frame(mr_peaks)
mr_report$SYMBOL <- anno_df$SYMBOL
symbol2ensembl <- mapIds(org.Rn.eg.db, keys = mr_report$SYMBOL,
                         column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
mr_report$ENSEMBL <- symbol2ensembl[mr_report$SYMBOL]

# ==== Prepare data for IPA ====
ipa_df <- mr_report %>%
  select(ENSEMBL, SYMBOL, Fold, p.value, FDR) %>%
  rename(`Ensembl ID` = ENSEMBL,
         `Gene Symbol` = SYMBOL,
         `log2 Fold Change` = Fold,
         `p-value` = p.value,
         `FDR` = FDR) %>%
  filter(!is.na(`Gene Symbol`)) %>%
  mutate(absFC = abs(`log2 Fold Change`)) %>%
  group_by(`Gene Symbol`) %>%
  slice_max(absFC, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(-absFC)

write.csv(ipa_df, "/Volumes/BAM Files/MR_FS_vs_BLPM_IPA_Ready.csv", row.names = FALSE)
write.table(unique(ipa_df$`Gene Symbol`),
            "/Volumes/BAM Files/MR_FS_vs_BLPM_Gene_List.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# ==== Save visualisation plots as PNG ====
png("/Volumes/BAM Files/MR_PCA_plot.png", width = 800, height = 600)
dba.plotPCA(mr, contrast = 3, label = DBA_ID)
dev.off()

png("/Volumes/BAM Files/MR_Heatmap_All.png", width = 800, height = 600)
dba.plotHeatmap(mr)
dev.off()

png("/Volumes/BAM Files/MR_Heatmap_DB_Only.png", width = 800, height = 600)
dba.plotHeatmap(mr, contrast = 3, correlations = FALSE, scale = "row")
dev.off()

png("/Volumes/BAM Files/MR_MA_plot.png", width = 800, height = 600)
dba.plotMA(mr, contrast = 3)
dev.off()

png("/Volumes/BAM Files/MR_Volcano_plot.png", width = 800, height = 600)
dba.plotVolcano(mr, contrast = 3)
dev.off()

png("/Volumes/BAM Files/MR_Boxplot_Reads.png", width = 800, height = 600)
dba.plotBox(mr)
dev.off()

png("/Volumes/BAM Files/MR_Annotation_PieChart.png", width = 800, height = 600)
plotAnnoPie(mr_anno)
dev.off()

png("/Volumes/BAM Files/MR_TSS_Barplot.png", width = 800, height = 600)
plotDistToTSS(mr_anno, title = "MR FS vs BLPM: Distance to TSS")
dev.off()









# ==== Load required libraries ====
library(DiffBind)
library(ChIPseeker)
library(TxDb.Rnorvegicus.UCSC.rn6.refGene)
library(org.Rn.eg.db)
library(AnnotationDbi)
library(GenomeInfoDb)
library(dplyr)
library(ggplot2)

# ==== Load S10 object ====
s10 <- readRDS("/Volumes/BAM Files/s10_counted.rds")
s10 <- dba.contrast(s10, categories = DBA_CONDITION, reorderMeta = list(Condition = "FS"))
s10 <- dba.analyze(s10)

# ==== Report Peaks and Annotate ====
s10_peaks <- dba.report(s10, contrast = 3)
seqlevelsStyle(s10_peaks) <- "UCSC"
txdb <- TxDb.Rnorvegicus.UCSC.rn6.refGene
s10_anno <- annotatePeak(s10_peaks, TxDb = txdb, tssRegion = c(-3000, 3000), annoDb = "org.Rn.eg.db")
anno_df <- as.data.frame(s10_anno)
s10_report <- as.data.frame(s10_peaks)
s10_report$SYMBOL <- anno_df$SYMBOL
symbol2ensembl <- mapIds(org.Rn.eg.db, keys = s10_report$SYMBOL,
                         column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
s10_report$ENSEMBL <- symbol2ensembl[s10_report$SYMBOL]

# ==== Prepare data for IPA ====
ipa_df <- s10_report %>%
  select(ENSEMBL, SYMBOL, Fold, p.value, FDR) %>%
  rename(`Ensembl ID` = ENSEMBL,
         `Gene Symbol` = SYMBOL,
         `log2 Fold Change` = Fold,
         `p-value` = p.value,
         `FDR` = FDR) %>%
  filter(!is.na(`Gene Symbol`)) %>%
  mutate(absFC = abs(`log2 Fold Change`)) %>%
  group_by(`Gene Symbol`) %>%
  slice_max(absFC, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(-absFC)

write.csv(ipa_df, "/Volumes/BAM Files/S10_FS_vs_BLPM_IPA_Ready.csv", row.names = FALSE)
write.table(unique(ipa_df$`Gene Symbol`),
            "/Volumes/BAM Files/S10_FS_vs_BLPM_Gene_List.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# ==== Save visualisation plots as PNG ====
png("/Volumes/BAM Files/S10_PCA_plot.png", width = 800, height = 600)
dba.plotPCA(s10, contrast = 3, label = DBA_ID)
dev.off()

png("/Volumes/BAM Files/S10_Heatmap_All.png", width = 800, height = 600)
dba.plotHeatmap(s10)
dev.off()

png("/Volumes/BAM Files/S10_Heatmap_DB_Only.png", width = 800, height = 600)
dba.plotHeatmap(s10, contrast = 3, correlations = FALSE, scale = "row")
dev.off()

png("/Volumes/BAM Files/S10_MA_plot.png", width = 800, height = 600)
dba.plotMA(s10, contrast = 3)
dev.off()

png("/Volumes/BAM Files/S10_Volcano_plot.png", width = 800, height = 600)
dba.plotVolcano(s10, contrast = 3)
dev.off()

png("/Volumes/BAM Files/S10_Boxplot_Reads.png", width = 800, height = 600)
dba.plotBox(s10)
dev.off()

png("/Volumes/BAM Files/S10_Annotation_PieChart.png", width = 800, height = 600)
plotAnnoPie(s10_anno)
dev.off()

png("/Volumes/BAM Files/S10_TSS_Barplot.png", width = 800, height = 600)
plotDistToTSS(s10_anno, title = "S10 FS vs BLPM: Distance to TSS")
dev.off()

