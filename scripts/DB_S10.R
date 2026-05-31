# ==== 1. Load packages ====
library(DiffBind)

# ==== 2. Load sample sheet and create DBA object ====
s10 <- dba(sampleSheet = file.path("/Volumes", "BAM Files", "S10_samples.csv"))

# QC heatmap before counting (based on peak overlap)
plot(s10)

# ==== 3. Count reads for affinity matrix ====
# Standardises peaks to 401bp windows centred on summits
s10 <- dba.count(s10, summits = 200)

# Save counted object
saveRDS(s10, "s10_counted.rds")

# ==== 4. Normalisation ====
# Library-size normalisation based on total aligned reads
s10 <- dba.normalize(s10, normalize = DBA_NORM_LIB, library = DBA_LIBSIZE_FULL)

# View normalisation details
norm_s10 <- dba.normalize(s10, bRetrieve = TRUE)
print(norm_s10)

# ==== 5. Define contrast: FS vs BLPM ====
# Assumes 'Condition' column contains "FS" and "BLPM"
s10 <- dba.contrast(s10, categories = DBA_CONDITION, reorderMeta = list(Condition = "FS"))

# ==== 6. Run differential binding analysis ====
s10 <- dba.analyze(s10)

# ==== 7. Show contrast summary ====
dba.show(s10, bContrasts = TRUE)










# === 1. Load required libraries ===
library(DiffBind)

# === 2. Save the analyzed object ===
# Assumes you have already run dba.contrast() and dba.analyze()
saveRDS(s10, file = "/Volumes/BAM Files/S10_DiffBind_Analyzed.rds")

# === 3. Extract and save FS vs BLPM results (contrast 3) ===
db_fs_vs_blpm <- dba.report(s10, contrast = 3)
write.csv(as.data.frame(db_fs_vs_blpm),
          "/Volumes/BAM Files/S10_FS_vs_BLPM_All_Peaks.csv",
          row.names = FALSE)

# === 4. Save significant FS vs BLPM peaks (FDR ≤ 0.05) ===
sig_fs_vs_blpm <- db_fs_vs_blpm[db_fs_vs_blpm$FDR <= 0.05]
write.csv(as.data.frame(sig_fs_vs_blpm),
          "/Volumes/BAM Files/S10_FS_vs_BLPM_Significant_Peaks_FDR0.05.csv",
          row.names = FALSE)

# === 5. Save FS vs BLAM results (contrast 1) ===
db_fs_vs_blam <- dba.report(s10, contrast = 1)
write.csv(as.data.frame(db_fs_vs_blam),
          "/Volumes/BAM Files/S10_FS_vs_BLAM_All_Peaks.csv",
          row.names = FALSE)

# === 6. Save BLPM vs BLAM results (contrast 2) ===
db_blpm_vs_blam <- dba.report(s10, contrast = 2)
write.csv(as.data.frame(db_blpm_vs_blam),
          "/Volumes/BAM Files/S10_BLPM_vs_BLAM_All_Peaks.csv",
          row.names = FALSE)

# === 7. Save normalisation info ===
norm_info <- dba.normalize(s10, bRetrieve = TRUE)
write.csv(as.data.frame(norm_info),
          "/Volumes/BAM Files/S10_Normalisation_Info.csv",
          row.names = FALSE)

# === 8. Save sample metadata ===
sample_info <- dba.show(s10)
write.csv(as.data.frame(sample_info),
          "/Volumes/BAM Files/S10_Sample_Metadata.csv",
          row.names = FALSE)

# === 9. Volcano plot (FS vs BLPM) ===
png("/Volumes/BAM Files/S10_FS_vs_BLPM_VolcanoPlot.png", width = 1000, height = 800)
dba.plotVolcano(s10, contrast = 3)
dev.off()

# === 10. Heatmap (FS vs BLPM) ===
png("/Volumes/BAM Files/S10_FS_vs_BLPM_Heatmap.png", width = 1200, height = 1000)
dba.plotHeatmap(s10, contrast = 3, correlations = FALSE, scale = "row")
dev.off()

# === 11. PCA plot (FS vs BLPM) ===
png("/Volumes/BAM Files/S10_PCA.png", width = 1000, height = 800)
dba.plotPCA(s10, contrast = 3, label = DBA_ID)
dev.off()

# === 12. Save session info ===
sink("/Volumes/BAM Files/S10_SessionInfo.txt")
sessionInfo()
sink()



