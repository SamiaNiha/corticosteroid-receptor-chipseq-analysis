# ==== 1. Load packages ====
library(DiffBind)

# ==== 2. Load sample sheet and create DBA object ====
mr <- dba(sampleSheet = file.path("/Volumes", "BAM Files", "MR_samples.csv"))

#QC heatmap before counting (based on peak overlap)
plot(mr)

# ==== 3. Count reads for affinity matrix ====
# Standardises peaks to 401bp windows centred on summits
mr <- dba.count(mr, summits = 200)

# Save counted object
saveRDS(mr, file = "/Volumes/BAM Files/mr_counted.rds")

# ==== 4. Normalisation ====
# Library-size based normalisation
mr <- dba.normalize(mr, normalize = DBA_NORM_LIB, library = DBA_LIBSIZE_FULL)

# View normalisation details
norm_mr <- dba.normalize(mr, bRetrieve = TRUE)
print(norm_mr)

# ==== 5. Define contrast: FS vs BLPM ====
# Assumes 'Condition' column contains "FS" and "BLPM"
mr <- dba.contrast(mr, categories = DBA_CONDITION, reorderMeta = list(Condition = "FS"))

# ==== 6. Run differential binding analysis ====
mr <- dba.analyze(mr)

# ==== 7. Show contrast summary ====
dba.show(mr, bContrasts = TRUE)







