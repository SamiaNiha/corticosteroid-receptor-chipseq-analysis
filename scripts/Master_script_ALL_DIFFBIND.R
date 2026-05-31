# === 1. Load DiffBind ===
library(DiffBind)

# === 2. Load sample sheet ===
gr <- dba(sampleSheet = "/Volumes/BAM Files/GR_samples.csv")

# === 3. Count reads  ===
gr <- dba.count(gr, summits = 200, minOverlap = 4)

# === 4. Save counted object ===
saveRDS(gr, "/Volumes/BAM Files/GR_dba_counted.rds")

# === 5. Normalisation ===
gr <- dba.normalize(gr, normalize = DBA_NORM_LIB, library = DBA_LIBSIZE_FULL)

# === 6. Define all 3 contrasts from 'Condition' metadata (BLAM, BLPM, FS) ===
gr <- dba.contrast(gr, categories = DBA_CONDITION)

# === 7. Run differential analysis ===
gr <- dba.analyze(gr)

# === 8. Save full analysed object ===
saveRDS(gr, "/Volumes/BAM Files/GR_dba_analyzed.rds")

# === 9. Save number of consensus peaks to file ===
consensus_count <- length(dba.peakset(gr, bRetrieve = TRUE))
write(consensus_count, file = "/Volumes/BAM Files/GR_consensus_peak_count.txt")

# === 10. Match contrasts by name and export differential peak tables ===

# Get contrast table
contrast_info <- dba.show(gr, bContrasts = TRUE)

# Define helper function to match and save contrasts
save_contrast_outputs <- function(gr_obj, contrast_name, prefix) {
  contrast_idx <- which(contrast_info$Description == contrast_name)
  result <- dba.report(gr_obj, contrast = contrast_idx)
  
  # Save all peaks
  write.csv(as.data.frame(result),
            paste0("/Volumes/BAM Files/", prefix, "_", gsub(" ", "_", contrast_name), "_diffbind.csv"),
            row.names = FALSE)
  
  # Save only FDR ≤ 0.05
  sig_result <- result[result$FDR <= 0.05]
  write.csv(as.data.frame(sig_result),
            paste0("/Volumes/BAM Files/", prefix, "_", gsub(" ", "_", contrast_name), "_sigFDR05.csv"),
            row.names = FALSE)
}

# Save each contrast
save_contrast_outputs(gr, "BLPM - FS", "GR")
save_contrast_outputs(gr, "BLAM - BLPM", "GR")
save_contrast_outputs(gr, "BLAM - FS", "GR")

# === 1. Load DiffBind ===
library(DiffBind)

# === 2. Load sample sheet ===
mr <- dba(sampleSheet = "/Volumes/BAM Files/MR_samples.csv")

# === 3. Count reads with correct settings ===
mr <- dba.count(mr, summits = 200, minOverlap = 4)
saveRDS(mr, "/Volumes/BAM Files/MR_dba_counted.rds")

# === 4. Normalisation ===
mr <- dba.normalize(mr, normalize = DBA_NORM_LIB, library = DBA_LIBSIZE_FULL)

# === 5. Define all 3 contrasts from 'Condition' metadata ===
mr <- dba.contrast(mr, categories = DBA_CONDITION)

# === 6. Run analysis ===
mr <- dba.analyze(mr)
saveRDS(mr, "/Volumes/BAM Files/MR_dba_analyzed.rds")

# === 7. Save consensus peak count ===
consensus_count <- length(dba.peakset(mr, bRetrieve = TRUE))
write(consensus_count, file = "/Volumes/BAM Files/MR_consensus_peak_count.txt")

# === 8. Export reports ===
contrast_info <- dba.show(mr, bContrasts = TRUE)

save_contrast_outputs <- function(dba_obj, contrast_name, prefix) {
  contrast_idx <- which(contrast_info$Description == contrast_name)
  result <- dba.report(dba_obj, contrast = contrast_idx)
  
  write.csv(as.data.frame(result),
            paste0("/Volumes/BAM Files/", prefix, "_", gsub(" ", "_", contrast_name), "_diffbind.csv"),
            row.names = FALSE)
  
  sig_result <- result[result$FDR <= 0.05]
  write.csv(as.data.frame(sig_result),
            paste0("/Volumes/BAM Files/", prefix, "_", gsub(" ", "_", contrast_name), "_sigFDR05.csv"),
            row.names = FALSE)
}

save_contrast_outputs(mr, "BLPM - FS", "MR")
save_contrast_outputs(mr, "BLAM - BLPM", "MR")
save_contrast_outputs(mr, "BLAM - FS", "MR")

# === 1. Load DiffBind ===
library(DiffBind)

# === 2. Load sample sheet ===
s10 <- dba(sampleSheet = "/Volumes/BAM Files/S10_samples.csv")

# === 3. Count reads with correct settings ===
s10 <- dba.count(s10, summits = 200, minOverlap = 4)
saveRDS(s10, "/Volumes/BAM Files/S10_dba_counted.rds")

# === 4. Normalisation ===
s10 <- dba.normalize(s10, normalize = DBA_NORM_LIB, library = DBA_LIBSIZE_FULL)

# === 5. Define all 3 contrasts ===
s10 <- dba.contrast(s10, categories = DBA_CONDITION)

# === 6. Run analysis ===
s10 <- dba.analyze(s10)
saveRDS(s10, "/Volumes/BAM Files/S10_dba_analyzed.rds")

# === 7. Save consensus peak count ===
consensus_count <- length(dba.peakset(s10, bRetrieve = TRUE))
write(consensus_count, file = "/Volumes/BAM Files/S10_consensus_peak_count.txt")

# === 8. Export reports ===
contrast_info <- dba.show(s10, bContrasts = TRUE)

save_contrast_outputs <- function(dba_obj, contrast_name, prefix) {
  contrast_idx <- which(contrast_info$Description == contrast_name)
  result <- dba.report(dba_obj, contrast = contrast_idx)
  
  write.csv(as.data.frame(result),
            paste0("/Volumes/BAM Files/", prefix, "_", gsub(" ", "_", contrast_name), "_diffbind.csv"),
            row.names = FALSE)
  
  sig_result <- result[result$FDR <= 0.05]
  write.csv(as.data.frame(sig_result),
            paste0("/Volumes/BAM Files/", prefix, "_", gsub(" ", "_", contrast_name), "_sigFDR05.csv"),
            row.names = FALSE)
}

save_contrast_outputs(s10, "BLPM - FS", "S10")
save_contrast_outputs(s10, "BLAM - BLPM", "S10")
save_contrast_outputs(s10, "BLAM - FS", "S10")






library(DiffBind)                         # for differential binding analysis results
library(ChIPseeker)                       # for annotating ChIP-seq peaks to genes
library(GenomicRanges)                   # for working with genomic coordinates
library(TxDb.Rnorvegicus.UCSC.rn6.refGene)  # rat genome annotation (UCSC rn6)
library(org.Rn.eg.db)                     # gene annotation for Rattus norvegicus
library(openxlsx)                         # to write Excel files
library(VennDiagram)                      # for visualising overlaps in gene targets

# set up input/output paths
input_path <- "/Volumes/BAM Files/"
output_path <- file.path(input_path, "DiffBind_Output")
dir.create(output_path, showWarnings = FALSE)  # make output folder

# load analysed DiffBind objects for each target
dba_list <- list(
  S10 = readRDS(file.path(input_path, "S10_dba_analyzed.rds")),
  GR  = readRDS(file.path(input_path, "GR_dba_analyzed.rds")),
  MR  = readRDS(file.path(input_path, "MR_dba_analyzed.rds"))
)

# function to go through each contrast for one target
process_contrasts <- function(dba_obj, target_name) {
  contrast_info <- dba.show(dba_obj, bContrasts = TRUE)
  print(contrast_info)  # just to double check what contrasts exist
  
  for (i in seq_len(nrow(contrast_info))) {
    # get condition names for the contrast (e.g. BLAM vs FS)
    cond1 <- as.character(contrast_info$Group[i])
    cond2 <- as.character(contrast_info$Group2[i])
    
    # fallback names in case something is missing
    if (length(cond1) == 0 || is.na(cond1) || cond1 == "") cond1 <- "Unknown1"
    if (length(cond2) == 0 || is.na(cond2) || cond2 == "") cond2 <- "Unknown2"
    
    contrast_name <- paste0(cond1, "_vs_", cond2)
    cat("Working on:", target_name, "-", contrast_name, "\n")
    
    # pull out DB peaks for this contrast (FDR < 0.1 as per IPA input)
    db_peaks <- dba.report(dba_obj, contrast = i, th = 0.1)
    
    if (length(db_peaks) == 0) {
      message("No significant peaks for ", contrast_name, " in ", target_name)
      next
    }
    
    # make sure chromosome naming matches the annotation style
    seqlevelsStyle(db_peaks) <- "UCSC"
    
    # annotate peaks to genes (ENSEMBL, SYMBOL, etc)
    annotated <- annotatePeak(
      db_peaks,
      TxDb = TxDb.Rnorvegicus.UCSC.rn6.refGene,
      annoDb = "org.Rn.eg.db"
    )
    
    # convert to dataframe
    df <- as.data.frame(annotated)
    
    # only keep the columns I actually need
    keep <- c("seqnames", "start", "end", "geneId", "ENSEMBL", "SYMBOL",
              "log2FoldChange", "pvalue", "FDR", "annotation")
    keep <- keep[keep %in% colnames(df)]
    df <- df[, keep, drop = FALSE]
    colnames(df)[1:3] <- c("chr", "start", "end")  # rename for clarity
    
    # save to Excel file – named by target and contrast
    file_out <- file.path(output_path, paste0(target_name, "_", contrast_name, "_FDR0.1.xlsx"))
    write.xlsx(df, file = file_out)
    
    # check if FKBP5 shows up
    fkbp5 <- df[df$SYMBOL == "Fkbp5", ]
    if (nrow(fkbp5) > 0) {
      message("FKBP5 found in ", target_name, " (", contrast_name, ")")
      cols_to_show <- c("chr", "start", "end", "log2FoldChange", "FDR", "annotation")
      cols_to_show <- cols_to_show[cols_to_show %in% colnames(fkbp5)]
      print(fkbp5[, cols_to_show, drop = FALSE])
    } else {
      message("FKBP5 not found in ", target_name, " (", contrast_name, ")")
    }
  }
}

# run the function on S10, GR, and MR
for (target in names(dba_list)) {
  process_contrasts(dba_list[[target]], target)
}

message("\n Done – Excel files saved to: ", output_path)




















# === Load libraries ===
library(DiffBind)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(DESeq2)

# === Load your DiffBind object ===
dba_obj <- readRDS("/Volumes/BAM Files/GR_dba_analyzed.rds")  # Change to MR or S10 as needed
prefix <- "GR"  # Change to "MR" or "S10" accordingly

# === Set output directory ===
outdir <- "/Volumes/BAM Files/Figures/"
dir.create(outdir, showWarnings = FALSE)

# === PCA Plot ===
pdf(file = file.path(outdir, paste0(prefix, "_PCA.pdf")))
dba.plotPCA(dba_obj, contrast = 1, label = DBA_ID)
dev.off()

# === Heatmap Plot ===
pdf(file = file.path(outdir, paste0(prefix, "_heatmap.pdf")))
dba.plotHeatmap(dba_obj, contrast = 1)
dev.off()

# === MA & Volcano plots for all contrasts ===
contrast_info <- dba.show(dba_obj, bContrasts = TRUE)

for (i in seq_len(nrow(contrast_info))) {
  cname <- gsub(" ", "_", contrast_info$Description[i])
  
  # MA Plot
  pdf(file = file.path(outdir, paste0(prefix, "_", cname, "_MAplot.pdf")))
  dba.plotMA(dba_obj, contrast = i)
  dev.off()
  
  # Volcano Plot
  pdf(file = file.path(outdir, paste0(prefix, "_", cname, "_Volcano.pdf")))
  dba.plotVolcano(dba_obj, contrast = i, method = DBA_DESEQ2)
  dev.off()
}









# === Load libraries ===
library(DiffBind)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(DESeq2)

# === Load MR DiffBind object ===
dba_obj <- readRDS("/Volumes/BAM Files/MR_dba_analyzed.rds")
prefix <- "MR"

# === Set output directory ===
outdir <- "/Volumes/BAM Files/Figures/"
dir.create(outdir, showWarnings = FALSE)

# === PCA Plot ===
pdf(file = file.path(outdir, paste0(prefix, "_PCA.pdf")))
dba.plotPCA(dba_obj, contrast = 1, label = DBA_ID)
dev.off()

# === Heatmap ===
pdf(file = file.path(outdir, paste0(prefix, "_heatmap.pdf")))
dba.plotHeatmap(dba_obj, contrast = 1)
dev.off()

# === MA and Volcano Plots ===
contrast_info <- dba.show(dba_obj, bContrasts = TRUE)

for (i in seq_len(nrow(contrast_info))) {
  cname <- gsub(" ", "_", contrast_info$Description[i])
  
  # MA Plot
  pdf(file = file.path(outdir, paste0(prefix, "_", cname, "_MAplot.pdf")))
  dba.plotMA(dba_obj, contrast = i)
  dev.off()
  
  # Volcano Plot
  pdf(file = file.path(outdir, paste0(prefix, "_", cname, "_Volcano.pdf")))
  dba.plotVolcano(dba_obj, contrast = i, method = DBA_DESEQ2)
  dev.off()
}







# === Load libraries ===
library(DiffBind)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(DESeq2)

# === Load S10 DiffBind object ===
dba_obj <- readRDS("/Volumes/BAM Files/S10_dba_analyzed.rds")
prefix <- "S10"

# === Set output directory ===
outdir <- "/Volumes/BAM Files/Figures/"
dir.create(outdir, showWarnings = FALSE)

# === PCA Plot ===
pdf(file = file.path(outdir, paste0(prefix, "_PCA.pdf")))
dba.plotPCA(dba_obj, contrast = 1, label = DBA_ID)
dev.off()

# === Heatmap ===
pdf(file = file.path(outdir, paste0(prefix, "_heatmap.pdf")))
dba.plotHeatmap(dba_obj, contrast = 1)
dev.off()

# === MA and Volcano Plots ===
contrast_info <- dba.show(dba_obj, bContrasts = TRUE)

for (i in seq_len(nrow(contrast_info))) {
  cname <- gsub(" ", "_", contrast_info$Description[i])
  
  # MA Plot
  pdf(file = file.path(outdir, paste0(prefix, "_", cname, "_MAplot.pdf")))
  dba.plotMA(dba_obj, contrast = i)
  dev.off()
  
  # Volcano Plot
  pdf(file = file.path(outdir, paste0(prefix, "_", cname, "_Volcano.pdf")))
  dba.plotVolcano(dba_obj, contrast = i, method = DBA_DESEQ2)
  dev.off()
}








