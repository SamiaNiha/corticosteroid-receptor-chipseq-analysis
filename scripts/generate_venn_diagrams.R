
# Load necessary packages
library(DiffBind)
library(GenomicRanges)
library(VennDiagram)
library(gridExtra)

# Function to load peaksets and extract peaks
load_peaks <- function(rds_path) {
  obj <- readRDS(rds_path)
  dba.peakset(obj, bRetrieve = TRUE)
}

# Load peaksets for each factor under each condition
conditions <- c("BLAM", "FS", "BLPM")
factors <- c("MR", "GR", "S10")
peak_data <- list()

for (cond in conditions) {
  for (fac in factors) {
    label <- paste0(fac, "_", cond, "_3of4.rds")
    peaks <- load_peaks(label)
    key <- paste(fac, cond, sep = "_")
    peak_data[[key]] <- reduce(peaks)  # ensure non-overlapping GRanges
  }
}

# Helper: get overlaps between three GRanges
triple_overlap <- function(gr1, gr2, gr3) {
  ol12 <- findOverlaps(gr1, gr2)
  ol23 <- findOverlaps(gr2, gr3)
  ol13 <- findOverlaps(gr1, gr3)
  ol_all <- findOverlaps(gr1[queryHits(ol12)], gr3)
  length(unique(queryHits(ol_all)))
}

# Create Venn diagrams for each condition
pdf("Venn_Diagrams_by_Condition.pdf", width = 10, height = 10)
for (cond in conditions) {
  s10 <- peak_data[[paste0("S10_", cond)]]
  mr  <- peak_data[[paste0("MR_", cond)]]
  gr  <- peak_data[[paste0("GR_", cond)]]

  # Build sets
  s10_only <- setdiff(setdiff(s10, mr), gr)
  mr_only  <- setdiff(setdiff(mr, s10), gr)
  gr_only  <- setdiff(setdiff(gr, s10), mr)
  s10_mr   <- intersect(s10, mr)
  s10_gr   <- intersect(s10, gr)
  mr_gr    <- intersect(mr, gr)
  all_three <- Reduce(intersect, list(s10, mr, gr))

  venn.plot <- draw.triple.venn(
    area1 = length(s10),
    area2 = length(mr),
    area3 = length(gr),
    n12 = length(s10_mr),
    n23 = length(mr_gr),
    n13 = length(s10_gr),
    n123 = length(all_three),
    category = c("H3K9acS10p", "MR", "GR"),
    fill = c("gold", "skyblue", "salmon"),
    lty = "blank",
    cex = 2,
    cat.cex = 2,
    cat.col = c("gold", "skyblue", "salmon"),
    main = paste("Venn Diagram -", cond)
  )
}
dev.off()
