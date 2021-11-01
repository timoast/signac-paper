library(ArchR)
library(GenomicRanges)

# load metadata
metadata <- read.table("data/biccn/Supplementary Table 2 - Metatable of nuclei.tsv", sep="\t", skip=1)
rownames(metadata) <- metadata$V1
colnames(metadata) <- c("cell", "sample", "barcode", "logUM", "TSSe", "class", "MajorType", "SubType", "na")
cells <- metadata$cell

frags <- "data/biccn/fragments.bed.gz"
peaks <- read.table(file = "data/biccn/unified_peaks.bed", sep = "\t", header = TRUE)
peaks <- makeGRangesFromDataFrame(peaks)
# remove chrM
peaks <- peaks[seqnames(peaks) != "chrM"]
message("Using ", length(peaks), " peaks")

addArchRThreads(threads = 8)
addArchRGenome("mm10")

start.time <- Sys.time()
ArrowFiles <- createArrowFiles(
  inputFiles = frags,
  sampleNames = "BICCN",
  validBarcodes = cells,
  force = TRUE,
  minFrags = 1,
  addGeneScoreMat = FALSE,
  addTileMat = FALSE
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  copyArrows = TRUE,
  showLogo = FALSE
)

proj <- addPeakSet(ArchRProj = proj, peakSet = peaks, force = TRUE)
proj <- addPeakMatrix(ArchRProj = proj, force = TRUE)
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "PeakMatrix", force = TRUE)
proj <- addClusters(input = proj, force = TRUE, dimsToUse = 2:100)
proj <- addUMAP(ArchRProj = proj, force = TRUE, dimsToUse = 2:100)

elapsed <- as.numeric(Sys.time() - start.time, units = "secs")
writeLines(text = as.character(elapsed), con = "data/biccn/archr_total_runtime.txt")
