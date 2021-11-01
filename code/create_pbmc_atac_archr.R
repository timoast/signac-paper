library(ArchR)
library(GenomicRanges)

cells <- readLines("data/pbmc_atac/cells.txt")
frags <- "data/pbmc_atac/fragments.bed.gz"
peaks <- read.table(file = "data/pbmc_atac/peaks.bed", sep = "\t", header = TRUE)
peaks <- makeGRangesFromDataFrame(peaks)

addArchRThreads(threads = 8)
addArchRGenome("hg19")

start.time <- Sys.time()
ArrowFiles <- createArrowFiles(
  inputFiles = frags,
  sampleNames = "PBMC",
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
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "PeakMatrix", sampleCellsPre = NULL, force = TRUE)
proj <- addClusters(input = proj, force = TRUE)
proj <- addUMAP(ArchRProj = proj, force = TRUE)

elapsed <- as.numeric(Sys.time() - start.time, units = "secs")
writeLines(text = as.character(elapsed), con = "data/pbmc_atac/archr_total_runtime.txt")
