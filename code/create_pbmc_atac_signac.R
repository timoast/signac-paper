library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)

plan("multiprocess", workers = 8)
options(future.globals.maxSize = 50 * 1024 ^ 3)

annot <- readRDS("data/pbmc_atac/annotations.rds")
cells <- readLines("data/pbmc_atac/cells.txt")
frags <- "data/pbmc_atac/fragments.bed.gz"
peaks <- read.table(file = "data/pbmc_atac/peaks.bed", sep = "\t", header = TRUE)
peaks <- makeGRangesFromDataFrame(peaks)

start.time <- Sys.time()
fragments <- CreateFragmentObject(
  path = frags,
  cells = cells
)

# quantify
counts <- FeatureMatrix(
  fragments = fragments,
  features = peaks,
  cells = cells
)

# create object
assay <- CreateChromatinAssay(counts = counts, fragments = fragments, annotation = annot)
obj <- CreateSeuratObject(counts = assay, assay = "ATAC")

# QC
obj <- NucleosomeSignal(obj)
obj <- TSSEnrichment(obj)

# LSI
obj <- FindTopFeatures(obj)
obj <- RunTFIDF(obj)
obj <- RunSVD(obj)

# clustering
obj <- FindNeighbors(obj, reduction = "lsi", dims = 2:30)
obj <- FindClusters(obj)

# UMAP
obj <- RunUMAP(obj, reduction = "lsi", dims = 2:30)

elapsed <- as.numeric(Sys.time() - start.time, units = "secs")
writeLines(text = as.character(elapsed), con = "data/pbmc_atac/signac_total_runtime.txt")
