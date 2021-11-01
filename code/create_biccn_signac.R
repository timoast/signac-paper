library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)

plan("multicore", workers = 8)
options(future.globals.maxSize = +Inf)

annot <- readRDS("data/biccn/annotations.rds")
# load metadata
metadata <- read.table("data/biccn/Supplementary Table 2 - Metatable of nuclei.tsv", sep="\t", skip=1)
rownames(metadata) <- metadata$V1
colnames(metadata) <- c("cell", "sample", "barcode", "logUM", "TSSe", "class", "MajorType", "SubType", "na")
cells <- metadata$cell

frags <- "data/biccn/fragments.bed.gz"
peaks <- read.table(file = "data/biccn/unified_peaks.bed", sep = "\t", header = TRUE)
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

gc()

# QC
obj <- NucleosomeSignal(obj)
obj <- TSSEnrichment(obj)

# LSI
obj <- FindTopFeatures(obj)
obj <- RunTFIDF(obj)
obj <- RunSVD(obj)

# clustering
obj <- FindNeighbors(obj, reduction = "lsi", dims = 2:100)
obj <- FindClusters(obj)

# UMAP
obj <- RunUMAP(obj, reduction = "lsi", dims = 2:100)

elapsed <- as.numeric(Sys.time() - start.time, units = "secs")
writeLines(text = as.character(elapsed), con = "data/biccn/signac_total_runtime.txt")
