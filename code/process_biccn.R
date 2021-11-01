library(Signac)
library(Seurat)
library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)
library(Matrix)


annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "mm10"

# load metadata
metadata <- read.table("data/biccn/Supplementary Table 2 - Metatable of nuclei.tsv", sep="\t", skip=1)
rownames(metadata) <- metadata$V1
colnames(metadata) <- c("cell", "sample", "barcode", "logUM", "TSSe", "class", "MajorType", "SubType", "na")
cells <- metadata$cell

frags <- CreateFragmentObject(
  path = "data/biccn/fragments.bed.gz",
  cells = cells,
  validate.fragments = FALSE
)

# load peaks
unified.peaks <- read.table("data/biccn/unified_peaks.bed", sep = "\t", header = TRUE)
unified.peaks <- makeGRangesFromDataFrame(unified.peaks)

# quantify
counts <- FeatureMatrix(
  fragments = frags,
  features = unified.peaks,
  cells = cells
)

csums <- colSums(counts)
rsums <- rowSums(counts)
counts <- counts[rsums > 100, csums > 100]

# create seurat object
biccn_assay <- CreateChromatinAssay(
  counts = counts,
  fragments = frags,
  genome = "mm10",
  annotation = annotations
)

biccn <- CreateSeuratObject(
  counts = biccn_assay,
  assay = "ATAC",
  project = "BICCN",
  meta.data = metadata
)

rm(biccn_assay)
gc()

# QC
biccn <- NucleosomeSignal(biccn, n = 1e9)
biccn <- TSSEnrichment(biccn)

# Dim reduc
biccn <- FindTopFeatures(biccn)
biccn <- RunTFIDF(biccn)
biccn <- RunSVD(biccn, n = 100)
biccn <- RunUMAP(biccn, reduction = "lsi", dims = 2:100)

saveRDS(object = biccn, file = "objects/biccn.rds")
