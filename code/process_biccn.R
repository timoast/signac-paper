library(Signac)
library(Seurat)
library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)
library(Matrix)
library(future)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 50 * 1024 ^ 3)


annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "mm10"

fragment.counts <- CountFragments(
  fragments = "data/biccn/fragments.sort.bed.gz"
)

cells <- fragment.counts[fragment.counts$frequency_count > 1500, ]$CB

rm(fragment.counts)
gc()

frags <- CreateFragmentObject(
  path = "data/biccn/fragments.sort.bed.gz",
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

rsums <- rowSums(counts)
csums <- colSums(counts)
counts <- counts[rsums > 100, csums > 1000]
gc()

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
  project = "BICCN"
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

# Clustering
# biccn <- FindNeighbors(biccn, reduction = "lsi", dims = 2:100)
# biccn <- FindClusters(biccn, algorithm = 3)

# Add region information

regions <- c(
  "1A" = "MOs-1",
  "1B" = "ORB",
  "1C" = "MOB",
  "2A" = "LIMBIC-1",
  "2B" = "MOs-2",
  "2C" = "MOp-1",
  "2D" = "PIR-1",
  "2E" = "AON",
  "3A" = "LIMBIC-2",
  "3B" = "MOs-3",
  "3C" = "MOp-2",
  "3D" = "AI",
  "3E" = "PIR-2",
  "3F" = "ACB-1",
  "4A" = "LIMBIC-3",
  "4B" = "MOp-3",
  "4C" = "SSp-1",
  "4D" = "CP-1",
  "4E" = "ACB-2",
  "4F" = "PIR-3",
  "4G" = "LSX-1",
  "4H" = "PAL-1",
  "5A" = "LIMBIC-4",
  "5B" = "SSp-2",
  "5C" = "SSs-1",
  "5D" = "MOp-4",
  "5E" = "CP-2",
  "5F" = "ACB-3",
  "5G" = "PIR-4",
  "5H" = "PAL-2",
  "5J" = "LSX-2",
  "6A" = "LIMBIC-5",
  "6B" = "SSp-3",
  "6C" = "SSs-2",
  "6D" = "PIR-5",
  "7B" = "SSp-4",
  "8B" = "SSp-5",
  "8E" = "CA-1",
  "8J" = "DG-1",
  "9A" = "Other",
  "9B" = "Other",
  "9D" = "Other",
  "9H" = "CA-2",
  "9J" = "DG-2",
  "10C" = "Other",
  "10G" = "Other",
  "11B" = "Other",
  "10E" = "CA-3",
  "10F" = "DG-3",
  "11E" = "CA-4",
  "11F" = "DG-4"
)

# add region
biccn$region <- regions[biccn$orig.ident]

# add coarse-level region
biccn$broad_region <- sapply(biccn$region, function(x) unlist(strsplit(x, split = "-", fixed = TRUE))[1])

saveRDS(object = biccn, file = "objects/biccn.rds")
