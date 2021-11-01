library(Signac)
library(Seurat)
library(GenomicRanges)

frags <- "data/pbmc_atac/fragments.bed.gz"

fragment.counts <- CountFragments(frags)
cells.use <- fragment.counts[fragment.counts$frequency_count > 1000, "CB"]

fragments <- CreateFragmentObject(
  path = frags,
  cells = cells.use,
  validate.fragments = FALSE
)

peaks <- CallPeaks(fragments, macs2.path = "/home/stuartt/miniconda3/envs/signac/bin/macs2")
peaks <- subsetByOverlaps(peaks, blacklist_hg19, invert = TRUE)

counts <- FeatureMatrix(
  fragments = fragments,
  features = peaks,
  cells = cells.use
)

pbmc <- CreateSeuratObject(
  counts = CreateChromatinAssay(
    counts = counts,
    fragments = fragments
  ),
  assay = "ATAC"
)

pbmc <- pbmc[, pbmc$nCount_ATAC > 1000]

peaks <- granges(pbmc)
peaks <- as.data.frame(peaks)
write.table(x = peaks, file = "data/pbmc_atac/peaks.bed", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
writeLines(text = colnames(x = pbmc), con = "data/pbmc_atac/cells.txt")

# cluster and make UMAP
pbmc <- FindTopFeatures(pbmc, min.cutoff = 10)
pbmc <- RunTFIDF(pbmc)
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(pbmc, reduction = "lsi", dims = 2:30)
pbmc <- FindNeighbors(pbmc, reduction = "lsi", dims = 2:30)
pbmc <- FindClusters(pbmc, algorithm = 3, resolution = 0.5)

saveRDS(object = pbmc, file = "objects/pbmc_atac.rds")
