library(Signac)
library(Seurat)
library(future)
plan(strategy = "multiprocess", workers = 8)
options(future.globals.maxSize = 50 * 1024 ^ 3)

pbmc <- readRDS('objects/pbmc.rds')
DefaultAssay(pbmc) <- "ATAC"

# link peaks to genes
pbmc <- LinkPeaks(
  object = pbmc,
  peak.assay = "ATAC",
  expression.assay = "SCT"
)

saveRDS(object = Links(pbmc), file = "objects/pbmc_links.rds")
