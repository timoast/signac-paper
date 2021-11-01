library(Signac)
library(Seurat)
library(DropletUtils)

set.seed(1234)

atac.assay <- "ATAC"
pbmc <- readRDS("objects/pbmc.rds")
counts <- GetAssayData(pbmc, slot = "counts", assay = atac.assay)

# downsample counts
ds_level <- rev(seq(0.2, 1, 0.2))

for (d in ds_level) {
  counts_use <- downsampleMatrix(x = counts, prop = d)
  saveRDS(object = counts_use, file = paste0("data/pbmc/downsamples/", d, ".rds"))
}
