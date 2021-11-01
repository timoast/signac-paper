library(Signac)
library(Seurat)

pbmc <- readRDS("objects/pbmc.rds")
atac.assay <- "ATAC"
method_use <- c(1, 2, 3, 4)
DefaultAssay(pbmc) <- atac.assay
obj <- pbmc

# pbmc multiome
ds_level <- rev(seq(0.2, 1, 0.2))
for (d in ds_level) {
  counts_use <- readRDS(file = paste0("data/pbmc/downsamples/", d, ".rds"))
  obj <- SetAssayData(obj, slot = "counts", assay = atac.assay, new.data = counts_use)
  for (m in method_use) {
    key <- paste(d, m, sep = "_")
    message(key)
    obj <- RunTFIDF(object = obj, assay = atac.assay, method = m)
    time.start <- Sys.time()
    obj <- RunSVD(obj, features = rownames(x = obj))
    elapsed <- as.numeric(Sys.time() - time.start, unit = "secs")
    emb <- Embeddings(object = obj, reduction = "lsi")
    saveRDS(object = emb, file = paste0("data/pbmc/downsamples/lsi_", key, ".rds"))
    write(
      x = paste0(elapsed, "\t", d, "\t", m),
      file = "data/pbmc/downsamples/lsi_runtime.txt",
      append = TRUE
    )
  }
}

# simulated bone marrow
bm_datasets <- c("250", "500", "1000", "2500", "5000")
filepath <- "data/chen/scATAC-benchmarking-master/Synthetic_Data/BoneMarrow_cov"
for (d in bm_datasets) {
  counts_use <- readRDS(file = paste0(filepath, d, "/input/bonemarrow_cov", d, ".rds"))
  obj <- CreateSeuratObject(counts = counts_use, min.cells = -1, min.features = -1, assay = "ATAC")
  for (m in method_use) {
    key <- paste(d, m, sep = "_")
    message(key)
    obj <- RunTFIDF(object = obj, assay = atac.assay, method = m)
    time.start <- Sys.time()
    obj <- RunSVD(obj, features = rownames(x = obj))
    elapsed <- as.numeric(Sys.time() - time.start, unit = "secs")
    emb <- Embeddings(object = obj, reduction = "lsi")
    saveRDS(object = emb, file = paste0("data/chen/embeddings/lsi_", key, ".rds"))
    write(
      x = paste0(elapsed, "\t", d, "\t", m),
      file = "data/chen/embeddings/lsi_runtime.txt",
      append = TRUE
    )
  }
}