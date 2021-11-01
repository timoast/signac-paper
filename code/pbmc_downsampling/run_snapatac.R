library(Signac)
library(SnapATAC)

# PBMC dataset
ds_level <- rev(seq(0.2, 1, 0.2))
for (d in ds_level) {
  counts_use <- readRDS(file = paste0("data/pbmc/downsamples/", d, ".rds"))
  counts_use <- t(x = counts_use)
  snap <- createSnapFromBmat(
    mat = counts_use,
    barcodes = rownames(x = counts_use),
    bins = StringToGRanges(regions = colnames(x = counts_use))
  )
  snap <- makeBinary(snap, mat = "bmat")
  time.start <- Sys.time()
  snap <- runDiffusionMaps(
    obj = snap,
    input.mat = "bmat",
    num.eigs = 50
  )
  elapsed <- as.numeric(Sys.time() - time.start, unit = "secs")
  reducedMatrix <- snap@smat@dmat
  saveRDS(object = reducedMatrix, file = paste0("data/pbmc/downsamples/snapatac_", d, ".rds"))
  write(
    x = paste0(elapsed, "\t", d),
    file = "data/pbmc/downsamples/snapatac_runtime.txt",
    append = TRUE
  )
}

# Chen dataset
bm_datasets <- c("250", "500", "1000", "2500", "5000")
filepath <- "data/chen/scATAC-benchmarking-master/Synthetic_Data/BoneMarrow_cov"
for (d in bm_datasets) {
  counts_use <- readRDS(file = paste0(filepath, d, "/input/bonemarrow_cov", d, ".rds"))
  counts_use <- t(x = counts_use)
  snap <- createSnapFromBmat(
    mat = counts_use,
    barcodes = rownames(x = counts_use),
    bins = StringToGRanges(regions = colnames(x = counts_use), sep = c("_", "_"))
  )
  snap <- makeBinary(snap, mat = "bmat")
  time.start <- Sys.time()
  snap <- runDiffusionMaps(
    obj = snap,
    input.mat = "bmat",
    num.eigs = 50
  )
  elapsed <- as.numeric(Sys.time() - time.start, unit = "secs")
  reducedMatrix <- snap@smat@dmat
  saveRDS(object = reducedMatrix, file = paste0("data/chen/embeddings/snapatac_", d, ".rds"))
  write(
    x = paste0(elapsed, "\t", d),
    file = "data/chen/embeddings/snapatac_runtime.txt",
    append = TRUE
  )
}