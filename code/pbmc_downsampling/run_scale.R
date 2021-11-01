library(Matrix)

dir.create(file.path("data/pbmc/downsamples/pbmc_scale"), showWarnings = FALSE)
# PBMC dataset
ds_level <- rev(seq(0.2, 1, 0.2))
for (d in ds_level) {
  counts_use <- readRDS(file = paste0("data/pbmc/downsamples/", d, ".rds"))
  writeMM(obj = counts_use, file = "data/pbmc/downsamples/pbmc_scale/counts.mtx")
  peaks <- rownames(counts_use)
  peaks <- gsub("-", "_", peaks)
  write.table(
    x = peaks,
    file = "data/pbmc/downsamples/pbmc_scale/peaks.txt",
    append = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
  barcodes <- colnames(counts_use)
  write.table(
    x = barcodes,
    file = "data/pbmc/downsamples/pbmc_scale/barcodes.txt",
    append = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
  time.start <- Sys.time()
  cmd <- paste0("SCALE.py -d data/pbmc/downsamples/pbmc_scale --min_peaks 1 -o data/pbmc/downsamples/scale_", d)
  system(command = cmd, wait = TRUE, ignore.stderr = FALSE, ignore.stdout = FALSE)
  elapsed <- as.numeric(Sys.time() - time.start, unit = "secs")
  write(
    x = paste0(elapsed, "\t", d),
    file = "data/pbmc/downsamples/scale_runtime.txt",
    append = TRUE
  )
}

dir.create(file.path("data/chen/embeddings/scale"), showWarnings = FALSE)
# Chen dataset
bm_datasets <- c("250", "500", "1000", "2500", "5000")
filepath <- "data/chen/scATAC-benchmarking-master/Synthetic_Data/BoneMarrow_cov"
for (d in bm_datasets) {
  counts_use <- readRDS(file = paste0(filepath, d, "/input/bonemarrow_cov", d, ".rds"))
  writeMM(obj = counts_use, file = "data/chen/embeddings/scale/counts.mtx")
  peaks <- rownames(counts_use)
  write.table(
    x = peaks,
    file = "data/chen/embeddings/scale/peaks.txt",
    append = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
  barcodes <- colnames(counts_use)
  write.table(
    x = barcodes,
    file = "data/chen/embeddings/scale/barcodes.txt",
    append = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
  time.start <- Sys.time()
  cmd <- paste0("SCALE.py -d data/chen/embeddings/scale -o data/chen/embeddings/scale_", d)
  system(command = cmd, wait = TRUE, ignore.stderr = FALSE, ignore.stdout = FALSE)
  elapsed <- as.numeric(Sys.time() - time.start, unit = "secs")
  write(
    x = paste0(elapsed, "\t", d),
    file = "data/chen/embeddings/scale_runtime.txt",
    append = TRUE
  )
}