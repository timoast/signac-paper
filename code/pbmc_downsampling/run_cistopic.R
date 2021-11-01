library(Signac)
library(cisTopic)

# PBMC dataset
ds_level <- rev(seq(0.2, 1, 0.2))
for (d in ds_level) {
  counts_use <- readRDS(file = paste0("data/pbmc/downsamples/", d, ".rds"))
  rownames(counts_use) <- GRangesToString(StringToGRanges(rownames(counts_use)), sep = c(":", "-"))
  cisTopicObj <- createcisTopicObject(count.matrix = counts_use, project.name = 'PBMC')
  
  # CGS model
  time.start <- Sys.time()
  cgs <- runCGSModels(object = cisTopicObj)
  elapsed.cgs <- as.numeric(Sys.time() - time.start, unit = "secs")
  cgs <- selectModel(object = cgs, type = "maximum")
  
  # WarpLDA model
  time.start <- Sys.time()
  wrp <- runWarpLDAModels(object = cisTopicObj)
  elapsed.warp <- as.numeric(Sys.time() - time.start, unit = "secs")
  wrp <- selectModel(object = wrp, type = "derivative")
  
  # extract coordinates
  dimreduc_cgs <- cgs@selected.model$document_expects
  dimreduc_wrp <- wrp@selected.model$document_expects
  
  # save
  saveRDS(object = dimreduc_cgs, file = paste0("data/pbmc/downsamples/cistopic_cgs_", d, ".rds"))
  saveRDS(object = dimreduc_wrp, file = paste0("data/pbmc/downsamples/cistopic_warp_", d, ".rds"))
  
  write(
    x = paste0(elapsed.cgs, "\t", d),
    file = "data/pbmc/downsamples/cistopic_cgs_runtime.txt",
    append = TRUE
  )
  write(
    x = paste0(elapsed.warp, "\t", d),
    file = "data/pbmc/downsamples/cistopic_warp_runtime.txt",
    append = TRUE
  )
}

# Chen dataset
bm_datasets <- c("250", "500", "1000", "2500", "5000")
filepath <- "data/chen/scATAC-benchmarking-master/Synthetic_Data/BoneMarrow_cov"
for (d in bm_datasets) {
  counts_use <- readRDS(file = paste0(filepath, d, "/input/bonemarrow_cov", d, ".rds"))
  rownames(counts_use) <- GRangesToString(StringToGRanges(rownames(counts_use), sep = c("_", "_")), sep = c(":", "-"))
  cisTopicObj <- createcisTopicObject(count.matrix = counts_use, project.name = 'Simulated')
  
  # CGS model
  time.start <- Sys.time()
  cgs <- runCGSModels(object = cisTopicObj)
  elapsed.cgs <- as.numeric(Sys.time() - time.start, unit = "secs")
  cgs <- selectModel(object = cgs, type = "maximum")
  
  # WarpLDA model
  time.start <- Sys.time()
  wrp <- runWarpLDAModels(object = cisTopicObj)
  elapsed.warp <- as.numeric(Sys.time() - time.start, unit = "secs")
  wrp <- selectModel(object = wrp, type = "derivative")
  
  # extract coordinates
  dimreduc_cgs <- cgs@selected.model$document_expects
  dimreduc_wrp <- wrp@selected.model$document_expects
  
  # save
  saveRDS(object = dimreduc_cgs, file = paste0("data/chen/embeddings/cistopic_cgs_", d, ".rds"))
  saveRDS(object = dimreduc_wrp, file = paste0("data/chen/embeddings/cistopic_warp_", d, ".rds"))
  
  write(
    x = paste0(elapsed.cgs, "\t", d),
    file = "data/chen/embeddings/cistopic_cgs_runtime.txt",
    append = TRUE
  )
  write(
    x = paste0(elapsed.warp, "\t", d),
    file = "data/chen/embeddings/cistopic_warp_runtime.txt",
    append = TRUE
  )
}