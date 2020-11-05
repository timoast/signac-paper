library(Signac)
library(Seurat)
library(DropletUtils)
library(ggplot2)
library(RANN)

set.seed(1234)

atac.assay <- "ATAC"
pbmc <- readRDS("objects/pbmc.rds")
counts <- GetAssayData(pbmc, slot = "counts", assay = atac.assay)

# downsample counts
ds_level <- rev(seq(0.2, 1, 0.2))
method_use <- c(1, 2)

dimplot_list <- list()
hist_list <- list()
purity_list <- list()
order.use <- colnames(pbmc)

# use cell types defined by RNA
clustering.use <- "celltype"
clusters <- pbmc[[clustering.use]][[1]]

knn_purity <- function(embeddings, clusters, k = 100) {
  nn <- nn2(data = embeddings, k = k + 1)$nn.idx[, 2:k] # remove self-neighbor
  # find percentage of neighbors that are of the same cluster
  nn_purity <- vector(mode = "numeric", length = length(x = clusters))
  for (i in seq_len(length.out = nrow(x = nn))) {
    nn_purity[i] <- sum(clusters[nn[i, ]] == clusters[i]) / k
  }
  return(nn_purity)
}

DefaultAssay(pbmc) <- atac.assay
obj <- pbmc
for (d in ds_level) {
  counts_use <- downsampleMatrix(x = counts, prop = d)
  obj <- SetAssayData(obj, slot = "counts", assay = atac.assay, new.data = counts_use)
  for (m in method_use) {
    key <- paste(d, m, sep = "_")
    message(key)
    obj <- RunTFIDF(object = obj, assay = atac.assay, method = m)
    nonzero_vals <- GetAssayData(object = obj, assay = atac.assay, slot = "data")@x
    nz_hist <- ggplot(data = data.frame(nzv = nonzero_vals),
                      mapping = aes(x = nzv)) +
      geom_histogram(bins = 100) +
      theme_bw() +
      xlab("Non-zero TF-IDF values")
    obj <- RunSVD(obj)
    obj <- RunUMAP(obj, reduction = "lsi", dims = 2:30, reduction.name = "dsumap")
    dp <- DimPlot(obj, group.by = clustering.use, reduction = "dsumap")
    emb <- Embeddings(object = obj, reduction = "lsi")[order.use, 2:30]
    nn <- knn_purity(embeddings = emb, clusters = clusters, k = 100)
    dimplot_list[[key]] <- dp
    hist_list[[key]] <- nz_hist
    purity_list[[key]] <- nn
  }
}

saveRDS(object = dimplot_list, file = "figures/dimplot_downsample_pbmc.rds")
saveRDS(object = hist_list, file = "figures/hist_downsample_pbmc.rds")
saveRDS(object = purity_list, file = "figures/knn_purity_pbmc.rds")
