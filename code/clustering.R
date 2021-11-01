library(Seurat)
library(Signac)
library(ggplot2)
library(mclust)


pbmc <- readRDS("objects/pbmc.rds")

# sweep clustering parameters
k.param <- seq(5, 50, 5)
dims.param <- seq(10, 50, 5)
celltypes <- pbmc$celltype

cluster.results <- data.frame()
for (i in seq_along(k.param)) {
  for (j in seq_along(dims.param)) {
    pbmc <- FindNeighbors(pbmc, reduction = "lsi", dims = 2:dims.param[[j]], k.param = k.param[[i]])
    pbmc <- FindClusters(pbmc, algorithm = 3, graph.name = "ATAC_snn")
    ari <- adjustedRandIndex(x = celltypes, y = pbmc$seurat_clusters)
    cluster.results <- rbind(cluster.results, data.frame(dims = dims.param[[j]],
                                       k = k.param[[i]],
                                       ari = ari))
  }
}

p <- ggplot(cluster.results, aes(dims, k, fill = ari)) +
  geom_tile() +
  scale_fill_viridis_c() +
  ylab("Number of nearest neighbors (k)") +
  xlab("LSI dimensions (2:n)") + 
  labs(fill = "Adjusted Rand Index") +
  theme_classic()

ggsave(filename = "figures/cluster_param_sweep.png", plot = p, height = 4, width = 7)
