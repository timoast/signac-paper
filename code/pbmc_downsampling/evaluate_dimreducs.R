library(Signac)
library(Seurat)
library(RANN)
library(ggplot2)
library(cluster)
library(dplyr)
library(patchwork)
library(paletteer)
library(SeuratDisk)

set.seed(1234)

pbmc <- readRDS("objects/pbmc.rds")
atac.assay <- "ATAC"
methods_keep <- c("LSI (Cusanovich2018)", "SnapATAC", "cisTopic CGS", "cisTopic Warp", "SCALE", "LSI (log-TF)", "LSI (Signac)")
colors.use <- paletteer_d("ggthemes::Tableau_10")
colors.use <- rev(colors.use[1:length(methods_keep)])

######## Data loading #########

read_lsi <- function(method, path ="data/pbmc/downsamples/") {
  methodstr <- paste0(as.character(method), ".rds")
  lsi <- paste0(path, c(
    paste0("lsi_1_", methodstr),
    paste0("lsi_0.8_", methodstr),
    paste0("lsi_0.6_", methodstr),
    paste0("lsi_0.4_", methodstr),
    paste0("lsi_0.2_", methodstr)
  ))
  lsi_obj <- lapply(X = lsi, readRDS)
  return(lsi_obj)
}

lsi_1 <- read_lsi(method = 1)
lsi_2 <- read_lsi(method = 2)
lsi_3 <- read_lsi(method = 3)
lsi_4 <- read_lsi(method = 4)

pbmc_ds <- c(1, 0.8, 0.6, 0.4, 0.2)

snap <- lapply(X = paste0("data/pbmc/downsamples/snapatac_", pbmc_ds, ".rds"), FUN = readRDS)
ct_cgs <- lapply(X = paste0("data/pbmc/downsamples/cistopic_cgs_", pbmc_ds, ".rds"), FUN = readRDS)
ct_warp <- lapply(X = paste0("data/pbmc/downsamples/cistopic_warp_", pbmc_ds, ".rds"), FUN = readRDS)

ct_cgs <- lapply(ct_cgs, t)
ct_warp <- lapply(ct_warp, t)

# convert h5ad to h5seurat
scale_pbmc_path <- lapply(
  X = pbmc_ds,
  function(x) {
    Convert(
      source = paste0("data/pbmc/downsamples/scale_", x, "/adata.h5ad"),
      dest = paste0("data/pbmc/downsamples/scale_", x, "/adata.h5seurat"),
      overwrite = TRUE
    )
  }
)

# load h5seurat
scale_pbmc_obj <- lapply(X = scale_pbmc_path, FUN = LoadH5Seurat)

# get embeddings
scale_pbmc <- lapply(X = scale_pbmc_obj, FUN = function(x) {
  Embeddings(x[["latent"]])
})

######## Determine dimensions to use ##########

seqdepth_pbmc <- pbmc$nCount_ATAC

# dim 1
lsi1_depth <- lapply(lsi_1, function(x) {
  which(abs(cor(x, seqdepth_pbmc)) > 0.9)
})
# no dims
lsi2_depth <- lapply(lsi_2, function(x) {
  which(abs(cor(x, seqdepth_pbmc)) > 0.9)
})
# no dims
lsi3_depth <- lapply(lsi_3, function(x) {
  which(abs(cor(x, seqdepth_pbmc)) > 0.9)
})
# no dims
lsi4_depth <- lapply(lsi_4, function(x) {
  which(abs(cor(x, seqdepth_pbmc)) > 0.9)
})
# no dims
cgs_depth <- lapply(ct_cgs, function(x) {
  which(abs(cor(x, seqdepth_pbmc)) > 0.9)
})
# no dims
warp_depth <- lapply(ct_warp, function(x) {
  which(abs(cor(x, seqdepth_pbmc)) > 0.9)
})
# dim 2
snap_depth <- lapply(snap, function(x) {
  which(abs(cor(x, seqdepth_pbmc)) > 0.9)
})
# no dims
scale_depth <- lapply(scale_pbmc, function(x) {
  which(abs(cor(x, seqdepth_pbmc)) > 0.9)
})

######## KNN #########

# use cell types defined by RNA
clustering.use <- "celltype"
clusters <- pbmc[[clustering.use]][[1]]

# first define neighbor graph using the RNA assay
k <- 100
rna.emb <- Embeddings(pbmc[["pca"]])
rna.nn <- nn2(data = rna.emb, k = k + 1)$nn.idx[, 2:k]

knn_purity <- function(embeddings, dims, clusters, rna.nn, k = 100) {
  nn <- nn2(data = embeddings[, dims], k = k + 1)$nn.idx[, 2:k]
  nn_purity <- vector(mode = "numeric", length = length(x = clusters))
  for (i in seq_len(length.out = nrow(x = nn))) {
    nn_purity[i] <- sum(clusters[nn[i, ]] == clusters[i]) / k
  }
  return(nn_purity)
}

get_knn_df <- function(emb_list, dims, clusters, rna_nn, method, ds_list, k) {
  # compute KNN purity for each dimension reduction
  knn_df <- data.frame()
  for (i in seq_along(along.with = emb_list)) {
    knn <- knn_purity(embeddings = emb_list[[i]], dims = dims, clusters = clusters, rna.nn = rna_nn, k = k)
    ds <- ds_list[[i]]
    kd <- data.frame(
      purity = knn,
      downsample = ds, method = method,
      celltype = clusters
    )
    knn_df <- rbind(knn_df, kd)
  }
  return(knn_df)
}

k <- 100
knn_lsi1 <- get_knn_df(
  emb_list = lsi_1,
  dims = 2:20,
  ds_list = pbmc_ds,
  clusters = clusters,
  rna_nn = rna.nn,
  method = "LSI (Signac)",
  k = k
)

knn_lsi2 <- get_knn_df(
  emb_list = lsi_2,
  dims = 1:20,
  ds_list = pbmc_ds,
  clusters = clusters,
  rna_nn = rna.nn,
  method = "LSI (Cusanovich2018)",
  k = k
)

knn_lsi3 <- get_knn_df(
  emb_list = lsi_3,
  dims = 1:20,
  ds_list = pbmc_ds,
  clusters = clusters,
  rna_nn = rna.nn,
  method = "LSI (log-TF)",
  k = k
)

knn_lsi4 <- get_knn_df(
  emb_list = lsi_4,
  dims = 1:20,
  ds_list = pbmc_ds,
  clusters = clusters,
  rna_nn = rna.nn,
  method = "LSI (Cellranger)",
  k = k
)

knn_snap <- get_knn_df(
  emb_list = snap,
  dims = c(1, 3:20),
  ds_list = pbmc_ds,
  clusters = clusters,
  rna_nn = rna.nn,
  method = "SnapATAC",
  k = k
)

knn_ct_cgs <- get_knn_df(
  emb_list = ct_cgs,
  dims = 1:20,
  ds_list = pbmc_ds,
  clusters = clusters,
  rna_nn = rna.nn,
  method = "cisTopic CGS",
  k = k
)

knn_ct_warp <- get_knn_df(
  emb_list = ct_warp,
  dims = 1:20,
  ds_list = pbmc_ds,
  clusters = clusters,
  rna_nn = rna.nn,
  method = "cisTopic Warp",
  k = k
)

knn_scale <- get_knn_df(
  emb_list = scale_pbmc,
  dims = 1:10,
  ds_list = pbmc_ds,
  clusters = clusters,
  rna_nn = rna.nn,
  method = "SCALE",
  k = k
)

knn_df <- rbind(knn_lsi1, knn_lsi2, knn_lsi3, knn_lsi4, knn_snap, knn_ct_cgs, knn_ct_warp, knn_scale)
knn_df$downsample <- factor(knn_df$downsample, levels = rev(pbmc_ds))

knn_plot <- knn_df[knn_df$method %in% methods_keep, ]
knn_plot$method <- factor(knn_plot$method, levels = methods_keep)

knn_plot<- knn_plot %>% 
  group_by(celltype, method, downsample) %>% 
  mutate(mn = mean(purity)) %>% 
  ungroup()

knn_plot <- knn_plot[, c("celltype", "method", "downsample", "mn")]
knn_plot <- unique(knn_plot)

p2 <- ggplot(knn_plot, aes(x = downsample, y = mn, fill = method)) +
  geom_boxplot(outlier.size = 0.1) +
  theme_bw() +
  ylab("Mean kNN celltype purity") +
  xlab("Fraction of counts retained") +
  scale_fill_manual(values = colors.use)

# test choice of K
knn_df$k <- k

k <- 150
knn_lsi1 <- get_knn_df(
  emb_list = lsi_1,
  dims = 2:20,
  ds_list = pbmc_ds,
  clusters = clusters,
  rna_nn = rna.nn,
  method = "LSI (Signac)",
  k = k
)

knn_lsi2 <- get_knn_df(
  emb_list = lsi_2,
  dims = 1:20,
  ds_list = pbmc_ds,
  clusters = clusters,
  rna_nn = rna.nn,
  method = "LSI (Cusanovich2018)",
  k = k
)

knn_lsi3 <- get_knn_df(
  emb_list = lsi_3,
  dims = 1:20,
  ds_list = pbmc_ds,
  clusters = clusters,
  rna_nn = rna.nn,
  method = "LSI (log-TF)",
  k = k
)

knn_lsi4 <- get_knn_df(
  emb_list = lsi_4,
  dims = 1:20,
  ds_list = pbmc_ds,
  clusters = clusters,
  rna_nn = rna.nn,
  method = "LSI (Cellranger)",
  k = k
)

knn_snap <- get_knn_df(
  emb_list = snap,
  dims = c(1, 3:20),
  ds_list = pbmc_ds,
  clusters = clusters,
  rna_nn = rna.nn,
  method = "SnapATAC",
  k = k
)

knn_ct_cgs <- get_knn_df(
  emb_list = ct_cgs,
  dims = 1:20,
  ds_list = pbmc_ds,
  clusters = clusters,
  rna_nn = rna.nn,
  method = "cisTopic CGS",
  k = k
)

knn_ct_warp <- get_knn_df(
  emb_list = ct_warp,
  dims = 1:20,
  ds_list = pbmc_ds,
  clusters = clusters,
  rna_nn = rna.nn,
  method = "cisTopic Warp",
  k = k
)

knn_scale <- get_knn_df(
  emb_list = scale_pbmc,
  dims = 1:10,
  ds_list = pbmc_ds,
  clusters = clusters,
  rna_nn = rna.nn,
  method = "SCALE",
  k = k
)

knn_df2 <- rbind(knn_lsi1, knn_lsi2, knn_lsi3, knn_lsi4, knn_snap, knn_ct_cgs, knn_ct_warp, knn_scale)
knn_df2$downsample <- factor(knn_df2$downsample, levels = rev(pbmc_ds))
knn_df2$k <- k

knn_df <- rbind(knn_df, knn_df2)

k <- 50
knn_lsi1 <- get_knn_df(
  emb_list = lsi_1,
  dims = 2:20,
  ds_list = pbmc_ds,
  clusters = clusters,
  rna_nn = rna.nn,
  method = "LSI (Signac)",
  k = k
)

knn_lsi2 <- get_knn_df(
  emb_list = lsi_2,
  dims = 1:20,
  ds_list = pbmc_ds,
  clusters = clusters,
  rna_nn = rna.nn,
  method = "LSI (Cusanovich2018)",
  k = k
)

knn_lsi3 <- get_knn_df(
  emb_list = lsi_3,
  dims = 1:20,
  ds_list = pbmc_ds,
  clusters = clusters,
  rna_nn = rna.nn,
  method = "LSI (log-TF)",
  k = k
)

knn_lsi4 <- get_knn_df(
  emb_list = lsi_4,
  dims = 1:20,
  ds_list = pbmc_ds,
  clusters = clusters,
  rna_nn = rna.nn,
  method = "LSI (Cellranger)",
  k = k
)

knn_snap <- get_knn_df(
  emb_list = snap,
  dims = c(1, 3:20),
  ds_list = pbmc_ds,
  clusters = clusters,
  rna_nn = rna.nn,
  method = "SnapATAC",
  k = k
)

knn_ct_cgs <- get_knn_df(
  emb_list = ct_cgs,
  dims = 1:20,
  ds_list = pbmc_ds,
  clusters = clusters,
  rna_nn = rna.nn,
  method = "cisTopic CGS",
  k = k
)

knn_ct_warp <- get_knn_df(
  emb_list = ct_warp,
  dims = 1:20,
  ds_list = pbmc_ds,
  clusters = clusters,
  rna_nn = rna.nn,
  method = "cisTopic Warp",
  k = k
)

knn_scale <- get_knn_df(
  emb_list = scale_pbmc,
  dims = 1:10,
  ds_list = pbmc_ds,
  clusters = clusters,
  rna_nn = rna.nn,
  method = "SCALE",
  k = k
)

knn_df2 <- rbind(knn_lsi1, knn_lsi2, knn_lsi3, knn_lsi4, knn_snap, knn_ct_cgs, knn_ct_warp, knn_scale)
knn_df2$downsample <- factor(knn_df2$downsample, levels = rev(pbmc_ds))
knn_df2$k <- k

knn_df <- rbind(knn_df, knn_df2)

k <- 10
knn_lsi1 <- get_knn_df(
  emb_list = lsi_1,
  dims = 2:20,
  ds_list = pbmc_ds,
  clusters = clusters,
  rna_nn = rna.nn,
  method = "LSI (Signac)",
  k = k
)

knn_lsi2 <- get_knn_df(
  emb_list = lsi_2,
  dims = 1:20,
  ds_list = pbmc_ds,
  clusters = clusters,
  rna_nn = rna.nn,
  method = "LSI (Cusanovich2018)",
  k = k
)

knn_lsi3 <- get_knn_df(
  emb_list = lsi_3,
  dims = 1:20,
  ds_list = pbmc_ds,
  clusters = clusters,
  rna_nn = rna.nn,
  method = "LSI (log-TF)",
  k = k
)

knn_lsi4 <- get_knn_df(
  emb_list = lsi_4,
  dims = 1:20,
  ds_list = pbmc_ds,
  clusters = clusters,
  rna_nn = rna.nn,
  method = "LSI (Cellranger)",
  k = k
)

knn_snap <- get_knn_df(
  emb_list = snap,
  dims = c(1, 3:20),
  ds_list = pbmc_ds,
  clusters = clusters,
  rna_nn = rna.nn,
  method = "SnapATAC",
  k = k
)

knn_ct_cgs <- get_knn_df(
  emb_list = ct_cgs,
  dims = 1:20,
  ds_list = pbmc_ds,
  clusters = clusters,
  rna_nn = rna.nn,
  method = "cisTopic CGS",
  k = k
)

knn_ct_warp <- get_knn_df(
  emb_list = ct_warp,
  dims = 1:20,
  ds_list = pbmc_ds,
  clusters = clusters,
  rna_nn = rna.nn,
  method = "cisTopic Warp",
  k = k
)

knn_scale <- get_knn_df(
  emb_list = scale_pbmc,
  dims = 1:10,
  ds_list = pbmc_ds,
  clusters = clusters,
  rna_nn = rna.nn,
  method = "SCALE",
  k = k
)

knn_df2 <- rbind(knn_lsi1, knn_lsi2, knn_lsi3, knn_lsi4, knn_snap, knn_ct_cgs, knn_ct_warp, knn_scale)
knn_df2$downsample <- factor(knn_df2$downsample, levels = rev(pbmc_ds))
knn_df2$k <- k

knn_df <- rbind(knn_df, knn_df2)

knn_df$k <- paste0("k=", as.character(knn_df$k))
knn_df$k <- factor(knn_df$k, levels = c("k=10", "k=50", "k=100", "k=150"))

knn_plot <- knn_df[knn_df$method %in% methods_keep, ]
knn_plot$method <- factor(knn_plot$method, levels = methods_keep)

knn_plot<- knn_plot %>% 
  group_by(celltype, method, downsample, k) %>% 
  mutate(mn = mean(purity)) %>% 
  ungroup()

knn_plot <- knn_plot[, c("celltype", "method", "downsample", "mn", "k")]
knn_plot <- unique(knn_plot)

knn_sensitivity <- ggplot(knn_plot, aes(x = downsample, y = mn, fill = method)) +
  geom_boxplot(outlier.size = 0.1) +
  theme_bw() +
  facet_wrap(~k, ncol = 1) + 
  ylab("Mean kNN celltype purity") +
  xlab("Fraction of counts retained") +
  scale_fill_manual(values = colors.use)

ggsave(filename = "figures/knn_sensitivity.png", plot = knn_sensitivity, height = 8, width = 10, dpi = 300)

###### UMAP ########

umaps_lsi_1 <- lapply(lsi_1, function(x) RunUMAP(x[, 2:20]))
umaps_lsi_2 <- lapply(lsi_2, function(x) RunUMAP(x[, 1:20]))
umaps_lsi_3 <- lapply(lsi_3, function(x) RunUMAP(x[, 1:20]))
umaps_lsi_4 <- lapply(lsi_4, function(x) RunUMAP(x[, 1:20]))
umaps_snap <- lapply(snap, function(x) {
  rownames(x) <- colnames(pbmc)
  dr <- RunUMAP(x[, c(1, 3:20)])
  dr
  }
)
umaps_ct_cgs <- lapply(ct_cgs, function(x) {
  rownames(x) <- colnames(pbmc)
  dr <- RunUMAP(x[, 1:20])
  dr
  }
)
umaps_ct_warp <- lapply(ct_warp, function(x) {
  rownames(x) <- colnames(pbmc)
  dr <- RunUMAP(x[, 1:20])
  dr
}
)
umaps_scale <- lapply(scale_pbmc, function(x) RunUMAP(x[, 1:10]))

######### Runtimes ##########

runtime_lsi <- read.table("data/pbmc/downsamples/lsi_runtime.txt", sep = "\t")
colnames(runtime_lsi) <- c("Seconds", "Downsample", "Method")
runtime_lsi[runtime_lsi$Method == 1, "Method"] <- "LSI (Signac)"
runtime_lsi[runtime_lsi$Method == 2, "Method"] <- "LSI (Cusanovich2018)"
runtime_lsi[runtime_lsi$Method == 3, "Method"] <- "LSI (log-TF)"
runtime_lsi[runtime_lsi$Method == 4, "Method"] <- "LSI (Cellranger)"

runtime_snap <- read.table("data/pbmc/downsamples/snapatac_runtime.txt", sep = "\t")
colnames(runtime_snap) <- c("Seconds", "Downsample")
runtime_snap$Method <- "SnapATAC"

runtime_cistopic_cg <- read.table("data/pbmc/downsamples/cistopic_cgs_runtime.txt", sep = "\t")
colnames(runtime_cistopic_cg) <- c("Seconds", "Downsample")
runtime_cistopic_cg$Method <- "cisTopic CGS"

runtime_cistopic_warp <- read.table("data/pbmc/downsamples/cistopic_warp_runtime.txt", sep = "\t")
colnames(runtime_cistopic_warp) <- c("Seconds", "Downsample")
runtime_cistopic_warp$Method <- "cisTopic Warp"

runtime_scale <- read.table("data/pbmc/downsamples/scale_runtime.txt", sep = "\t")
colnames(runtime_scale) <- c("Seconds", "Downsample")
runtime_scale$Method <- "SCALE"

runtimes <- rbind(runtime_cistopic_cg, runtime_cistopic_warp, runtime_lsi, runtime_snap, runtime_scale)
runtimes <- runtimes[runtimes$Method %in% methods_keep, ]
runtimes$Method <- factor(runtimes$Method, levels = methods_keep)

###### Silhouette ######

# # use cell types defined by RNA (executed above)
clustering.use <- "celltype"
clusters <- pbmc[[clustering.use]][[1]]

get_silhouette <- function(embeddings.list, dims, clusters, method, ds) {
  df <- data.frame()
  for (i in seq_along(along.with = embeddings.list)) {
    dist.matrix <- dist(x = embeddings.list[[i]][, dims])
    sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
    res <- data.frame(
      "celltype" = clusters,
      "silhouette" = sil[, 3],
      "method" = method,
      "downsample" = ds[[i]]
    )
    df <- rbind(df, res)
  }
  return(df)
}

lsi1_sil <- get_silhouette(lsi_1, 2:20, clusters, "LSI (Signac)", pbmc_ds)
lsi2_sil <- get_silhouette(lsi_2, 1:20, clusters, "LSI (Cusanovich2018)", pbmc_ds)
lsi3_sil <- get_silhouette(lsi_3, 1:20, clusters, "LSI (log-TF)", pbmc_ds)
lsi4_sil <- get_silhouette(lsi_4, 1:20, clusters, "LSI (Cellranger)", pbmc_ds)
snap_sil <- get_silhouette(snap, c(1, 3:20), clusters, "SnapATAC", pbmc_ds)
cgs_sil <- get_silhouette(ct_cgs, 1:20, clusters, "cisTopic CGS", pbmc_ds)
warp_sil <- get_silhouette(ct_warp, 1:20, clusters, "cisTopic Warp", pbmc_ds)
scale_sil <- get_silhouette(scale_pbmc, 1:10, clusters, "SCALE", pbmc_ds)

sil_pbmc <- rbind(lsi1_sil, lsi2_sil, lsi3_sil, lsi4_sil,
                  snap_sil, cgs_sil, warp_sil, scale_sil)

sil_pbmc$downsample <- factor(sil_pbmc$downsample)

sil_pbmc_plot <- sil_pbmc[sil_pbmc$method %in% methods_keep, ]
sil_pbmc_plot$method <- factor(sil_pbmc_plot$method, levels =  methods_keep)

sil_pbmc_plot<- sil_pbmc_plot %>% 
  group_by(celltype, method, downsample) %>% 
  mutate(mn = mean(silhouette)) %>% 
  ungroup()

sil_pbmc_plot <- sil_pbmc_plot[, c("celltype", "method", "downsample", "mn")]
sil_pbmc_plot <- unique(sil_pbmc_plot)

sil_plot <- ggplot(sil_pbmc_plot, aes(x = downsample, y = mn, fill = method)) +
  geom_boxplot(outlier.size = 0.1) +
  scale_fill_manual(values = colors.use) +
  xlab("Fraction of counts retained") +
  ylab("Mean Silhouette") +
  theme_bw()

###### Figure ######

create_plot <- function(dr, object) {
  object[['dr']] <- dr
  p <- DimPlot(object, reduction = "dr", group.by = "celltype", pt.size = 0.1) +
    ggtitle("") + ylab("") + xlab("") +
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 2)))
}

# umaps
umap.use <- c(1, 5)
lsi_plot <- lapply(umaps_lsi_1[umap.use], create_plot, object = pbmc)
lsi2_plot <- lapply(umaps_lsi_2[umap.use], create_plot, object = pbmc)
scale_plot <- lapply(umaps_scale[umap.use], create_plot, object = pbmc)
snap_plot <- lapply(umaps_snap[umap.use], create_plot, object = pbmc)

lsi_plot[[1]] <- lsi_plot[[1]] + ylab("Full dataset") + ggtitle("LSI (Signac)")
lsi_plot[[2]] <- lsi_plot[[2]] + ylab("20% counts")
scale_plot[[1]] <- scale_plot[[1]] + ggtitle("SCALE")
snap_plot[[1]] <- snap_plot[[1]] + ggtitle("SnapATAC")

umaps <- wrap_plots(
  list(lsi_plot[[1]], scale_plot[[1]], snap_plot[[1]],
    lsi_plot[[2]], scale_plot[[2]], snap_plot[[2]]),
  ncol = 3,
  guides = "collect"
)

bs <- 16

p3 <- ggplot(runtimes[runtimes$Downsample == 1, ], aes(y = Seconds/60, x = Method, fill = Method)) +
  geom_bar(stat = "identity") +
  scale_y_log10() +
  ggtitle("Total run time") +
  theme_bw(base_size = bs) +
  scale_fill_manual(values = colors.use) +
  ylab("Time (minutes)") +
  xlab("") +
  theme(legend.position = 'none', axis.text.x = element_text(size = 8, angle = 25, vjust = 1, hjust=1))

sil_plot <- sil_plot + theme(legend.position = "none") + theme_bw(base_size = bs)
umaps <- umaps & theme_bw(base_size = bs)
p2 <- p2 + theme_bw(base_size = bs)

fig <- (umaps | p3) + plot_layout(widths = c(3, 1))
metrics <- (p2 / sil_plot) + plot_layout(guides = "collect")
pbmc_fig <- (fig / metrics) & theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
pbmc_fig + ggsave(filename = "figures/dimreduc_pbmc.png", height = 12, width = 16)

######## Chen #########

chen_levels <- c(250, 500, 1000, 2500, 5000)

read_lsi <- function(method, path ="data/chen/embeddings/", levels = chen_levels) {
  methodstr <- paste0(as.character(method), ".rds")
  lsi <- paste0(path, "lsi_", chen_levels, "_", methodstr)
  lsi_obj <- lapply(X = lsi, readRDS)
  return(lsi_obj)
}

lsi_chen_1 <- read_lsi(method = 1)
lsi_chen_2 <- read_lsi(method = 2)
lsi_chen_3 <- read_lsi(method = 3)
lsi_chen_4 <- read_lsi(method = 4)

snap_chen <- lapply(X = paste0("data/chen/embeddings/snapatac_", chen_levels, ".rds"), FUN = readRDS)
ct_cgs_chen <- lapply(X = paste0("data/chen/embeddings/cistopic_cgs_", chen_levels, ".rds"), FUN = readRDS)
ct_warp_chen <- lapply(X = paste0("data/chen/embeddings/cistopic_warp_", chen_levels, ".rds"), FUN = readRDS)

ct_cgs_chen <- lapply(ct_cgs_chen, t)
ct_warp_chen <- lapply(ct_warp_chen, t)

# convert h5ad to h5seurat
chen_scale_paths <- lapply(
  X = chen_levels,
  function(x) {
    Convert(
      source = paste0("data/chen/embeddings/scale_", x, "/adata.h5ad"),
      dest = paste0("data/chen/embeddings/scale_", x, "/adata.h5seurat"),
      overwrite = TRUE
    )
  }
)

# load h5seurat
scale_chen_obj <- lapply(X = chen_scale_paths, FUN = LoadH5Seurat)

# extract embeddings
scale_chen <- lapply(X = scale_chen_obj, FUN = function(x) {
  Embeddings(x[["latent"]])
})

counts <- readRDS("data/chen/scATAC-benchmarking-master/Synthetic_Data/BoneMarrow_cov5000/input/bonemarrow_cov5000.rds")
chen_obj <- CreateSeuratObject(counts = counts)
chen_obj$celltype <- chen_obj$orig.ident

######## Determine dimensions to use ##########

seqdepth_chen <- chen_obj$nCount_RNA

lsi1_depth <- lapply(lsi_chen_1, function(x) {
  which(abs(cor(x, seqdepth_chen)) > 0.9)
})
lsi2_depth <- lapply(lsi_chen_2, function(x) {
  which(abs(cor(x, seqdepth_chen)) > 0.9)
})
lsi3_depth <- lapply(lsi_chen_3, function(x) {
  which(abs(cor(x, seqdepth_chen)) > 0.9)
})
lsi4_depth <- lapply(lsi_chen_4, function(x) {
  which(abs(cor(x, seqdepth_chen)) > 0.9)
})
cgs_depth <- lapply(ct_cgs_chen, function(x) {
  which(abs(cor(x, seqdepth_chen)) > 0.9)
})
warp_depth <- lapply(ct_warp_chen, function(x) {
  which(abs(cor(x, seqdepth_chen)) > 0.9)
})
snap_depth <- lapply(snap_chen, function(x) {
  which(abs(cor(x, seqdepth_chen)) > 0.9)
})
scale_depth <- lapply(scale_chen, function(x) {
  which(abs(cor(x, seqdepth_chen)) > 0.9)
})

##### neighbors #####

get_knn_df <- function(emb_list, dims, clusters, rna_nn, method, ds_list) {
  # compute KNN purity for each dimension reduction
  knn_df <- data.frame()
  for (i in seq_along(along.with = emb_list)) {
    knn <- knn_purity(embeddings = emb_list[[i]], dims = dims, clusters = clusters, rna.nn = rna_nn)
    ds <- ds_list[[i]]
    kd <- data.frame(purity = knn, downsample = ds, method = method, celltype = clusters)
    knn_df <- rbind(knn_df, kd)
  }
  return(knn_df)
}

clusters <- unlist(lapply(strsplit(x = rownames(lsi_chen_1[[1]]), split = "_"), FUN = `[[`, 1))

knn_lsi1_chen <- get_knn_df(
  emb_list = lsi_chen_1,
  dims = 1:5,
  ds_list = chen_levels,
  clusters = clusters,
  rna_nn = NULL,
  method = "LSI (Signac)"
)

knn_lsi2_chen <- get_knn_df(
  emb_list = lsi_chen_2,
  dims = 1:5,
  ds_list = chen_levels,
  clusters = clusters,
  rna_nn = NULL,
  method = "LSI (Cusanovich2018)"
)

knn_lsi3_chen <- get_knn_df(
  emb_list = lsi_chen_3,
  dims = 1:5,
  ds_list = chen_levels,
  clusters = clusters,
  rna_nn = NULL,
  method = "LSI (log-TF)"
)

knn_lsi4_chen <- get_knn_df(
  emb_list = lsi_chen_4,
  dims = 1:5,
  ds_list = chen_levels,
  clusters = clusters,
  rna_nn = NULL,
  method = "LSI (Cellranger)"
)

knn_ct_cgs_chen <- get_knn_df(
  emb_list = ct_cgs_chen,
  dims = 1:5,
  ds_list = chen_levels,
  clusters = clusters,
  rna_nn = NULL,
  method = "cisTopic CGS"
)

knn_ct_warp_chen <- get_knn_df(
  emb_list = ct_warp_chen,
  dims = 1:5,
  ds_list = chen_levels,
  clusters = clusters,
  rna_nn = NULL,
  method = "cisTopic Warp"
)

knn_snap_chen <- get_knn_df(
  emb_list = snap_chen,
  dims = 1:5,
  ds_list = chen_levels,
  clusters = clusters,
  rna_nn = NULL,
  method = "SnapATAC"
)

knn_scale_chen <- get_knn_df(
  emb_list = scale_chen,
  dims = 1:10,
  ds_list = chen_levels,
  clusters = clusters,
  rna_nn = NULL,
  method = "SCALE"
)

knn_df_chen <- rbind(knn_lsi1_chen, knn_lsi2_chen, knn_lsi3_chen, knn_lsi4_chen,
                     knn_ct_cgs_chen, knn_ct_warp_chen, knn_snap_chen, knn_scale_chen)

knn_df_chen_plot <- knn_df_chen[knn_df_chen$method %in% methods_keep, ]
knn_df_chen_plot$method <- factor(knn_df_chen_plot$method, levels =  methods_keep)

knn_df_chen_plot <- knn_df_chen_plot %>% 
  group_by(downsample, method, celltype) %>% 
  mutate(mn = mean(purity)) %>% 
  ungroup()

knn_df_chen_plot <- knn_df_chen_plot[, c("mn", "method", "downsample")]
knn_df_chen_plot <- unique(knn_df_chen_plot)

# compute UMAP for each
umaps_lsi_1_chen <- lapply(lsi_chen_1, function(x) RunUMAP(x[, 1:5]))
umaps_lsi_2_chen <- lapply(lsi_chen_2, function(x) RunUMAP(x[, 1:5]))
umaps_lsi_3_chen <- lapply(lsi_chen_3, function(x) RunUMAP(x[, 1:5]))
umaps_lsi_4_chen <- lapply(lsi_chen_4, function(x) RunUMAP(x[, 1:5]))
umaps_snap_chen <- lapply(snap_chen, function(x) {
  rownames(x) <- colnames(chen_obj)
  dr <- RunUMAP(x[, 1:5])
  dr
}
)
umaps_ct_cgs_chen<- lapply(ct_cgs_chen, function(x) {
  rownames(x) <- colnames(chen_obj)
  dr <- RunUMAP(x[, 1:5])
  dr
}
)
umaps_ct_warp_chen <- lapply(ct_warp_chen, function(x) {
  rownames(x) <- colnames(chen_obj)
  dr <- RunUMAP(x[, 1:5])
  dr
}
)
umaps_scale_chen <- lapply(scale_chen, function(x) RunUMAP(x[, 1:10]))

## Silhouette ##

lsi1_sil_chen <- get_silhouette(lsi_chen_1, 1:5, clusters, "LSI (Signac)", chen_levels)
lsi2_sil_chen <- get_silhouette(lsi_chen_2, 1:5, clusters, "LSI (Cusanovich2018)", chen_levels)
lsi3_sil_chen <- get_silhouette(lsi_chen_3, 1:5, clusters, "LSI (log-TF)", chen_levels)
lsi4_sil_chen <- get_silhouette(lsi_chen_4, 1:5, clusters, "LSI (Cellranger)", chen_levels)
snap_sil_chen <- get_silhouette(snap_chen, 1:5, clusters, "SnapATAC", chen_levels)
cgs_sil_chen <- get_silhouette(ct_cgs_chen, 1:5, clusters, "cisTopic CGS", chen_levels)
warp_sil_chen <- get_silhouette(ct_warp_chen, 1:5, clusters, "cisTopic Warp", chen_levels)
scale_sil_chen <- get_silhouette(scale_chen, 1:10, clusters, "SCALE", chen_levels)

sil_chen <- rbind(lsi1_sil_chen, lsi2_sil_chen, lsi3_sil_chen, lsi4_sil_chen,
                  snap_sil_chen, cgs_sil_chen, warp_sil_chen, scale_sil_chen)

sil_chen$downsample <- factor(sil_chen$downsample)

sil_chen_plot <- sil_chen[sil_chen$method %in% methods_keep, ]
sil_chen_plot$method <- factor(sil_chen_plot$method, levels =  methods_keep)

sil_chen_plot <- sil_chen_plot %>% 
  group_by(downsample, method, celltype) %>% 
  mutate(mn = mean(silhouette)) %>% 
  ungroup()
sil_chen_plot <- sil_chen_plot[, c("downsample", "method", "mn")]
sil_chen_plot <- unique(sil_chen_plot)

sil_chen <- ggplot(sil_chen_plot, aes(x = downsample, y = mn, fill = method)) +
  geom_boxplot(outlier.size = 0.1) +
  scale_fill_manual(values = colors.use) +
  xlab("Average counts per cell") +
  ylab("Mean Silhouette") +
  theme_bw()

######### Runtimes ##########

runtime_lsi <- read.table("data/chen/embeddings/lsi_runtime.txt", sep = "\t")
colnames(runtime_lsi) <- c("Seconds", "Downsample", "Method")
runtime_lsi[runtime_lsi$Method == 1, "Method"] <- "LSI (Signac)"
runtime_lsi[runtime_lsi$Method == 2, "Method"] <- "LSI (Cusanovich2018)"
runtime_lsi[runtime_lsi$Method == 3, "Method"] <- "LSI (log-TF)"
runtime_lsi[runtime_lsi$Method == 4, "Method"] <- "LSI (Cellranger)"

runtime_snap <- read.table("data/chen/embeddings/snapatac_runtime.txt", sep = "\t")
colnames(runtime_snap) <- c("Seconds", "Downsample")
runtime_snap$Method <- "SnapATAC"

runtime_cistopic_cg <- read.table("data/chen/embeddings/cistopic_cgs_runtime.txt", sep = "\t")
colnames(runtime_cistopic_cg) <- c("Seconds", "Downsample")
runtime_cistopic_cg$Method <- "cisTopic CGS"

runtime_cistopic_warp <- read.table("data/chen/embeddings/cistopic_warp_runtime.txt", sep = "\t")
colnames(runtime_cistopic_warp) <- c("Seconds", "Downsample")
runtime_cistopic_warp$Method <- "cisTopic Warp"

runtime_scale <- read.table("data/chen/embeddings/scale_runtime.txt", sep = "\t")
colnames(runtime_scale) <- c("Seconds", "Downsample")
runtime_scale$Method <- "SCALE"

runtimes <- rbind(runtime_cistopic_cg, runtime_cistopic_warp, runtime_lsi, runtime_snap, runtime_scale)
runtimes <- runtimes[runtimes$Method %in% methods_keep, ]
runtimes$Method <- factor(runtimes$Method, levels = methods_keep)

chen_runtimes <- ggplot(runtimes[runtimes$Downsample == 5000, ], aes(y = Seconds, x = Method, fill = Method)) +
  geom_bar(stat = "identity") +
  scale_y_log10() +
  theme_bw() +
  scale_fill_manual(values = colors.use) +
  ylab("Time (seconds)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

#### Figure #####

# create figure with ARI and UMAPs for chen dataset (supplementary figure)

lsi_plot <- lapply(umaps_lsi_1_chen, create_plot, object = chen_obj)
lsi2_plot <- lapply(umaps_lsi_2_chen, create_plot, object = chen_obj)
lsi3_plot <- lapply(umaps_lsi_3_chen, create_plot, object = chen_obj)
scale_plot <- lapply(umaps_scale_chen, create_plot, object = chen_obj)
cgs_plot <- lapply(umaps_ct_cgs_chen, create_plot, object = chen_obj)
warp_plot <- lapply(umaps_ct_warp_chen, create_plot, object = chen_obj)
snap_plot <- lapply(umaps_snap_chen, create_plot, object = chen_obj)

# set axis grid names
lsi_plot[[1]] <- lsi_plot[[1]] + ylab("LSI (Signac)") + ggtitle("250")
lsi_plot[[2]] <- lsi_plot[[2]] + ggtitle("500")
lsi_plot[[3]] <- lsi_plot[[3]] + ggtitle("1000")
lsi_plot[[4]] <- lsi_plot[[4]] + ggtitle("2500")
lsi_plot[[5]] <- lsi_plot[[5]] + ggtitle("5000")

lsi2_plot[[1]] <- lsi2_plot[[1]] + ylab("LSI (Cusanovich2018)")
lsi3_plot[[1]] <- lsi3_plot[[1]] + ylab("LSI (log-TF)")
scale_plot[[1]] <- scale_plot[[1]] + ylab("SCALE")
cgs_plot[[1]] <- cgs_plot[[1]] + ylab("cisTopic CGS")
warp_plot[[1]] <- warp_plot[[1]] + ylab("cisTopic Warp")
snap_plot[[1]] <- snap_plot[[1]] + ylab("SnapATAC")

p1 <- wrap_plots(
  c(lsi_plot, lsi2_plot, lsi3_plot, scale_plot, cgs_plot, warp_plot, snap_plot),
  nrow = 7,
  guides = "collect"
)

p1 <- p1 & theme(plot.margin =  unit(c(1,1,1,1), "mm"))

p2 <- ggplot(knn_df_chen_plot, aes(x = as.factor(downsample), y = mn, fill = method)) +
  geom_boxplot(outlier.size = 0.1) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_fill_manual(values = colors.use) +
  ylab("Mean kNN celltype purity") +
  xlab("Average counts per cell")

metrics <- (p2 / sil_chen) + plot_layout(guides = "collect")

pp <- (p1 / metrics) + plot_layout(heights = c(3, 1)) & theme(axis.text = element_text(size=8))
ggsave(filename = "figures/dimreduc_chen.png", plot = pp, width = 10, height = 16, units = "in", dpi = 400)
