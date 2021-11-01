library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(BSgenome.Hsapiens.UCSC.hg38)

pbmc <- readRDS("objects/pbmc.rds")
lnk <- readRDS("objects/pbmc_links.rds")

DefaultAssay(pbmc) <- "ATAC"
Links(pbmc) <- lnk

# ----- QC plots -----

nucleosome_plot <- FragmentHistogram(pbmc, group.by = "orig.ident") +
  ggtitle("Nucleosome signal")

tss_plot <- TSSPlot(pbmc, assay = "cellranger", group.by = "orig.ident")

saveRDS(object = nucleosome_plot, file = "figures/nucleosome_signal.rds")
saveRDS(object = tss_plot, file = "figures/tss_enrichment.rds")

# ----- Dim plots ----- 

rna_dimplot <- DimPlot(pbmc, reduction = "umap.rna", label = TRUE, repel = TRUE)
atac_dimplot <- DimPlot(pbmc, reduction = "umap.atac", label = TRUE, repel = TRUE)

saveRDS(object = rna_dimplot, file = "figures/pbmc_rna_dimplot.rds")
saveRDS(object = atac_dimplot, file = "figures/pbmc_atac_dimplot.rds")

# ----- Markers -----

markers <- FindMarkers(
  object = pbmc,
  ident.1 = "CD8 TEM",
  ident.2 = "CD8 Naive",
  test.use = "LR",
  latent.vars = "nCount_ATAC",
  only.pos = TRUE
)

top.markers <- markers[markers$p_val_adj < 0.01 & markers$avg_log2FC > 0.4, ]

motifs <- FindMotifs(
  object = pbmc,
  features = rownames(top.markers),
  features.match = c("GC.percent", "count", "sequence.length")
)

# EOMES, TBX21, TBX2 equally enriched in effector T cell peaks
# look at RNA data to see which is expressed
# compare with chromvar deviations

DefaultAssay(pbmc) <- "RNA"
tf_use <- c("EOMES", "TBX21", "TBX2")

tf_expression <- VlnPlot(
  object = pbmc,
  features = tf_use,
  idents = c("CD8 TEM", "CD8 Naive"),
  pt.size = 0
) & ylim(c(0, 2.5)) & ylab("RNA expression") & ggtitle("") & xlab("")

DefaultAssay(pbmc) <- "chromvar"

tf_chromvar <- lapply(X = seq_along(tf_use), function(x) {
  VlnPlot(
    object = pbmc,
    features = motifs$motif[x],
    idents = c("CD8 TEM", "CD8 Naive"),
    pt.size = 0
  ) + ggtitle(tf_use[x]) + NoLegend() + ylab("chromVAR deviation") +
    xlab("") + theme(axis.text.x = element_blank())
})

tf_chromvar <- wrap_plots(tf_chromvar, ncol = 3)

DefaultAssay(pbmc) <- "ATAC"

pbmc <- Footprint(
  object = pbmc,
  motif.name = c("EOMES", "TBX21"),
  genome = BSgenome.Hsapiens.UCSC.hg38
)

fp <- PlotFootprint(
  object = pbmc,
  features = c("EOMES", "TBX21"),
  idents = c("CD8 TEM", "CD8 Naive")
) & NoLegend() & plot_layout(ncol = 1)

mp <- MotifPlot(pbmc, head(motifs$motif, 3))

saveRDS(object = mp, file = "figures/motifplot.rds")
saveRDS(object = tf_chromvar, file = "figures/chromvar_vln.rds")
saveRDS(object = tf_expression, file = "figures/tf_rna_vln.rds")
saveRDS(object = fp, file = "figures/footprint.rds")

# ----- Coverage plot -----

covplot <- CoveragePlot(
  object = pbmc,
  idents = c("CD4 Naive", "CD4 TCM", "CD8 Naive",
             "CD8 TEM", "MAIT", "NK", "Treg"),
  region = "CD8A",
  features = "CD8A",
  expression.assay = "RNA",
  extend.upstream = 2000,
  extend.downstream = 2000,
  links = FALSE
)

saveRDS(object = covplot, file = "figures/coverage_plot.rds")

# ----- Link analysis -----

# ratio of positive to negative links
sum(lnk$score < 0) / length(lnk) * 100

# total over 100 kb
sum(width(lnk) > 100000) / length(lnk)

# number of links per gene (regulatory complexity)
# compare cell-type-specific vs houskeeping genes

# for each gene, find the number of linked peaks
link.df <- as.data.frame(lnk)

links_per_gene <- link.df %>% 
  mutate(pos_link = score > 0) %>% 
  group_by(gene) %>% 
  summarise(positive_links = sum(pos_link), negative_links = sum(!pos_link))

mean(links_per_gene$positive_links + links_per_gene$negative_links)
# 6.373724

sd(links_per_gene$positive_links + links_per_gene$negative_links)
# 7.110643

# total links per gene
link_per_gene_plot <- links_per_gene %>%
  group_by(positive_links, negative_links) %>%
  summarise(count = n()) %>% 
  ggplot(data = ., aes(x = positive_links, y = negative_links, fill = log10(count+1))) +
  geom_tile() +
  theme_bw() +
  scale_fill_viridis_c() +
  ylab("Total negative links") +
  xlab("Total positive links") +
  ggtitle("Number of linked peaks per gene")

# number of linked genes per peak
genes_per_link <- link.df %>% 
  mutate(pos_link = score > 0) %>% 
  group_by(peak) %>% 
  summarise(positive_links = sum(pos_link), negative_links = sum(!pos_link))

mean(genes_per_link$positive_links + genes_per_link$negative_links)
# 1.578854

sd(genes_per_link$positive_links + genes_per_link$negative_links)
# 1.259847

# total links per gene
gene_per_link_plot <- genes_per_link %>%
  group_by(positive_links, negative_links) %>%
  summarise(count = n()) %>% 
  ggplot(data = ., aes(x = positive_links, y = negative_links, fill = log10(count+1))) +
  geom_tile() +
  theme_bw() +
  scale_x_continuous(breaks = 0:10) +
  scale_y_continuous(breaks = 0:10) +
  scale_fill_viridis_c() +
  ylab("Total negative links") +
  xlab("Total positive links") +
  ggtitle("Number of linked genes per peak")

# distance from peak to tss
p1 <- ggplot(data = link.df[link.df$score > 0, ], aes(x = width)) +
  geom_histogram(bins = 100) +
  theme_classic() +
  xlab("") +
  ylab("Count") +
  ggtitle("Positive gene associations")

p2 <- ggplot(data = link.df[link.df$score < 0, ], aes(x = width)) +
  geom_histogram(bins = 100) +
  theme_classic() +
  xlab("Distance to gene TSS (bp)") +
  ylab("Count") +
  ggtitle("Negative gene associations")

p3 <- ggplot(data = link.df, mapping = aes(x = pvalue)) +
  geom_histogram(bins = 100) +
  theme_classic() +
  xlab("p-value") +
  ylab("Count") +
  ggtitle("p-value distribution")

saveRDS(object = p1, file = "figures/distance_positive.rds")
saveRDS(object = p2, file = "figures/distance_negative.rds")
saveRDS(object = p3, file = "figures/link_pvals.rds")
saveRDS(object = gene_per_link_plot, file = "figures/genes_per_link_plot.rds")
saveRDS(object = link_per_gene_plot, file = "figures/link_per_gene_plot.rds")

# ----- Link plots ----- 

linked_1 <- CoveragePlot(
  object = pbmc,
  region = "MS4A1",
  features = "MS4A1",
  idents = c("B naive", "B intermediate", "B memory", "CD14 Mono", "CD16 Mono", "CD8 TEM", "CD8 Naive"),
  extend.upstream = 500,
  extend.downstream = 10000
)

linked_2 <- CoveragePlot(
  object = pbmc,
  region = "LYZ",
  features = "LYZ",
  idents = c("B naive", "B intermediate", "B memory", "CD14 Mono", "CD16 Mono", "CD8 TEM", "CD8 Naive"),
  extend.upstream = 5000,
  extend.downstream = 5000
)

saveRDS(object = linked_1, file = "figures/linked_covplot1.rds")
saveRDS(object = linked_2, file = "figures/linked_covplot2.rds")

# ----- Peak calling comparison ----- 

DefaultAssay(pbmc) <- "cellranger"

# call MACS2 peaks on pseudobulk
pks <- CallPeaks(
  object = pbmc,
  macs2.path = "/home/stuartt/miniconda3/envs/signac/bin/macs2",
  additional.args = "--max-gap 50"
)

# example where cellranger incorrectly merges peaks
cp <- CoveragePlot(
  object = pbmc,
  region = "CD8A",
  ranges = pks,
  ranges.title = "MACS2",
  extend.upstream = 2000,
  extend.downstream = 2000,
  links = FALSE
)
saveRDS(object = cp, file = "figures/cellranger_peakcalling.rds")

DefaultAssay(pbmc) <- "ATAC"

# example where celltype specific peaks missed (need to also run MACS2 on bulk and compare)
# find markers for a rare population

mrk_cd56 <- FindMarkers(
  object = pbmc,
  ident.1 = "NK_CD56bright",
  ident.2 = "NK",
  latent.vars = "nCount_ATAC",
  test.use = "LR",
  only.pos = TRUE
)

all.markers <- FindAllMarkers(
  object = pbmc,
  test.use = "LR",
  latent.vars = "nCount_ATAC",
  only.pos = TRUE
)

all.markers$isunique <- Biobase::isUnique(all.markers$gene)
all.unique.markers <- all.markers[all.markers$isunique, ]

n_celltype <- table(Idents(pbmc))
fraction_recovered <- vector(mode = 'numeric', length = length(n_celltype))
for (i in seq_along(n_celltype)) {
  celltype <- names(n_celltype)[[i]]
  markers.use <- all.unique.markers[all.unique.markers$cluster == celltype, ]
  markers.ranges <- StringToGRanges(markers.use$gene)
  frac_recovered <- sum(countOverlaps(query = markers.ranges, subject = pks) > 0) / length(markers.ranges)
  fraction_recovered[[i]] <- frac_recovered
}

df <- data.frame(n_cells = n_celltype, fraction_recovered = fraction_recovered)

missed_peak_count <- ggplot(df, aes(n_celltype, fraction_recovered)) +
  geom_point() +
  theme_classic() +
  ylab("Fraction of cell-type-specific peaks identified") +
  xlab("Number of cells")

missed <- CoveragePlot(
  object = pbmc,
  region = "chr19-3805000-3806000",
  ranges.title = "Bulk",
  ranges = pks,
  extend.upstream = 5000,
  extend.downstream = 5000,
  links = FALSE
)

saveRDS(object = missed_peak_count, file = "figures/missed_peak_count.rds")
saveRDS(object = missed, file = "figures/macs2_pseudobulk.rds")

# get average number of macs2 peaks overlapped by cellranger peak
olap <- findOverlaps(query = pks, subject = pbmc[["cellranger"]])

sum(table(subjectHits(olap)) > 1) # 13751
sum(table(queryHits(olap)) > 1)   # 2

# call peaks per DNA accessibility cluster
cluster_peaks <- CallPeaks(
  object = pbmc,
  group.by = "cellranger_snn_res.0.8",
  additional.args = "--max-gap 50",
  macs2.path = "/home/stuartt/miniconda3/envs/signac/bin/macs2"
)

# check overlap with cell-type specific peaks
olap <- findOverlaps(query = pbmc[["ATAC"]], subject = cluster_peaks)
sum(table(subjectHits(olap)) >= 1)  / length(granges(pbmc[["ATAC"]])) # 0.9270874

sum(table(subjectHits(olap)) > 1) # 1059
sum(table(queryHits(olap)) > 1)   # 2497
length(cluster_peaks) # 152473
nrow(pbmc[["ATAC"]])  # 155611

olap <- findOverlaps(query = pbmc[["ATAC"]], subject = pks)
sum(table(subjectHits(olap)) >= 1)  / length(granges(pbmc[["ATAC"]])) # 0.7859085

sum(table(subjectHits(olap)) > 1) # 3529
sum(table(queryHits(olap)) > 1)   # 1310
length(pks) # 134195
nrow(pbmc[["ATAC"]])  # 155611

# DimPlot showing ATAC clusters
cluster_dimplot <- DimPlot(pbmc, group.by = "cellranger_snn_res.0.8", label = TRUE, repel = TRUE, reduction = "umap.atac") +
  ggtitle("scATAC-seq cell clusters") + theme(legend.position = "none") +
  xlab("UMAP 1") + ylab("UMAP 2")

ggsave(filename = "figures/atac_cluster_dimplot.png", plot = cluster_dimplot, height = 8, width = 8, dpi = 400)
saveRDS(cluster_dimplot, "figures/atac_cluster_dimplot.rds")

## ------- Multimodal label transfer ----------

pbmc.atac <- readRDS("objects/multimodal_label_transfer.rds")

# high-res prediction
pbmc.atac$predicted.id <- factor(pbmc.atac$predicted.id, levels = levels(pbmc.atac$gt))
pbmc.atac$annotation_correct <- pbmc.atac$predicted.id == pbmc.atac$gt
p1 <- DimPlot(pbmc.atac, group.by = "gt", label = TRUE, repel = TRUE, reduction = "umap.atac") + NoLegend() + ggtitle("Ground-truth annotation")
p2 <- DimPlot(pbmc.atac, group.by = "predicted.id", label = TRUE, repel = TRUE, reduction = "umap.atac") + NoLegend() + ggtitle("Predicted annotation")

predictions <- table(pbmc.atac$gt, pbmc.atac$predicted.id)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions <- as.data.frame(predictions)
p3 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  scale_fill_viridis_c() +
  xlab("Annotated cell type (RNA)") +
  ylab("Predicted cell type (ATAC)") + 
  labs(fill = "Fraction of cells") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

incorrect <- length(which(pbmc.atac$gt != pbmc.atac$predicted.id))
1 - (incorrect / ncol(pbmc.atac)) # 0.870151

data <- FetchData(pbmc.atac, vars = c("prediction.score.max", "annotation_correct"))
p4 <- ggplot(data, aes(prediction.score.max, fill = annotation_correct, colour = annotation_correct)) + 
  geom_histogram() +
  facet_wrap(~annotation_correct) +
  xlab("Prediction Score") +
  theme_bw()

ggsave(filename = "figures/label_transfer_accuracy.pdf", plot = p3, height = 4.5, width = 6)
saveRDS(object = p3, file = "figures/label_transfer_accuracy.rds")
(p1 | p2 | p4) + ggsave("figures/multimodal_label_transfer.png", height = 8, width = 18)

# coarse prediction
incorrect <- length(which(pbmc.atac$coarse_celltype != pbmc.atac$coarse_predicted))
1 - (incorrect / ncol(pbmc.atac)) # 0.925473
