library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(EnsDb.Hsapiens.v75)
library(dplyr)

# load counts and metadata from cellranger-atac
counts <- Read10X_h5(filename = "data/mito/CRC_v12-mtMask_mgatk.filtered_peak_bc_matrix.h5")

metadata <- read.csv(
  file = "data/mito/CRC_v12-mtMask_mgatk.singlecell.csv",
  header = TRUE,
  row.names = 1
)

# load gene annotations from Ensembl
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg19"

# create object
crc_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  annotation = annotations,
  min.cells = 10,
  genome = "hg19",
  fragments = 'data/mito/CRC_v12-mtMask_mgatk.fragments.tsv.gz'
)
crc <- CreateSeuratObject(
  counts = crc_assay,
  assay = 'peaks',
  meta.data = metadata
)

# Augment QC metrics that were computed by cellranger-atac
crc$pct_reads_in_peaks <- crc$peak_region_fragments / crc$passed_filters * 100
crc$pct_reads_in_DNase <- crc$DNase_sensitive_region_fragments / crc$passed_filters * 100
crc$blacklist_ratio <- crc$blacklist_region_fragments / crc$peak_region_fragments

# compute TSS enrichment score and nucleosome banding pattern
crc <- TSSEnrichment(crc)
crc <- NucleosomeSignal(crc)

# remove low-quality cells
crc <- subset(
  x = crc,
  subset = nCount_peaks > 1000 &
    nCount_peaks < 50000 &
    pct_reads_in_DNase > 40 &
    blacklist_ratio < 0.05 &
    TSS.enrichment > 3 &
    nucleosome_signal < 4
)

# load mGATK output
mito.data <- ReadMGATK(dir = "data/mito/")

# create an assay
mito <- CreateAssayObject(counts = mito.data$counts)

# Subset to cell present in the scATAC-seq assay
mito <- subset(mito, cells = colnames(crc))

# add assay and metadata to the seurat object
crc[["mito"]] <- mito
crc <- AddMetaData(crc, metadata = mito.data$depth, col.name = "mtDNA_depth")

# filter cells based on mitochondrial depth
crc <- subset(crc, mtDNA_depth >= 10)

crc <- RunTFIDF(crc)
crc <- FindTopFeatures(crc, min.cutoff = 10)
crc <- RunSVD(crc)
crc <- RunUMAP(crc, reduction = "lsi", dims = 2:50)
crc <- FindNeighbors(crc, reduction = "lsi", dims = 2:50)
crc <- FindClusters(crc, resolution = 0.5, algorithm = 3)

crc <- RenameIdents(
  object = crc,
  '0' = 'Epithelial',
  '1' = 'Epithelial',
  '2' = 'Basophil',
  '3' = 'Myeloid_1',
  '4' = 'Myeloid_2',
  '5' = 'T-cell'
)
crc$celltype <- Idents(crc)

dp <- DimPlot(crc, label = TRUE, group.by = "celltype", pt.size = 0.1) + NoLegend()
saveRDS(object = dp, file = "figures/mito_dimplot.rds")

variable.sites <- IdentifyVariants(crc, assay = "mito", refallele = mito.data$refallele)
varplot <- VariantPlot(variants = variable.sites)
saveRDS(object = varplot, file = "figures/mito_varplot.rds")

# Establish a filtered data frame of variants based on this processing
high.conf <- subset(
  variable.sites, subset = n_cells_conf_detected >= 5 &
    strand_correlation >= 0.65 &
    vmr > 0.01
)

crc <- AlleleFreq(
  object = crc,
  variants = high.conf$variant,
  assay = "mito"
)

DefaultAssay(crc) <- "alleles"
alleles.view <- c("16147C>T", "824T>C")

allele.plot <- FeaturePlot(
  object = crc,
  features = alleles.view,
  order = TRUE,
  cols = c("grey", "darkred"),
  ncol = 1,
  pt.size = 0.5
) + plot_layout(guides = "collect")

saveRDS(object = allele.plot, file = "figures/mito_allele_plot.rds")

crc <- FindClonotypes(object = crc)

prop <- as.data.frame(table(crc$celltype, Idents(crc)))
colnames(prop) <- c("Celltype", "Clone", "Count")

prop <- prop %>%
  group_by(Clone) %>%
  mutate(n.pred = sum(Count)) %>%
  mutate("Fraction cells" = Count / n.pred) %>% 
  group_by(Celltype) %>% 
  mutate(n.celltype = sum(Count)) %>% 
  mutate("Normalized cell\nfraction" = `Fraction cells` / n.celltype)

prop$Clone <- as.character(prop$Clone)

hm <- ggplot(as.data.frame(prop), aes(Celltype, Clone, fill = `Normalized cell\nfraction`)) +
  geom_tile() +
  theme_classic() +
  scale_fill_viridis_c(option = "B")

saveRDS(object = hm, file = "figures/mito_clone_hm.rds")

# find DA peaks between clones within the same cell type
DefaultAssay(crc) <- "peaks"

# find markers for each epithelial clone
epi_clones <- c(1, 2, 4)

# subset to epithelial cells to find differences only among epithelial clones
epithelial <- subset(x = crc, subset = celltype == "Epithelial")

all.markers <- list()
for (i in seq_along(epi_clones)) {
  all.markers[[i]] <- FindMarkers(
    object = epithelial,
    ident.1 = epi_clones[[i]],
    test.use = "LR",
    latent.vars = "nCount_peaks",
    only.pos = TRUE
  )
}

covplots <- list()
for (i in seq_along(all.markers)) {
  covplots[[i]] <- CoveragePlot(
    object = crc,
    region = rownames(all.markers[[i]])[1:2],
    ncol = 1,
    group.by = "alleles_snn_res.1",
    extend.upstream = 500,
    extend.downstream = 500,
    idents = c(1,2,4),
    peaks = FALSE,
    annotation = FALSE
  )
}

covplots <- wrap_plots(covplots, ncol = 3)

saveRDS(object = covplots, file = "figures/mito_covplot.rds")
saveRDS(object = crc, file = "objects/mito.rds")