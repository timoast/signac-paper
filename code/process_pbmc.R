library(Signac)
library(Seurat)
library(SeuratDisk)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)

set.seed(1234)

counts <- Read10X_h5("data/pbmc/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
fragments <- "data/pbmc/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"

annotations <- GetGRangesFromEnsDb(EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"

# create object
pbmc <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  project = "coassay",
  assay = "RNA"
)

# create initial assay using cellranger counts
counts <- counts$Peaks
grange.counts <- StringToGRanges(rownames(counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
counts <- counts[as.vector(grange.use), ]

pbmc[["cellranger"]] <- CreateChromatinAssay(
  counts = counts,
  genome = "hg38",
  sep = c(":", "-"),
  fragments = fragments,
  annotation = annotations
)

DefaultAssay(pbmc) <- "cellranger"

# QC
pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc, fast = FALSE)

pbmc <- subset(
  x = pbmc,
  subset = nCount_cellranger < 100000 &
    nCount_RNA < 25000 &
    nCount_cellranger > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)

# Process RNA and annotate cell types
DefaultAssay(pbmc) <- "RNA"

pbmc <- NormalizeData(pbmc)
pbmc <- SCTransform(pbmc)
pbmc <- RunPCA(pbmc)
pbmc <- RunUMAP(pbmc, dims = 1:50, reduction.name = "umap.rna")
pbmc <- FindNeighbors(pbmc, dims = 1:50)
pbmc <- FindClusters(pbmc, algorithm = 3, resolution = 1.5)

# map to PBMC reference
reference <- LoadH5Seurat("data/pbmc/pbmc_multimodal.h5seurat")

transfer_anchors <- FindTransferAnchors(
  reference = reference,
  query = pbmc,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)

predictions <- TransferData(
  anchorset = transfer_anchors, 
  refdata = reference$celltype.l2,
  weight.reduction = pbmc[['pca']],
  dims = 1:50
)

pbmc <- AddMetaData(
  object = pbmc,
  metadata = predictions
)

# for spurious cell annotations, assign to most abundant classification in cells cluster
ct.remove <- names(which(table(pbmc$predicted.id) < 10))

# reassign erythrocytes since we know these are nuclei
ct.remove <- c(ct.remove, "Eryth")

cell.reassign <- WhichCells(pbmc, expression = predicted.id %in% ct.remove)
pbmc$celltype <- pbmc$predicted.id
nn.graph <- pbmc@graphs$SCT_nn
for (i in cell.reassign) {
  # most frequent predicted.id of neighboring cells
  nn.cells <- names(which(nn.graph[i, ] > 0))
  nn.cells <- setdiff(nn.cells, cell.reassign)
  celltype.reassign <- names(sort(table(pbmc$predicted.id[nn.cells]), decreasing = TRUE))[[1]]
  pbmc$celltype[i] <- celltype.reassign
}

# rename cDC2 to cDC
pbmc$celltype <- ifelse(pbmc$celltype == "cDC2", "cDC", pbmc$celltype)

# call peaks for each cell type using MACS2
DefaultAssay(pbmc) <- "cellranger"
peaks <- CallPeaks(
  object = pbmc,
  group.by = "celltype",
  additional.args = "--max-gap 50"
)

peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")

# remove peaks in blacklist regions
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify peaks
peakcounts <- FeatureMatrix(
  fragments = Fragments(pbmc),
  features = peaks,
  cells = colnames(pbmc)
)

pbmc[["ATAC"]] <- CreateChromatinAssay(
  counts = peakcounts,
  min.cells = 5,
  genome = "hg38",
  fragments = fragments,
  annotation = annotations
)

DefaultAssay(pbmc) <- "ATAC"

pbmc <- FindTopFeatures(pbmc, min.cutoff = 10)
pbmc <- RunTFIDF(pbmc)
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(pbmc, reduction = "lsi", dims = 2:40, reduction.name = "umap.atac")
pbmc <- FindNeighbors(pbmc, reduction = "lsi", dims = 2:40)
pbmc <- FindClusters(pbmc, algorithm = 3)

# # subcluster CD8 Naive
cd8_stress <- pbmc[, pbmc$celltype == "CD8 Naive"]
cd8_stress <- WhichCells(cd8_stress, expression = seurat_clusters == 11)
pbmc$celltype <- ifelse(colnames(pbmc) %in% cd8_stress, "CD8 Naive JUN-", pbmc$celltype)

# joint clustering
pbmc <- FindMultiModalNeighbors(
  object = pbmc,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

pbmc <- RunUMAP(
  object = pbmc,
  nn.name = "wknn",
  reduction.name = "jumap",
  reduction.key = "JUMAP_",
  assay = "RNA",
  verbose = TRUE
)

Idents(pbmc) <- "celltype"

# add motif information
library(TFBSTools)
library(JASPAR2020)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BiocParallel)
register(SerialParam())

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

pbmc <- AddMotifs(
  object = pbmc,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm,
  assay = "ATAC",
  verbose = TRUE
)

pbmc <- RunChromVAR(
  object = pbmc,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# do the same for cellranger assay
pbmc <- AddMotifs(
  object = pbmc,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm,
  assay = "cellranger",
  verbose = TRUE
)

# add gene activities
DefaultAssay(pbmc) <- "ATAC"

ga <- GeneActivity(object = pbmc)
pbmc[["GA"]] <- CreateAssayObject(counts = ga)

levels(pbmc) <- c("CD4 Naive", "CD4 TCM", "CD4 CTL", "CD4 TEM", "CD8 Naive",
                 "CD8 Naive JUN-", "CD8 TEM", "CD8 TCM", "MAIT", "NK", "NK_CD56bright",
                 "NK Proliferating", "gdT",
                 "Treg", "B naive", "B intermediate", "B memory", "Plasmablast",
                 "CD14 Mono", "CD16 Mono",
                 "cDC", "pDC", "HSPC")

saveRDS(object = pbmc, file = "objects/pbmc.rds")
