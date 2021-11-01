library(Signac)
library(Seurat)

pbmc <- readRDS("objects/pbmc.rds")

# label transfer
# create separate object from the RNA assay
pbmc.rna <- CreateSeuratObject(
  counts = GetAssayData(pbmc, assay = "RNA", slot = "counts"),
  meta.data = pbmc[[]]
)
pbmc.rna <- NormalizeData(pbmc.rna)
pbmc.rna <- FindVariableFeatures(pbmc.rna, nfeatures = 3000)
pbmc.rna <- ScaleData(pbmc.rna)

# Identify anchors
transfer.anchors <- FindTransferAnchors(
  reference = pbmc.rna,
  query = pbmc,
  features = VariableFeatures(object = pbmc.rna), 
  reference.assay = "RNA",
  query.assay = "GA",
  reduction = "cca",
  dims = 1:30
)

pbmc.rna$ct <- pbmc$celltype
pbmc.rna$coarse_celltype <- pbmc$coarse_celltype

celltype.predictions <- TransferData(
  anchorset = transfer.anchors,
  refdata = pbmc.rna$ct,
  weight.reduction = pbmc[["lsi"]],
  dims = 2:30
)

coarse.predictions <- TransferData(
  anchorset = transfer.anchors,
  refdata = pbmc.rna$coarse_celltype,
  weight.reduction = pbmc[["lsi"]],
  dims = 2:30
)

pbmc.atac <- pbmc
pbmc.atac <- AddMetaData(pbmc.atac, metadata = celltype.predictions)
pbmc.atac$coarse_predicted <- coarse.predictions$predicted.id

# remove unneeded assays
pbmc.atac[["RNA"]] <- NULL
pbmc.atac[["SCT"]] <- NULL
pbmc.atac[["GA"]] <- NULL

pbmc.atac$gt <- as.factor(pbmc.rna$ct)
saveRDS(object = pbmc.atac, file = "objects/multimodal_label_transfer.rds")
