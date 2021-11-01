library(ArchR)
library(Seurat)
library(Signac)
set.seed(1234)

# create object from each downsampled fragment file
options(scipen=999)
downsamples_biccn <- c(50000, 100000, 200000, 300000, 400000, 500000, 600000, 700000)
downsamples_pbmc <- seq(from = 1000, to = 26000, by = 2000)

addArchRThreads(threads = 1)
addArchRGenome("hg19")

for (i in downsamples_pbmc) {
  # load the signac fragment file to get list of cells to include
  frags <- readRDS(paste0("/scratch/tim/pbmc_atac/downsampling/", i, ".rds"))
  cells <- Cells(frags)
  start.time <- Sys.time()
  ArrowFiles <- createArrowFiles(
    inputFiles = paste0("/scratch/tim/pbmc_atac/downsampling/", i, ".bed.gz"),
    sampleNames = paste0("pbmc_", i),
    validBarcodes = cells,
    excludeChr = "",
    force = TRUE,
    minFrags = 1,
    addGeneScoreMat = FALSE,
    addTileMat = FALSE
  )
  elapsed.arrow <- as.numeric(Sys.time() - start.time, units = "secs")
  start.time <- Sys.time()
  proj <- ArchRProject(
    ArrowFiles = ArrowFiles,
    copyArrows = FALSE,
    showLogo = FALSE
  )
  elapsed.proj <- as.numeric(Sys.time() - start.time, units = "secs")
  saveArchRProject(ArchRProj = proj, outputDirectory = paste0("archr_pbmc/", i))
  # save timing
  write(
    x = paste0(elapsed.arrow, "\tArrow\tPBMC\t", i, "\n", elapsed.proj, "\tProject\tPBMC\t", i),
    file = "data/pbmc_atac/benchmarks/archr_object_creation.tsv",
    append = TRUE
  )
}

addArchRGenome("mm10")

for (i in downsamples_biccn) {
  # load the signac fragment file to get list of cells to include
  frags <- readRDS(paste0("/scratch/tim/biccn/downsampling/", i, ".rds"))
  cells <- Cells(frags)
  start.time <- Sys.time()
  ArrowFiles <- createArrowFiles(
    inputFiles = paste0("/scratch/tim//biccn/downsampling/", i, ".bed.gz"),
    sampleNames = paste0("biccn_", i),
    validBarcodes = cells,
    force = TRUE,
    excludeChr = "",
    minFrags = 1,
    addGeneScoreMat = FALSE,
    addTileMat = FALSE
  )
  elapsed.arrow <- as.numeric(Sys.time() - start.time, units = "secs")
  start.time <- Sys.time()
  proj <- ArchRProject(
    ArrowFiles = ArrowFiles,
    copyArrows = FALSE,
    showLogo = FALSE
  )
  elapsed.proj <- as.numeric(Sys.time() - start.time, units = "secs")
  saveArchRProject(ArchRProj = proj, outputDirectory = paste0("archr_biccn/", i))
  # save timing
  write(
    x = paste0(elapsed.arrow, "\tArrow\tBICCN\t", i, "\n", elapsed.proj, "\tProject\tBICCN\t", i),
    file = "data/biccn/benchmarks/archr_object_creation.tsv",
    append = TRUE
  )
}