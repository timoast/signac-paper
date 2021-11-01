library(Signac)
set.seed(1234)

# load the full dataset
biccn <- readRDS("objects/biccn.rds")
pbmc <- readRDS("objects/pbmc_atac.rds")

# randomly sample different numbers of cells
downsamples_biccn <- c(50000, 100000, 200000, 300000, 400000, 500000, 600000, 700000)
downsamples_pbmc <- seq(from = 1000, to = ncol(pbmc), by = 2000)

cells.biccn <- sapply(X = downsamples_biccn, FUN = function(x) {
  sample(x = colnames(x = biccn), replace = FALSE, size = x)
})

cells.pbmc <- sapply(X = downsamples_pbmc, FUN = function(x) {
  sample(x = colnames(x = pbmc), replace = FALSE, size = x)
})

downsample_fragments <- function(fragpath, downsamples, cells, outpath, timepath, project) {
  for (i in seq_along(along.with = downsamples)) {
    ds <- format(x = downsamples[[i]], scientific = FALSE)
    outfile <- paste0(outpath, ds, ".bed.gz")
    frag.dest <- paste0(outpath, ds, ".rds")
    FilterCells(
      fragments = fragpath,
      cells = cells[[i]],
      outfile = outfile,
      verbose = TRUE
    )
    time.start <- Sys.time()
    frags <- CreateFragmentObject(
      path = outfile,
      cells = cells[[i]]
    )
    elapsed <- as.numeric(Sys.time() - time.start, units = "secs")
    saveRDS(object = frags, file = frag.dest, version = 2)
    write(
      x = paste0(elapsed, "\tSignac\t", project, "\t", ds),
      file = timepath,
      append = TRUE
    )
  }
}

downsample_fragments(
  fragpath = "data/biccn/fragments.bed.gz",
  outpath = "/scratch/tim/biccn/downsampling/",
  timepath = "data/biccn/benchmarks/signac_object_creation.tsv",
  project = "BICCN",
  cells = cells.biccn,
  downsamples = downsamples_biccn
)

downsample_fragments(
  fragpath = "data/pbmc_atac/fragments.bed.gz",
  outpath = "/scratch/tim/pbmc_atac/downsampling/",
  timepath = "data/pbmc_atac/benchmarks/signac_object_creation.tsv",
  project = "PBMC",
  cells = cells.pbmc,
  downsamples = downsamples_pbmc
)
