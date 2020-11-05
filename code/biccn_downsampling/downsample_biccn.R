library(Signac)
library(Rsamtools)
set.seed(1234)

# load the full BICCN dataset
biccn <- readRDS("objects/biccn.rds")

# randomly sample different numbers of cells
downsamples <- c(50000, 100000, 200000, 300000, 400000, 500000, 600000, 700000)

cells.select <- sapply(X = downsamples, FUN = function(x) {
  sample(x = colnames(x = biccn), replace = FALSE, size = x)
})

# create downsampled fragment files
for (i in seq_along(along.with = downsamples)) {
  ds <- format(x = downsamples[[i]], scientific = FALSE)
  outfile <- paste0("data/biccn/downsampling/", ds, ".bed.gz")
  frag.dest <- paste0("data/biccn/downsampling/", ds, ".rds")
  FilterCells(
    fragments = "data/biccn/fragments.sort.bed.gz",
    cells = cells.select[[i]],
    outfile = outfile,
    verbose = TRUE
  )
  frags <- CreateFragmentObject(
    path = outfile,
    cells = cells.select[[i]],
    validate.fragments = FALSE
  )
  saveRDS(object = frags, file = frag.dest, version = 2)
}