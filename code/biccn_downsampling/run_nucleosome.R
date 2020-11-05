library(Signac)
library(Seurat)

args = commandArgs(trailingOnly = TRUE)
countpath <- args[1]
fragpath <- args[2]
annotpath <- args[3]
nrep <- args[4]
timefile <- args[5]
outfile <- args[6]

# load the count matrix and fragment object
counts <- readRDS(file = countpath)
frags <- readRDS(file = fragpath)
annot <- readRDS(file = annotpath)

# create Seurat object
obj <- CreateSeuratObject(
  counts = CreateChromatinAssay(
    counts = counts,
    fragments = frags,
    annotation = annot,
    validate.fragments = FALSE
  ),
  assay = "peaks",
  project = "BICCN"
)

# run nucleosome signal n times
invisible(gc())
timings <- c()
for (i in seq_len(length.out = as.numeric(x = nrep))) {
  start.time <- Sys.time()
  obj <- NucleosomeSignal(obj)
  elapsed <- as.numeric(Sys.time() - start.time, units = "secs")
  timings <- c(timings, elapsed)
  invisible(gc())
}

# save timings
writeLines(text = sapply(X = timings, FUN = as.character), con = timefile, sep = "\n")

# save object for next step
saveRDS(object = obj, file = outfile, version = 2)
