library(Seurat)
library(Signac)

args = commandArgs(trailingOnly = TRUE)
objpath <- args[1]
nrep <- args[2]
timefile <- args[3]
outfile <- args[4]

# load object
obj <- readRDS(file = objpath)

# run svd n times
invisible(gc())
timings <- c()
for (i in seq_len(length.out = as.numeric(x = nrep))) {
  start.time <- Sys.time()
  obj <- RunSVD(obj, features = rownames(obj))
  elapsed <- as.numeric(Sys.time() - start.time, units = "secs")
  timings <- c(timings, elapsed)
  invisible(gc())
}

# save timings
writeLines(text = sapply(X = timings, FUN = as.character), con = timefile, sep = "\n")

# save object for next step
saveRDS(object = obj, file = outfile, version = 2)