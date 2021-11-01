library(Seurat)
library(Signac)

args = commandArgs(trailingOnly = TRUE)
objpath <- args[1]
ncore <- as.numeric(args[2])
nrep <- args[3]
timefile <- args[4]

# load object
obj <- readRDS(file = objpath)


message("Using ", ncore, " cores")
if (ncore > 1) {
  library(future)
  plan("multicore", workers = ncore)
  options(future.globals.maxSize = +Inf)
}

# run svd n times
invisible(gc())
timings <- c()
for (i in seq_len(length.out = as.numeric(x = nrep))) {
  start.time <- Sys.time()
  rna <- GeneActivity(obj, process_n = 1000)
  elapsed <- as.numeric(Sys.time() - start.time, units = "secs")
  timings <- c(timings, elapsed)
  invisible(gc())
}

# save timings
writeLines(text = sapply(X = timings, FUN = as.character), con = timefile, sep = "\n")