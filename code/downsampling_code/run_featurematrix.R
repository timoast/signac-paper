library(Signac)
library(GenomicRanges)

args = commandArgs(trailingOnly = TRUE)
ncore <- as.numeric(args[1])
fragpath <- args[2]
peakpath <- args[3]
nrep <- args[4]
timefile <- args[5]
outfile <- args[6]

message("Using ", ncore, " cores")
if (ncore > 1) {
  library(future)
  plan("multicore", workers = ncore)
  options(future.globals.maxSize = 100 * 1024^3)
}

# load fragment object
frags <- readRDS(file = fragpath)

# load peaks
peaks <- read.table(file = peakpath, sep = "\t", header = TRUE)
peaks <- makeGRangesFromDataFrame(peaks)
message("Using ", length(peaks), " peaks")

# run featurematrix n times
invisible(gc())
timings <- c()
for (i in seq_len(length.out = as.numeric(x = nrep))) {
  start.time <- Sys.time()
  fmat <- FeatureMatrix(fragments = frags, features = peaks, cells = Cells(frags))
  elapsed <- as.numeric(Sys.time() - start.time, units = "secs")
  timings <- c(timings, elapsed)
  invisible(gc())
}

# save timings
writeLines(text = sapply(X = timings, FUN = as.character), con = timefile, sep = "\n")

# save last count matrix for next step
saveRDS(object = fmat, file = outfile, version = 2)