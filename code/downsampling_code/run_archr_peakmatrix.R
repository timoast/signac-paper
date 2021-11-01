library(ArchR)
library(GenomicRanges)

args = commandArgs(trailingOnly = TRUE)
ncore <- as.numeric(args[1])
arrowfile <- args[2]
peakpath <- args[3]
nrep <- args[4]
timefile <- args[5]
genome <- args[6]

# load archr project for downsampling level
proj <- loadArchRProject(path = arrowfile)

addArchRGenome(genome)

# load peaks
peaks <- read.table(file = peakpath, sep = "\t", header = TRUE)
peaks <- makeGRangesFromDataFrame(peaks)
# remove chrM
peaks <- peaks[seqnames(peaks) != "chrM"]

message("Using ", length(peaks), " peaks")

# set threads
addArchRThreads(threads = ncore)

# run featurematrix n times
invisible(gc())
timings <- c()
for (i in seq_len(length.out = as.numeric(x = nrep))) {
  start.time <- Sys.time()
  proj <- addPeakSet(ArchRProj = proj, peakSet = peaks, force = TRUE)
  proj <- addPeakMatrix(ArchRProj = proj, force = TRUE)
  elapsed <- as.numeric(Sys.time() - start.time, units = "secs")
  timings <- c(timings, elapsed)
  invisible(gc())
}

# save timings
writeLines(text = sapply(X = timings, FUN = as.character), con = timefile, sep = "\n")
