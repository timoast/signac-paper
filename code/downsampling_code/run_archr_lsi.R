library(ArchR)
library(GenomicRanges)

args = commandArgs(trailingOnly = TRUE)
arrowfile <- args[1]
peakpath <- args[2]
nrep <- args[3]
timefile <- args[4]
genome <- args[5]

# load archr project for downsampling level
proj <- loadArchRProject(path = arrowfile)

# load peaks
peaks <- read.table(file = peakpath, sep = "\t", header = TRUE)
peaks <- makeGRangesFromDataFrame(peaks)
# remove chrM
peaks <- peaks[seqnames(peaks) != "chrM"]
message("Using ", length(peaks), " peaks")

# set threads
addArchRThreads(threads = 1)
addArchRGenome(genome)

proj <- addPeakSet(ArchRProj = proj, peakSet = peaks, force = TRUE)
proj <- addPeakMatrix(ArchRProj = proj, force = TRUE)

# run featurematrix n times
invisible(gc())
timings <- c()
for (i in seq_len(length.out = as.numeric(x = nrep))) {
  start.time <- Sys.time()
  proj <- addIterativeLSI(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    force = TRUE,
    sampleCellsPre = NULL
  )
  elapsed <- as.numeric(Sys.time() - start.time, units = "secs")
  timings <- c(timings, elapsed)
  invisible(gc())
}

# save timings
writeLines(text = sapply(X = timings, FUN = as.character), con = timefile, sep = "\n")
