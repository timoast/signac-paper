library(ArchR)

args = commandArgs(trailingOnly = TRUE)
ncore <- as.numeric(args[1])
arrowfile <- args[2]
nrep <- args[3]
timefile <- args[4]
genome <- args[5]

# load archr project for downsampling level
proj <- loadArchRProject(path = arrowfile)

# set threads
addArchRThreads(threads = ncore)
addArchRGenome(genome)

# run gene activity n times
invisible(gc())
timings <- c()
for (i in seq_len(length.out = as.numeric(x = nrep))) {
  start.time <- Sys.time()
  proj <- addGeneScoreMatrix(input = proj, force = TRUE)
  elapsed <- as.numeric(Sys.time() - start.time, units = "secs")
  timings <- c(timings, elapsed)
  invisible(gc())
}

# save timings
writeLines(text = sapply(X = timings, FUN = as.character), con = timefile, sep = "\n")