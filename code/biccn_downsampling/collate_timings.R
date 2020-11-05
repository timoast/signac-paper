downsamples <- c('50000', '100000', '200000', '300000', '400000', '500000', '600000', '700000')
cores <- c(1, 2, 4, 8)

results_df <- data.frame()
for (i in downsamples) {
  for (j in cores) {
    max_resident_mem <- readLines(con = paste0("data/biccn/benchmarks/featmat_mem_", i, "_", j, ".txt"))[10]
    max_resident_mem <- strsplit(max_resident_mem, "\tMaximum resident set size (kbytes): ", fixed = TRUE)[[1]][[2]]
    runtime <- readLines(con = paste0("data/biccn/benchmarks/featmat_runtime_", i, "_", j, ".txt"))
    runtime <- sapply(runtime, as.numeric, USE.NAMES = FALSE)
    result <- data.frame(
      "Cells" = i,
      "Cores" = j,
      "Step" = "FeatureMatrix",
      "Memory" = max_resident_mem,
      "Runtime" = runtime
    )
    results_df <- rbind(results_df, result)
  }
}

for (i in downsamples) {
  max_resident_mem <- readLines(con = paste0("data/biccn/benchmarks/nucleosome_mem_", i, ".txt"))[10]
  max_resident_mem <- strsplit(max_resident_mem, "\tMaximum resident set size (kbytes): ", fixed = TRUE)[[1]][[2]]
  runtime <- readLines(con = paste0("data/biccn/benchmarks/nucleosome_runtime_", i, ".txt"))
  runtime <- sapply(runtime, as.numeric, USE.NAMES = FALSE)
  result <- data.frame(
    "Cells" = i,
    "Cores" = 1,
    "Step" = "NucleosomeSignal",
    "Memory" = max_resident_mem,
    "Runtime" = runtime
  )
  results_df <- rbind(results_df, result)
}

for (i in downsamples) {
  max_resident_mem <- readLines(con = paste0("data/biccn/benchmarks/tss_mem_", i, ".txt"))[10]
  max_resident_mem <- strsplit(max_resident_mem, "\tMaximum resident set size (kbytes): ", fixed = TRUE)[[1]][[2]]
  runtime <- readLines(con = paste0("data/biccn/benchmarks/tss_runtime_", i, ".txt"))
  runtime <- sapply(runtime, as.numeric, USE.NAMES = FALSE)
  result <- data.frame(
    "Cells" = i,
    "Cores" = 1,
    "Step" = "TSSEnrichment",
    "Memory" = max_resident_mem,
    "Runtime" = runtime
  )
  results_df <- rbind(results_df, result)
}

for (i in downsamples) {
  max_resident_mem <- readLines(con = paste0("data/biccn/benchmarks/tfidf_mem_", i, ".txt"))[10]
  max_resident_mem <- strsplit(max_resident_mem, "\tMaximum resident set size (kbytes): ", fixed = TRUE)[[1]][[2]]
  runtime <- readLines(con = paste0("data/biccn/benchmarks/tfidf_runtime_", i, ".txt"))
  runtime <- sapply(runtime, as.numeric, USE.NAMES = FALSE)
  result <- data.frame(
    "Cells" = i,
    "Cores" = 1,
    "Step" = "RunTFIDF",
    "Memory" = max_resident_mem,
    "Runtime" = runtime
  )
  results_df <- rbind(results_df, result)
}

for (i in downsamples) {
  max_resident_mem <- readLines(con = paste0("data/biccn/benchmarks/svd_mem_", i, ".txt"))[10]
  max_resident_mem <- strsplit(max_resident_mem, "\tMaximum resident set size (kbytes): ", fixed = TRUE)[[1]][[2]]
  runtime <- readLines(con = paste0("data/biccn/benchmarks/svd_runtime_", i, ".txt"))
  runtime <- sapply(runtime, as.numeric, USE.NAMES = FALSE)
  result <- data.frame(
    "Cells" = i,
    "Cores" = 1,
    "Step" = "RunSVD",
    "Memory" = max_resident_mem,
    "Runtime" = runtime
  )
  results_df <- rbind(results_df, result)
}

results_df$Cells <- as.numeric(results_df$Cells)
results_df$Cores <- as.factor(results_df$Cores)
results_df$Memory <- as.numeric(results_df$Memory)

write.table(x = results_df, file = "data/biccn/timings.tsv")
