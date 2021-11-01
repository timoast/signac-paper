downsamples <- seq(from = 1000, to = 26000, by = 2000)
cores <- c(1, 2, 4, 8)

results_df <- data.frame()

runtime <- read.table(file = "data/pbmc_atac/benchmarks/signac_object_creation.tsv", sep = "\t")
runtime_archr <- read.table(file = "data/pbmc_atac/benchmarks/archr_object_creation.tsv", sep = "\t")
runtime_archr <- runtime_archr[runtime_archr$V2 == "Arrow", ]
result <- data.frame(
  "Cells" = c(runtime$V4, runtime_archr$V4),
  "Cores" = 1,
  "Step" = "Create",
  "Runtime" = c(runtime$V1, runtime_archr$V1),
  "Method" = c(rep("Signac", nrow(runtime)), rep("ArchR", nrow(runtime_archr)))
)
results_df <- rbind(results_df, result)

for (i in downsamples) {
  for (j in cores) {
    runtime <- readLines(con = paste0("data/pbmc_atac/benchmarks/featmat_runtime_", i, "_", j, ".txt"))
    rt_archr <- readLines(con = paste0("data/pbmc_atac/benchmarks/archr_featmat_runtime_", i, "_", j, ".txt"))
    runtime <- sapply(runtime, as.numeric, USE.NAMES = FALSE)
    rt_archr <- sapply(rt_archr, as.numeric, USE.NAMES = FALSE)
    result <- data.frame(
      "Cells" = i,
      "Cores" = j,
      "Step" = "FeatureMatrix",
      "Runtime" = c(runtime, rt_archr),
      "Method" = c(rep("Signac", length(runtime)), rep("ArchR", length(rt_archr)))
    )
    results_df <- rbind(results_df, result)
  }
}

for (i in downsamples) {
  runtime <- readLines(con = paste0("data/pbmc_atac/benchmarks/nucleosome_runtime_", i, ".txt"))
  runtime <- sapply(runtime, as.numeric, USE.NAMES = FALSE)
  result <- data.frame(
    "Cells" = i,
    "Cores" = 1,
    "Step" = "NucleosomeSignal",
    "Runtime" = runtime,
    "Method" = "Signac"
  )
  results_df <- rbind(results_df, result)
}

for (i in downsamples) {
  for (j in cores) {
    runtime <- readLines(con = paste0("data/pbmc_atac/benchmarks/ga_runtime_", i, "_", j, ".txt"))
    rt_archr <- readLines(con = paste0("data/pbmc_atac/benchmarks/archr_geneactivity_runtime_", i, "_", j, ".txt"))
    runtime <- sapply(runtime, as.numeric, USE.NAMES = FALSE)
    rt_archr <- sapply(rt_archr, as.numeric, USE.NAMES = FALSE)
    result <- data.frame(
      "Cells" = i,
      "Cores" = j,
      "Step" = "GeneActivity",
      "Runtime" = c(runtime, rt_archr),
      "Method" = c(rep("Signac", length(runtime)), rep("ArchR", length(rt_archr)))
    )
    results_df <- rbind(results_df, result)
  }
}

for (i in downsamples) {
  for (j in cores) {
    runtime <- readLines(con = paste0("data/pbmc_atac/benchmarks/tss_runtime_", i, "_", j, ".txt"))
    runtime <- sapply(runtime, as.numeric, USE.NAMES = FALSE)
    result <- data.frame(
      "Cells" = i,
      "Cores" = j,
      "Step" = "TSSEnrichment",
      "Runtime" = runtime,
      "Method" = "Signac"
    )
    results_df <- rbind(results_df, result)
  }
}

for (i in downsamples) {
  runtime <- readLines(con = paste0("data/pbmc_atac/benchmarks/tfidf_runtime_", i, ".txt"))
  runtime <- sapply(runtime, as.numeric, USE.NAMES = FALSE)
  result <- data.frame(
    "Cells" = i,
    "Cores" = 1,
    "Step" = "RunTFIDF",
    "Runtime" = runtime,
    "Method" = "Signac"
  )
  results_df <- rbind(results_df, result)
}

for (i in downsamples) {
  runtime <- readLines(con = paste0("data/pbmc_atac/benchmarks/svd_runtime_", i, ".txt"))
  runtime <- sapply(runtime, as.numeric, USE.NAMES = FALSE)
  result <- data.frame(
    "Cells" = i,
    "Cores" = 1,
    "Step" = "RunSVD",
    "Runtime" = runtime,
    "Method" = "Signac"
  )
  results_df <- rbind(results_df, result)
}

for (i in downsamples) {
  runtime <- readLines(con = paste0("data/pbmc_atac/benchmarks/archr_lsi_runtime_", i, ".txt"))
  runtime <- sapply(runtime, as.numeric, USE.NAMES = FALSE)
  result <- data.frame(
    "Cells" = i,
    "Cores" = 1,
    "Step" = "LSI",
    "Runtime" = runtime,
    "Method" = "ArchR"
  )
  results_df <- rbind(results_df, result)
}

for (i in downsamples) {
  runtime <- readLines(con = paste0("data/pbmc_atac/benchmarks/archr_est_lsi_runtime_", i, ".txt"))
  runtime <- sapply(runtime, as.numeric, USE.NAMES = FALSE)
  result <- data.frame(
    "Cells" = i,
    "Cores" = 1,
    "Step" = "estLSI",
    "Runtime" = runtime,
    "Method" = "ArchR"
  )
  results_df <- rbind(results_df, result)
}

# add LSI
tfidf <- results_df[results_df$Step == "RunTFIDF", ]
runsvd <- results_df[results_df$Step == "RunSVD", ]
lsi <- runsvd
lsi$Runtime <- runsvd$Runtime + tfidf$Runtime
lsi$Step <- "LSI"
results_df <- rbind(results_df, lsi)

results_df$Cells <- as.numeric(results_df$Cells)
results_df$Cores <- as.factor(results_df$Cores)

write.table(x = results_df, file = "data/pbmc_atac/timings.tsv")
