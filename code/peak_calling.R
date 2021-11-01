library(Seurat)
library(Signac)
library(ggplot2)
library(patchwork)
library(GenomicRanges)


pbmc <- readRDS("objects/pbmc.rds")
ident.use <- "CD14 Mono"

obj <- pbmc[, Idents(pbmc) == ident.use]

ds.level <- seq(50, 2850, 100)
pk.list <- list()
cp.list <- list()
for (i in seq_along(ds.level)) {
  set.seed(1234)
  cells.use <- sample(x = colnames(obj), size = ds.level[[i]], replace = FALSE)
  obj.ds <- obj[, cells.use]
  obj.ds$ds <- paste0(ds.level[[i]], " cells")
  pk <- CallPeaks(
    object = obj.ds,
    group.by = "celltype",
    additional.args = "--max-gap 50"
  )
  cp.ds <- CoveragePlot(
    object = obj.ds,
    region = "LYZ",
    group.by = "ds",
    extend.upstream = 6000,
    extend.downstream = 8000,
    peaks = FALSE,
    ranges = pk,
    ymax = 260
  )
  pk.list[[as.character(i)]] <- pk
  cp.list[[as.character(i)]] <- cp.ds
}

# find overlaps with highest sampling
pk.highest <- pk.list[[length(pk.list)]]

olap <- c()
for (i in seq_along(pk.list)) {
  ol <- sum(countOverlaps(query = pk.list[[i]], pk.highest)) / length(pk.highest)
  olap[[i]] <- ol
}

df <- data.frame(x = unlist(olap), cells = ds.level)
p <- ggplot(data = df, mapping = aes(x = cells, y = x)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  xlab("Number of cells") +
  ylab("Fraction of peaks recovered") +
  theme_bw() +
  ylim(c(0, 1)) +
  scale_x_continuous(breaks=seq(0, 2850, 200))

cp.use <- cp.list[c(29, 10, 1)]
p2 <- wrap_plots(cp.use, ncol = 1)
fig <- (p | p2) + plot_layout(heights = c(1, 3))
ggsave(filename = "figures/peakcalls.png", plot = fig, height = 10, width = 16)
