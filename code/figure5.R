library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)

biccn <- readRDS("objects/biccn.rds")
timings <- read.table("data/biccn/timings.tsv")
timings$Cores <- as.factor(timings$Cores)

# featurematrix timings
runtime_fmat <- ggplot(data = timings[timings$Step == "FeatureMatrix", ], mapping = aes(x = Cells, y = Runtime / 60, color = Cores)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  ylab("Runtime (minutes)") +
  scale_x_continuous(labels = scales::comma, breaks = seq(0, 700000, 200000)) +
  scale_y_continuous(breaks = seq(0, 300, 60)) +
  theme_bw() +
  ggtitle(label = "FeatureMatrix", subtitle = "315,334 peaks")

mem_fmat <- ggplot(data = timings[timings$Step == "FeatureMatrix", ], mapping = aes(x = Cells, y = Memory / 1024^2, color = Cores)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  ylab("Max memory (Gb)") +
  scale_x_continuous(labels = scales::comma, breaks = seq(0, 700000, 200000)) +
  scale_y_continuous(breaks = seq(0, 150, 25)) +
  theme_bw()

fmat <- (runtime_fmat | mem_fmat) + plot_layout(guides = "collect")

# nucleosome signal timings
runtime_ns <- ggplot(data = timings[timings$Step == "NucleosomeSignal", ], mapping = aes(x = Cells, y = Runtime / 60)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  ylab("Runtime (minutes)") +
  scale_x_continuous(labels = scales::comma, breaks = seq(0, 700000, 200000)) +
  theme_bw() +
  ggtitle("NucleosomeSignal")

mem_ns <- ggplot(data = timings[timings$Step == "NucleosomeSignal", ], mapping = aes(x = Cells, y = Memory / 1024^2)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  ylab("Max memory (Gb)") +
  scale_x_continuous(labels = scales::comma, breaks = seq(0, 700000, 100000)) +
  theme_bw()

ns <- runtime_ns | mem_ns

# tss enrichment timings
runtime_tss <- ggplot(data = timings[timings$Step == "TSSEnrichment", ], mapping = aes(x = Cells, y = Runtime / 60)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  ylab("Runtime (minutes)") +
  scale_x_continuous(labels = scales::comma, breaks = seq(0, 700000, 200000)) +
  theme_bw() +
  ggtitle("TSSEnrichment")

mem_tss <- ggplot(data = timings[timings$Step == "TSSEnrichment", ], mapping = aes(x = Cells, y = Memory / 1024^2)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  ylab("Max memory (Gb)") +
  scale_x_continuous(labels = scales::comma, breaks = seq(0, 700000, 200000)) +
  theme_bw()

tss <- runtime_tss | mem_tss

# tf-idf timings
runtime_tfidf <- ggplot(data = timings[timings$Step == "RunTFIDF", ], mapping = aes(x = Cells, y = Runtime / 60)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  ylab("Runtime (minutes)") +
  scale_x_continuous(labels = scales::comma, breaks = seq(0, 700000, 200000)) +
  theme_bw() +
  ggtitle("RunTFIDF")

mem_tfidf <- ggplot(data = timings[timings$Step == "RunTFIDF", ], mapping = aes(x = Cells, y = Memory / 1024^2)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  ylab("Max memory (Gb)") +
  scale_x_continuous(labels = scales::comma, breaks = seq(0, 700000, 200000)) +
  theme_bw()

tfidf <- runtime_tfidf | mem_tfidf

# svd timings
runtime_svd <- ggplot(data = timings[timings$Step == "RunSVD", ], mapping = aes(x = Cells, y = Runtime / 60)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  ylab("Runtime (minutes)") +
  scale_x_continuous(labels = scales::comma, breaks = seq(0, 700000, 200000)) +
  theme_bw() +
  ggtitle("RunSVD")

mem_svd <- ggplot(data = timings[timings$Step == "RunSVD", ], mapping = aes(x = Cells, y = Memory / 1024^2)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  ylab("Max memory (Gb)") +
  scale_x_continuous(labels = scales::comma, breaks = seq(0, 700000, 200000)) +
  theme_bw()

runsvd <- runtime_svd | mem_svd

# max memory use for steps other than FeatureMatrix just refects the difference in
# object size for different number of cells, rather than scaling of the function
# itself

# collate all runtime panels
runtimes <- fmat / (runtime_ns | runtime_tss) / (runtime_tfidf | runtime_svd)

# biccn dimplot
dp <- DimPlot(biccn, group.by = "broad_region", pt.size = 0.1)

ggsave(filename = "figures/figure5a.png", plot = dp, height = 12, width = 14, dpi = 200)
ggsave(filename = "figures/figure5b.png", plot = runtimes, height = 8, width = 6, dpi = 300)
