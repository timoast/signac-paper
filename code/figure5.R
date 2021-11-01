library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(paletteer)

biccn <- readRDS("objects/biccn.rds")
pbmc <- readRDS("objects/pbmc_atac.rds")

timings_biccn <- read.table("data/biccn/timings.tsv")
timings_pbmc <- read.table("data/pbmc_atac/timings.tsv")

timings_pbmc$Cores <- as.factor(timings_pbmc$Cores)
timings_biccn$Cores <- as.factor(timings_biccn$Cores)

colors.use <- paletteer_d("ggthemes::Tableau_10")
signac_color <- colors.use[[1]]
archr_color <- colors.use[[9]]

# dimplots
dp <- DimPlot(biccn, group.by = "MajorType", label = TRUE, repel = TRUE, raster = FALSE, pt.size = 0.1) +
  theme_classic() +
  NoLegend() +
  ggtitle(label = "Adult mouse brain", subtitle = "734,000 nuclei")

dp_batch <- DimPlot(biccn, group.by = "orig.ident", label = FALSE, raster = FALSE, pt.size = 0.1, shuffle = TRUE) +
  theme_classic() +
  ggtitle(label = "Adult mouse brain", subtitle = "734,000 nuclei")

pbmc$all <- "PBMC"
dp_pbmc <- DimPlot(pbmc, group.by = "all", label = FALSE) +
  theme_classic() +
  NoLegend() +
  ggtitle(label = "Human PBMCs", subtitle = "26,579 nuclei")

########## PBMC ##########
# object creation
runtime_create <- ggplot(data = timings_pbmc[timings_pbmc$Step == "Create", ], mapping = aes(x = Cells, y = Runtime/60, color = Cores)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  facet_wrap(~Method) +
  ylab("Runtime (minutes)") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_x_continuous(labels = scales::comma, breaks = seq(0, 25000, 5000)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ggtitle("Object creation")

# featurematrix timings
runtime_fmat <- ggplot(data = timings_pbmc[timings_pbmc$Step == "FeatureMatrix", ], mapping = aes(x = Cells, y = Runtime / 60, color = Cores)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  facet_wrap(~Method) +
  ylab("Runtime (minutes)") +
  theme_bw() +
  scale_x_continuous(labels = scales::comma, breaks = seq(0, 25000, 5000)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ggtitle(label = "FeatureMatrix", subtitle = "160,906 peaks")

# gene activity
runtime_ga <- ggplot(data = timings_pbmc[timings_pbmc$Step == "GeneActivity", ], mapping = aes(x = Cells, y = Runtime / 60, color = Cores)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  facet_wrap(~Method) +
  ylab("Runtime (minutes)") +
  theme_bw() +
  scale_x_continuous(labels = scales::comma, breaks = seq(0, 25000, 5000)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ggtitle(label = "GeneActivity")

# nucleosome signal timings
runtime_ns <- ggplot(data = timings_pbmc[timings_pbmc$Step == "NucleosomeSignal", ], mapping = aes(x = Cells, y = Runtime / 60)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  ylab("Runtime (minutes)") +
  theme_bw() +
  scale_x_continuous(labels = scales::comma, breaks = seq(0, 25000, 5000)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ggtitle("NucleosomeSignal")

# QC timings
runtime_qc <- ggplot(
  data = timings_pbmc[timings_pbmc$Step %in% c("NucleosomeSignal", "TSSEnrichment"), ],
  mapping = aes(x = Cells, y = Runtime / 60, color = Cores)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  facet_wrap(~Step, scales = "free_y") +
  ylab("Runtime (minutes)") +
  theme_bw() +
  scale_x_continuous(labels = scales::comma, breaks = seq(0, 25000, 5000)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ggtitle("Signac quality control metrics")

# LSI timings
runtime_lsi <- ggplot(data = timings_pbmc[timings_pbmc$Step == "LSI", ], mapping = aes(x = Cells, y = Runtime / 60, color = Cores)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  facet_wrap(~Method) +
  ylab("Runtime (minutes)") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_x_continuous(labels = scales::comma, breaks = seq(0, 25000, 5000)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ggtitle("LSI")

# total runtime
archr_pbmc_total <- as.numeric(readLines(con = "data/pbmc_atac/archr_total_runtime.txt"))
signac_pbmc_total <- as.numeric(readLines(con = "data/pbmc_atac/signac_total_runtime.txt"))

total_pbmc <- ggplot(data = data.frame(Method = c("Signac", "ArchR"), runtime = c(signac_pbmc_total/60, archr_pbmc_total/60)),
                     mapping = aes(y = runtime, x = Method, fill = Method)) +
  geom_bar(stat = "identity") +
  ylab("Runtime (minutes)") +
  theme_bw() +
  ggtitle(label = "Total runtime", subtitle = "26,579 nuclei; 8 cores") +
  scale_fill_manual(values = c(archr_color, signac_color)) +
  theme(legend.position = "none")

# collate all runtime panels
runtimes_pbmc <- ((runtime_create / runtime_fmat) | (runtime_ga / runtime_qc) | (runtime_lsi / total_pbmc)) + plot_layout(guides = "collect")

######### BICCN ###########
# object creation
runtime_create <- ggplot(data = timings_biccn[timings_biccn$Step == "Create", ], mapping = aes(x = Cells, y = Runtime/60, color = Cores)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  facet_wrap(~Method) +
  scale_x_continuous(labels = scales::comma, breaks = seq(0, 700000, 200000)) +
  ylab("Runtime (minutes)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(legend.position = "none") +
  ggtitle("Object creation")

# featurematrix timings
runtime_fmat <- ggplot(data = timings_biccn[timings_biccn$Step == "FeatureMatrix", ], mapping = aes(x = Cells, y = Runtime / 60, color = Cores)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  facet_wrap(~Method) +
  ylab("Runtime (minutes)") +
  scale_x_continuous(labels = scales::comma, breaks = seq(0, 700000, 200000)) +
  scale_y_continuous(breaks = seq(0, 300, 60)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ggtitle(label = "FeatureMatrix", subtitle = "263,815 peaks")

# gene activity
runtime_ga <- ggplot(data = timings_biccn[timings_biccn$Step == "GeneActivity", ], mapping = aes(x = Cells, y = Runtime / 60, color = Cores)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  facet_wrap(~Method) +
  ylab("Runtime (minutes)") +
  scale_x_continuous(labels = scales::comma, breaks = seq(0, 700000, 200000)) +
  scale_y_continuous(breaks = seq(0, 300, 60)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ggtitle(label = "GeneActivity")

# QC step
runtime_qc <- ggplot(
  data = timings_biccn[timings_biccn$Step %in% c("NucleosomeSignal", "TSSEnrichment"), ],
  mapping = aes(x = Cells, y = Runtime / 60, color = Cores)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  facet_wrap(~Step) +
  ylab("Runtime (minutes)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_x_continuous(labels = scales::comma, breaks = seq(0, 700000, 200000)) +
  ggtitle("Signac quality control metrics")

# LSI timings
runtime_lsi <- ggplot(data = timings_biccn[timings_biccn$Step == "LSI", ], mapping = aes(x = Cells, y = Runtime / 60, color = Cores)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  facet_wrap(~Method) +
  ylab("Runtime (minutes)") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_x_continuous(labels = scales::comma, breaks = seq(0, 700000, 200000)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ggtitle("LSI")

# total runtime
archr_biccn_total <- as.numeric(readLines(con = "data/biccn/archr_total_runtime.txt"))
signac_biccn_total <- as.numeric(readLines(con = "data/biccn/signac_total_runtime.txt"))

total_biccn <- ggplot(data = data.frame(Method = c("Signac", "ArchR"), runtime = c(signac_biccn_total/60/60, archr_biccn_total/60/60)),
                     mapping = aes(y = runtime, x = Method, fill = Method)) +
  geom_bar(stat = "identity") +
  ylab("Runtime (hours)") +
  theme_bw() +
  scale_y_continuous(breaks = c(0,1,2,3,4,5)) +
  ggtitle(label = "Total runtime", subtitle = "734,000 nuclei; 8 cores") +
  scale_fill_manual(values = c(archr_color, signac_color)) +
  theme(legend.position = "none")

# collate all runtime panels
runtimes_biccn <- ((runtime_create / runtime_fmat) | (runtime_ga / runtime_qc) | (runtime_lsi / total_biccn)) + plot_layout(guides = "collect")

# save plots
ggsave(filename = "figures/figure5a.png", plot = runtimes_pbmc, height = 6, width = 12, dpi = 500)
ggsave(filename = "figures/figure5b.png", plot = runtimes_biccn, height = 6, width = 12, dpi = 500)
ggsave(filename = "figures/biccn_dimplot_batch.png", plot = dp_batch, height = 12, width = 18, dpi = 300)
ggsave(filename = "figures/biccn_dimplot.png", plot = dp, height = 6, width = 5, dpi = 300)
ggsave(filename = "figures/pbmc_atac_dimplot.png", plot = dp_pbmc, height = 6, width = 5, dpi = 300)
