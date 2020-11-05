library(ggplot2)
library(patchwork)

dimplot.mito <- readRDS("figures/mito_dimplot.rds") & theme(text = element_text(size = 10))
varplot.mito <- readRDS("figures/mito_varplot.rds") + theme(text = element_text(size = 10))
featplot.mito <- readRDS("figures/mito_allele_plot.rds") & labs(color = "Allele frequency") & theme(text = element_text(size = 10))
heatmap.mito <- readRDS("figures/mito_clone_hm.rds") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(text = element_text(size = 10))
covplots <- readRDS("figures/mito_covplot.rds") & theme(text = element_text(size = 10))

covplots <- covplots & theme(axis.text.x = element_text(size=4))
fig4 <- (dimplot.mito | varplot.mito | featplot.mito) / ((heatmap.mito | covplots) + plot_layout(widths = c(1, 2)))

ggsave(filename = "figures/figure4.png", plot = fig4, height = 8, width = 12, units = "in")