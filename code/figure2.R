library(ggplot2)
library(patchwork)

# load figures
tss_plot <- readRDS("figures/tss_enrichment.rds")
nucleosome_plot <- readRDS("figures/nucleosome_signal.rds") + ggtitle("Nucleosome signal")
joint_dimplot <- readRDS("figures/pbmc_join_dimplot.rds") + NoLegend() + xlab("UMAP 1") + ylab("UMAP 2")
mp <- readRDS("figures/motifplot.rds")
tf_chromvar <- readRDS("figures/chromvar_vln.rds")
tf_expression <- readRDS("figures/tf_rna_vln.rds")
fp <- readRDS("figures/footprint.rds")
covplot <- readRDS("figures/coverage_plot.rds")
gene_per_link_plot <- readRDS("figures/genes_per_link_plot.rds")
link_per_gene_plot <- readRDS("figures/link_per_gene_plot.rds")
distplot_positive <- readRDS("figures/distance_positive.rds")
distplot_negative <- readRDS("figures/distance_negative.rds")

tf_chromvar <- tf_chromvar & theme(text = element_text(size = 12), axis.text = element_text(size = 12))
tf_expression <- tf_expression & theme(text = element_text(size = 12), axis.text = element_text(size = 12))
lnkplot <- gene_per_link_plot / link_per_gene_plot
distances <- distplot_positive / distplot_negative
top.panel <- ((nucleosome_plot / tss_plot) | joint_dimplot | covplot) + plot_layout(widths = c(1, 3, 4)) &
  theme(text = element_text(size = 12), axis.text = element_text(size = 12))
panel1 <- (mp / wrap_plots(tf_chromvar) / tf_expression & xlab("")) & theme(text = element_text(size = 12))

ggsave(filename = "figures/figure2.png", plot = top.panel, width = 16, height = 5, units = "in")
ggsave(filename = "figures/figure2_2.png", plot = panel1, width = 6, height = 6, units = "in")
ggsave(filename = "figures/figure2_3.png", plot = fp, width = 4, height = 6, units = 'in')
ggsave(filename = "figures/figure2_4.png", plot = lnkplot, width = 5, height = 6, units = 'in')
ggsave(filename = "figures/figure2_5.png", plot = distances, width = 4, height = 6, units = 'in')

linked_1 <- readRDS("figures/linked_covplot1.rds")
linked_2 <- readRDS("figures/linked_covplot2.rds")

lower.panel <- (linked_1 | linked_2) &
  theme(text = element_text(size = 10), axis.text = element_text(size = 10))

ggsave(filename = "figures/figure2_6.png", plot = lower.panel, width = 16, height = 5, units = "in")

## supplementary figure 1
cp <- readRDS("figures/cellranger_peakcalling.rds")
missed <- readRDS("figures/macs2_pseudobulk.rds")
missed_peak_count <- readRDS("figures/missed_peak_count.rds")
supfig1 <- missed | (missed_peak_count / plot_spacer()) | cp
ggsave(filename = "figures/figs1.png", plot = supfig1, height = 8, width = 15, units = "in")

## supplementary figure 2
atac_dimplot <- readRDS("figures/pbmc_atac_dimplot.rds") + NoLegend() + xlab("UMAP 1") + ylab("UMAP 2")
rna_dimplot <- readRDS("figures/pbmc_rna_dimplot.rds") + NoLegend() + xlab("UMAP 1") + ylab("UMAP 2")
supfig2 <- ((atac_dimplot + ggtitle("DNA accessibility")) | (rna_dimplot + ggtitle("Gene expression"))) +
  plot_layout(guides = "collect")
ggsave(filename = "figures/sup_figure_2.png", plot = supfig2, height = 8, width = 16)
