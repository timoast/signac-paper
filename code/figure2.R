library(ggplot2)
library(patchwork)

# load figures
qc_dist <- readRDS("figures/qc_dist.rds")
tss_plot <- readRDS("figures/tss_enrichment.rds")
nucleosome_plot <- readRDS("figures/nucleosome_signal.rds") + ggtitle("Nucleosome signal")
atac_dimplot <- readRDS("figures/pbmc_atac_dimplot.rds")
mp <- readRDS("figures/motifplot.rds")
tf_chromvar <- readRDS("figures/chromvar_vln.rds")
tf_expression <- readRDS("figures/tf_rna_vln.rds")
fp <- readRDS("figures/footprint.rds")
gene_per_link_plot <- readRDS("figures/genes_per_link_plot.rds")
link_per_gene_plot <- readRDS("figures/link_per_gene_plot.rds")
distplot_positive <- readRDS("figures/distance_positive.rds")
distplot_negative <- readRDS("figures/distance_negative.rds")
pval_dist <- readRDS("figures/link_pvals.rds")
label_transfer_accuracy <- readRDS("figures/label_transfer_accuracy.rds")

tf_chromvar <- tf_chromvar & theme(text = element_text(size = 12), axis.text = element_text(size = 12))
tf_expression <- tf_expression & theme(text = element_text(size = 12), axis.text = element_text(size = 12))
lnkplot <- gene_per_link_plot / link_per_gene_plot
distances <- pval_dist / distplot_positive / distplot_negative
qc <- (nucleosome_plot / tss_plot) & ggtitle("")
atac_dimplot <- atac_dimplot + theme(legend.position = "none") +
  xlab("UMAP 1") + ylab("UMAP 2")

qc_violin <- qc_dist[[2]] / qc_dist[[1]] & theme_bw() & theme(legend.position = "none",
                                                              axis.text.x = element_blank(),
                                                              axis.ticks = element_blank())

top.panel <- (qc_violin | qc | atac_dimplot) + plot_layout(widths = c(1, 1, 2.5)) &
  theme(text = element_text(size = 12), axis.text = element_text(size = 12))

panel1 <- (mp / wrap_plots(tf_chromvar) / tf_expression & xlab("")) & theme(text = element_text(size = 12))

ggsave(filename = "figures/figure2.png", plot = top.panel, width = 12, height = 5.5, units = "in", dpi = 400)
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
