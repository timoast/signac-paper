library(ggplot2)
library(patchwork)

dimplot_list <- readRDS("figures/dimplot_downsample_pbmc.rds")
knn_list <- readRDS("figures/knn_purity_pbmc.rds")

knn_df <- data.frame()
for (i in seq_along(knn_list)) {
  ds_level <- as.factor(unlist(strsplit(names(knn_list)[i], "_"))[1])
  method <- unlist(strsplit(names(knn_list)[i], "_"))[2]
  method <- ifelse(test = method == 1, "log(TF*IDF)", "TF*log(IDF)")
  frc <- knn_list[[i]]
  dft <- data.frame("downsample" = ds_level, "method" = method, "knn" = frc) 
  knn_df <- rbind(knn_df, dft)
}

knn_box <- ggplot(
  data = knn_df,
  mapping = aes(x = as.factor(downsample), y = knn, fill = method)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() + 
  scale_x_discrete(limits = levels(knn_df$downsample)) +
  ylab("KNN purity") +
  xlab("Fraction counts retained")

plots_use <- c("1_1", "0.6_1", "0.2_1", "1_2", "0.6_2", "0.2_2")

p <- wrap_plots(
  dimplot_list[plots_use],
  nrow = 2,
  guides = "collect"
) & xlab("UMAP1") & ylab("UMAP2")

fig2 <- p / knn_box + plot_layout(heights = c(3, 1))

ggsave(filename = "figures/figure3.png", plot = fig2, height = 10, width = 15, dpi = 300)
