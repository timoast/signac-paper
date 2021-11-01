library(Seurat)
library(Signac)
library(ggplot2)
library(patchwork)
library(ROCR)

set.seed(1234)

p.val.threshold <- 0.01
fc.threshold <- 0.4

pbmc <- readRDS("objects/pbmc.rds")
idents.use <- c("NK", "CD4 TCM")

de.all <-  FindMarkers(
  object = pbmc,
  group.by = "celltype",
  ident.1 = "NK",
  ident.2 = "CD4 TCM",
  test.use = "LR",
  min.pct = 0,
  logfc.threshold = 0,
  latent.vars = "nCount_ATAC"
)

gt <- de.all[de.all$p_val_adj < p.val.threshold & abs(de.all$avg_log2FC) > fc.threshold, ]

sample_groups <- function(obj, ident.1, ident.2, n, p1, p2) {
  cells1 <- WhichCells(object = obj, idents = ident.1)
  cells2 <- WhichCells(object = obj, idents = ident.2)
  
  # shuffle order
  cells1 <- sample(cells1)
  cells2 <- sample(cells2)
  
  cells1.group1 <- head(cells1, n = p1)
  cells2.group1 <- head(cells2, n = n-p1)
  
  cells1.remain <- setdiff(cells1, cells1.group1)
  cells2.remain <- setdiff(cells2, cells2.group1)
  
  cells1.group2 <- head(cells1.remain, n = p2)
  cells2.group2 <- head(cells2.remain, n = n-p2)
  
  return(list(c(cells1.group1, cells2.group1), c(cells1.group2, cells2.group2)))
}

de.results <- list()
p.list <- list()
for (i in seq(0, 100, 10)) {
  message("Adding ", i, " cells to contrast group")
  cells.use <- sample_groups(
    obj = pbmc,
    ident.1 = "NK",
    ident.2 = "CD4 TCM",
    n = 100,
    p1 = i,
    p2 = 100
  )
  pbmc$group <- ifelse(test = colnames(pbmc) %in% cells.use[[1]], yes = "A",
                       no = ifelse(colnames(pbmc) %in% cells.use[[2]], "B", NA))
  p1 <-  CoveragePlot(
    object = pbmc,
    region = "chr12-52566466-52567398",
    group.by = "group",
    idents = c("A", "B"),
    extend.upstream = 100,
    extend.downstream = 100,
    ymax = 54,
    annotation = FALSE
  )
  de.sampled <-  FindMarkers(
    object = pbmc,
    group.by = "group",
    ident.1 = "A",
    ident.2 = "B",
    test.use = "LR",
    min.pct = 0,
    logfc.threshold = 0,
    latent.vars = "nCount_ATAC"
  )
  de.results[[as.character(i)]] <- de.sampled
  p.list[[as.character(i)]] <- p1
}

for (i in seq_along(p.list)) {
  
  res1 <- de.results[[i]]["chr12-52566466-52567398", ]
  
  p.list[[i]][[1]] <- p.list[[i]][[1]] + ggtitle(
    paste0(
      (11 - i) * 10, "% cells from different pop.",
      "\nAdjusted p-value = ", format(res1$p_val_adj, digits = 3), "\nlog2 fold change = ",
      format(res1$avg_log2FC, digits = 3))
  )
  p.list[[i]][[2]] <- p.list[[i]][[2]] + theme(axis.text.x = element_text(size = 5))
  
}

# compute AUPRC
true.da <- as.numeric(rownames(de.results[[1]]) %in% rownames(gt))

compute_roc <- function(true.da, pred.da) {
  pred <- prediction(predictions = pred.da, labels = true.da)
  perf <- performance(pred, "tpr", "fpr")
  fpr <- perf@x.values[[1]]
  tpr <- perf@y.values[[1]]
  auc.tmp <- performance(pred, "auc")
  auc <- as.numeric(auc.tmp@y.values)
  return(list(data.frame(fpr = fpr, tpr = tpr), auc))
}

roc.df <- data.frame()
roc.auc <- c()
for (j in seq_along(de.results)) {
  pred.da <- abs(de.results[[j]]$avg_log2FC)
  rc <- compute_roc(true.da = true.da, pred.da = pred.da)
  rc[[1]]$cell_percent = (11-j) * 10
  roc.df <- rbind(roc.df, rc[[1]])
  roc.auc <- c(roc.auc, rc[[2]])
}

roc.df <- roc.df[roc.df$cell_percent > 0, ]
roc.df$cell_percent <- factor(roc.df$cell_percent, levels = seq(100, 10, -10))

p <- ggplot(roc.df, aes(x = fpr, y = tpr, color = cell_percent)) +
  geom_line() +
  xlab("False positive rate") +
  ylab("True positive rate") +
  theme_bw() +
  ggtitle(paste0("AUROC = ", format(roc.auc[[1]], digits = 3))) +
  labs(color = "Percent cells from\ndifferent population")

p1 <- wrap_plots(c(p.list, list(p))) & theme(plot.title = element_text(size = 10))
ggsave(filename = "figures/de_peaks_downsample_fig.png", plot = p1, height = 10, width = 12, dpi = 500)
