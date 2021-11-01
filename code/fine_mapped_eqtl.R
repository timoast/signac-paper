library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(tidyr)

set.seed(1234)

# ----- Data loading -----

pbmc <- readRDS("objects/pbmc.rds")
lnk <- readRDS("objects/pbmc_links.rds")

DefaultAssay(pbmc) <- "ATAC"
Links(pbmc) <- lnk

eqtl <- read.table(
  file = "data/gtex/GTEx_v8_finemapping_CAVIAR/CAVIAR_Results_v8_GTEx_LD_HighConfidentVariants.gz",
  header = TRUE
)

# add gene name
gene.annotations <- Annotation(pbmc)
gene.annotations <- as.data.frame(gene.annotations, row.names = NULL)[, c("gene_name", "gene_id")]
gene.annotations <- unique(gene.annotations)

gene.conversion <- gene.annotations$gene_name
names(gene.conversion) <- gene.annotations$gene_id

eqtl$gene_id_mod <- sapply(strsplit(x = eqtl$GENE, split = ".", fixed = TRUE), `[[`, 1)
eqtl$gene_name <- gene.conversion[eqtl$gene_id_mod]

eqtl <- makeGRangesFromDataFrame(
  df = eqtl,
  keep.extra.columns = TRUE,
  seqnames.field = "CHROM",
  start.field = "POS",
  end.field = "POS"
)

seqlevelsStyle(eqtl) <- "UCSC"
genome(eqtl) <- "hg38"

# subset to whole blood
eqtl_blood <- eqtl[eqtl$TISSUE == "Whole_Blood"]

# remove genes that don't have gene name or were not included in PBMC dataset
eqtl_blood <- eqtl_blood[!is.na(eqtl_blood$gene_name)]
eqtl_blood <- eqtl_blood[eqtl_blood$gene_name %in% rownames(pbmc[["SCT"]])]

# get transcript coordinates
transcripts <- Signac:::CollapseToLongestTranscript(ranges = Annotation(pbmc))

# get TSS coordinate
tss <- resize(transcripts, width = 1, fix = "start")

# look up gene start coord for each eQTL gene
start_coords <- vector()
for (i in seq_len(length.out = length(eqtl_blood))) {
  start_coords[[i]] <- start(tss[tss$gene_name == eqtl_blood$gene_name[[i]]])
}

# add to eqtls
eqtl_blood$tss_position <- start_coords
eqtl_blood$distance_to_gene <- abs(start(eqtl_blood) - eqtl_blood$tss_position)

# remove eQTLs more than 500 kb from linked gene TSS
eqtl_blood <- eqtl_blood[eqtl_blood$distance_to_gene <= 500000]

# ----- Intersection with links -----

linked_peak_coords <- StringToGRanges(regions = lnk$peak)
linked_peak_coords$gene <- lnk$gene
linked_peak_coords$score <- lnk$score

check_overlap <- function(
  signac_links,
  validated_links,
  gene_column = "gene_name"
) {
  # find eQTLs overlapping a linked peak
  olap <- findOverlaps(query = validated_links, subject = signac_links)
  
  # find the peak overlapped
  olap_signac <- signac_links[subjectHits(olap)]
  
  # find the eQTL overlapped
  olap_validated <- validated_links[queryHits(olap)]
  
  # combine information
  olap_signac$eqtl_gene <- mcols(olap_validated)[[gene_column]]
  olap_signac$posterior_prob <- mcols(olap_validated)[["Probability"]]
  olap_signac$eqtl <- mcols(olap_validated)[["eQTL"]]
  
  olap_signac <- as.data.frame(olap_signac)
  olap_signac <- unite(olap_signac, col = "peak", c("seqnames", "start", "end"), remove = FALSE)
  
  # for any peaks that overlap multiple eQTL variants linked to the same gene,
  # keep only one eQTL (don't count twice)
  # same peak, different eqtl, same linked gene
  overlap_counts <- olap_signac %>%
    group_by(peak) %>%
    summarise(n_eqtl = length(unique(eqtl_gene)),
              count = sum(unique(eqtl_gene) %in% unique(gene)))
  
  return(list(sum(overlap_counts$count), sum(overlap_counts$n_eqtl)))
}

linked_peak_overlap_eqtl <- check_overlap(
  signac_links = linked_peak_coords,
  validated_links = eqtl_blood,
  gene_column = 'gene_name'
)

# 52.6%

# find number of peaks intersecting an eqtl

olap <- subsetByOverlaps(granges(pbmc), eqtl_blood)
length(olap) / nrow(pbmc) * 100

ol2 <- subsetByOverlaps(StringToGRanges(lnk$peak), eqtl_blood)
# 2054
# 1.329415 %

# 2.328742

# ----- Permute links -----

# what's the expected number of correctly linked genes
# choose same number of linked peaks for each gene from peaks within 500 kb at random
gene.coords <- Signac:::CollapseToLongestTranscript(ranges = Annotation(pbmc))
genes <- rownames(x = pbmc[["SCT"]])
gene.coords.use <- gene.coords[gene.coords$gene_name %in% genes, ]
distances <- Signac:::DistanceToTSS(peaks = granges(pbmc), genes = gene.coords.use, distance = 5e+5)

# for each gene, get the number of linked peaks

lnk.df <- as.data.frame(lnk)
npeaks <- lnk.df %>% 
  group_by(gene) %>% 
  summarize(count = n())

n_linked_peaks <- npeaks$count
names(n_linked_peaks) <- npeaks$gene

# for each gene select n peaks at random from the set within 500 kb
df <- data.frame()
linked_genes <- unique(linked_peak_coords$gene)
for (i in seq_along(along.with = linked_genes)) {
  gene.use <- linked_genes[[i]]
  # find all peaks within 500 kb
  peak.choice <- distances[, gene.use]
  peak.choice <- names(peak.choice[peak.choice > 0])
  
  # work out how many peaks to choose
  nchoose <- n_linked_peaks[[gene.use]]
  
  # sample n
  peaks.selected <- sample(x = peak.choice, size = nchoose, replace = FALSE)
  
  # create dataframe
  subdf <- data.frame(peak = peaks.selected, gene = gene.use)
  df <- rbind(df, subdf)
}

df <- separate(df, col = peak, into = c("chromosome", "start", "end"), remove = FALSE)
permuted_links <- makeGRangesFromDataFrame(df = df, keep.extra.columns = TRUE)

# ----- eQTL overlap with permuted links -----

linked_permuted_overlap_eqtl <- check_overlap(
  signac_links = permuted_links,
  validated_links = eqtl_blood,
  gene_column = 'gene_name'
)
linked_permuted_overlap_eqtl

# 13.4%
