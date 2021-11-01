library(Signac)
library(EnsDb.Mmusculus.v79)
library(EnsDb.Hsapiens.v75)
library(GenomeInfoDb)

# extract gene annotations from EnsDb
annotations.mm <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
annotations.hg <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)

# change to UCSC style since the data was mapped to hg19/mm10
seqlevelsStyle(annotations.hg) <- 'UCSC'
genome(annotations.hg) <- "hg19"

seqlevelsStyle(annotations.mm) <- 'UCSC'
genome(annotations.mm) <- "mm10"

# save
saveRDS(object = annotations.mm, file = "data/biccn/annotations.rds")
saveRDS(object = annotations.hg, file = "data/pbmc_atac/annotations.rds")
