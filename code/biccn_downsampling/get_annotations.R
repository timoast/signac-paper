library(Signac)
library(EnsDb.Mmusculus.v79)
library(GenomeInfoDb)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

# save
saveRDS(object = annotations, file = "data/biccn/annotations.rds")