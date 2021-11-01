options(repos = c("CRAN" = "https://cran.rstudio.com/"))
options(Ncpus = 4)

install.packages(
  pkgs = c("remotes", "BiocManager", "tidyr", "dplyr", "RANN", "cluster", "ROCR",
           "patchwork", "mclust", "paletteer", "ggthemes", "dplyr", "arrow")
)
BiocManager::install()
BiocManager::install(pkgs = c("GenomeInfoDbData", "HSMMSingleCell", "GO.db", "DelayedArray"))
setRepositories(ind = 1:2)

install.packages("Seurat", dependencies = TRUE)
install.packages("Signac", dependencies = TRUE)
remotes::install_github(repo = "jlmelville/uwot")
remotes::install_github(repo = "mojaveazure/seurat-disk")

BiocManager::install(
  pkgs = c("EnsDb.Mmusculus.v79",
           "BSgenome.Mmusculus.UCSC.mm10",
           "TFBSTools",
           "JASPAR2020",
           "EnsDb.Hsapiens.v86",
           "BSgenome.Hsapiens.UCSC.hg38",
           "EnsDb.Hsapiens.v75",
           "BSgenome.Hsapiens.UCSC.hg19",
           "DropletUtils",
           "chromVAR",
           "HDF5Array",
           "DelayedMatrixStats",
           "batchelor",
           "scater"
           )
  )

# snapatac
install.packages(c("doSNOW", "plot3D"))
devtools::install_github("r3fang/SnapATAC")

# archr
devtools::install_github("GreenleafLab/ArchR", ref="release_1.0.1", repos = BiocManager::repositories())

# cistopic
devtools::install_github("aertslab/RcisTarget")
devtools::install_github("aertslab/AUCell")
devtools::install_github("aertslab/cisTopic")
