options(repos = c("CRAN" = "https://cran.rstudio.com/"))
options(Ncpus = 4)

install.packages(
  pkgs = c("remotes", "BiocManager", "tidyr", "dplyr", "patchwork")
)
BiocManager::install()
BiocManager::install(pkgs = c("GenomeInfoDbData", "HSMMSingleCell", "GO.db"))
setRepositories(ind = 1:2)

remotes::install_github(repo = "satijalab/seurat", ref = "release/4.0.0", dependencies = TRUE)
remotes::install_github(repo = "timoast/signac", ref = "develop", dependencies = TRUE)
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
           "DropletUtils",
           "chromVAR",
           "HDF5Array",
           "DelayedMatrixStats",
           "batchelor",
           "scater"
           )
  )

remotes::install_github('satijalab/seurat-wrappers')
remotes::install_github('cole-trapnell-lab/leidenbase')
remotes::install_github('cole-trapnell-lab/monocle3')
