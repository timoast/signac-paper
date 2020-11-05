datasets = ["pbmc", "biccn", "mito"]

rule all:
    input:
        "figures/figure2.png",
        "figures/figure3.png",
        "figures/figure4.png",
        "figures/figure5a.png",
        "figures/figure5b.png"

# ------- Install -------

rule install:
    output: "install.done"
    threads: 4
    shell:
        """
        Rscript install_r_packages.R
        touch install.done
        """

# ------- Data download -------

rule download:
    input:
        "datasets/{dset}.txt"
    output:
        "data/{dset}/done.txt"
    message: "Download datasets"
    shell:
        """
        wget -i {input} -P data/{wildcards.dset}
        touch data/{wildcards.dset}/done.txt
        """

rule download_gtex:
    input:
        "datasets/gtex.txt"
    output:
        "data/gtex/GTEx_v8_finemapping_CAVIAR/CAVIAR_Results_v8_GTEx_LD_HighConfidentVariants.gz"
    shell:
        """
        wget -i {input} -P data/gtex
        cd data/gtex
        tar -xvf GTEx_v8_finemapping_CAVIAR.tar
        rm GTEx_v8_finemapping_CAVIAR.tar
        """

# ------- Data processing -------

rule process:
    input:
        "data/{dset}/done.txt",
        "install.done"
    output:
        "objects/{dset}.rds"
    message: "Process data"
    shell: "Rscript code/process_{wildcards.dset}.R"

rule gather_benchmark_timings:
    input:
        "data/biccn/downsampling/svd_50000.rds"
    output:
        "data/biccn/timings.tsv"
    message: "Collating benchmark data"
    shell: "Rscript code/biccn_downsampling/collate_timings.R"

# ------- Downsampling -------

rule downsample:
    input:
        "objects/biccn.rds"
    output:
        "data/biccn/downsampling/50000.rds"
    message: "Downsample BICCN fragment file"
    threads: 1
    shell: "Rscript code/biccn_downsampling/downsample_biccn.R"

rule create_annotations:
    input: "install.done"
    output: "data/biccn/annotations.rds"
    message: "Extracting annotations"
    threads: 1
    shell: "Rscript code/biccn_downsampling/get_annotations.R"

rule benchmark:
    input:
        "data/biccn/downsampling/50000.rds",
        "data/biccn/annotations.rds"
    output:
        "data/biccn/downsampling/svd_50000.rds"
    message: "Running benchmarks"
    threads: 8
    shell: "bash code/biccn_downsampling/benchmark.sh"

rule tfidf_downsampling:
    input:
        "objects/pbmc.rds"
    output:
        "figures/dimplot_downsample_pbmc.rds"
    shell:
        """
        Rscript code/pbmc_downsampling/run_pbmc_downsample.R
        """

# ------- Links -------

rule links:
    input:
        "objects/pbmc.rds"
    output:
        "objects/pbmc_links.rds"
    threads: 8
    shell:
        """
        Rscript code/link_peaks.R
        """

# ------- Analysis -------

rule analyze_eqtl:
    input:
        "objects/pbmc.rds",
        "objects/pbmc_links.rds",
        "data/gtex/GTEx_v8_finemapping_CAVIAR/CAVIAR_Results_v8_GTEx_LD_HighConfidentVariants.gz"
    output:
        "eqtl.done"
    shell:
        """
        Rscript code/analyze_pbmc.R
        touch eqtl.done
        """

rule analyze_pbmc:
    input:
        "objects/pbmc.rds",
        "objects/pbmc_links.rds",
        "eqtl.done"
    output:
        "figures/tss_enrichment.rds"
    shell:
        """
        Rscript code/analyze_pbmc.R
        """

# ------- Figures -------

rule fig2:
    input:
        "figures/tss_enrichment.rds"
    output:
        "figures/figure2.png"
    shell:
        """
        Rscript code/figure2.R
        """

rule fig3:
    input:
        "figures/dimplot_downsample_pbmc.rds"
    output:
        "figures/figure3.png"
    shell:
        """
        Rscript code/figure3.R
        """

rule fig4:
    input:
        "objects/mito.rds"
    output:
        "figures/figure4.png"
    shell:
        """
        Rscript code/figure4.R
        """

rule fig5:
    input:
        "objects/biccn.rds",
        "data/biccn/timings.tsv"
    output:
        "figures/figure5a.png",
        "figures/figure5b.png"
    shell:
        """
        Rscript code/figure5.R
        """
