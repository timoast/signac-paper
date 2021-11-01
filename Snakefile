datasets = ["pbmc", "biccn", "mito"]

rule all:
    input:
        "figures/figure2.png",
        "figures/dimreduc_pbmc.png",
        "figures/figure4.png",
        "figures/figure5a.png",
        "figures/figure5b.png",
        "figures/peakcalls.png",
        "figures/cluster_param_sweep.png",
        "figures/multimodal_label_transfer.png"

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
    threads: 1
    shell:
        """
        wget -i {input} -P data/{wildcards.dset}
        touch data/{wildcards.dset}/done.txt
        """

rule download_pbmc_atac:
    input:
        "datasets/pbmc_atac.txt"
    output:
        "data/pbmc_atac/done.txt"
    message: "Download PBMC scATAC fragment files"
    threads: 1
    shell:
        """
        wget -i {input} -P data/pbmc_atac
        touch data/pbmc_atac/done.txt
        """

rule download_gtex:
    input:
        "datasets/gtex.txt"
    output:
        "data/gtex/GTEx_v8_finemapping_CAVIAR/CAVIAR_Results_v8_GTEx_LD_HighConfidentVariants.gz"
    message: "Download GTEx data"
    threads: 1
    shell:
        """
        wget -i {input} -P data/gtex
        cd data/gtex
        tar -xvf GTEx_v8_finemapping_CAVIAR.tar
        rm GTEx_v8_finemapping_CAVIAR.tar
        """

rule download_chen:
    input: "datasets/chen.txt"
    output: "data/chen/scATAC-benchmarking-master/Synthetic_Data/BoneMarrow_cov250/input/bonemarrow_cov250.rds"
    message: "Downloading Chen et al. dataset"
    threads: 1
    shell:
        """
        wget -i {input} -P data/chen
        cd data/chen
        unzip master.zip
        rm master.zip
        """

# ------- Data processing -------

rule combine_fragments:
    input:
        "data/pbmc_atac/done.txt"
    output:
        "data/pbmc_atac/fragments.bed.gz"
    threads: 4
    message: "Combine PBMC scATAC fragment files"
    shell:
        """
        cd data/pbmc_atac
        gzip -d *.tsv.gz
        awk 'BEGIN {{FS=OFS="\\t"}} {{print $1,$2,$3,"10kng_"$4,$5}}' atac_pbmc_10k_nextgem_fragments.tsv > 1.bed
        awk 'BEGIN {{FS=OFS="\\t"}} {{print $1,$2,$3,"10k_"$4,$5}}' atac_pbmc_10k_v1_fragments.tsv > 2.bed
        awk 'BEGIN {{FS=OFS="\\t"}} {{print $1,$2,$3,"5kng_"$4,$5}}' atac_pbmc_5k_nextgem_fragments.tsv > 3.bed
        awk 'BEGIN {{FS=OFS="\\t"}} {{print $1,$2,$3,"5k_"$4,$5}}' atac_pbmc_5k_v1_fragments.tsv > 4.bed
        cat *.bed > frags.bed
        sort -k1,1 -k2,2n frags.bed > fragments.bed
        bgzip -@ {threads} fragments.bed
        tabix -p bed fragments.bed.gz
        rm *.bed *.tsv
        """

rule process_pbmc_atac:
    input:
        "data/pbmc_atac/fragments.bed.gz"
    output:
        "objects/pbmc_atac.rds"
    message: "Processing PBMC scATAC-seq data"
    threads: 1
    shell: "Rscript code/process_pbmc_atac.R"

rule process:
    input:
        "data/{dset}/done.txt",
        "install.done"
    output:
        "objects/{dset}.rds"
    message: "Process data"
    shell: "Rscript code/process_{wildcards.dset}.R"

# ------- Cell Downsampling -------

rule downsample_cells:
    input:
        "objects/biccn.rds",
        "objects/pbmc_atac.rds"
    output:
        "/scratch/tim/biccn/downsampling/50000.rds",
        "/scratch/tim/pbmc_atac/downsampling/25000.rds"
    message: "Downsample fragment files"
    threads: 1
    shell: "Rscript code/downsampling_code/downsample.R"

rule create_archr_files:
    input:
        "/scratch/tim/biccn/downsampling/50000.rds",
        "/scratch/tim/pbmc_atac/downsampling/25000.rds"
    output:
        "archr_pbmc/1000/Save-ArchR-Project.rds",
        "archr_biccn/700000/Save-ArchR-Project.rds"
    message: "Creating ArchR projects"
    threads: 1
    shell: "code/downsampling_code/downsample_archr.R"

rule create_annotations:
    input: "install.done"
    output:
        "data/biccn/annotations.rds",
        "data/pbmc_atac/annotations.rds"
    message: "Extracting annotations"
    threads: 1
    shell: "Rscript code/downsampling_code/get_annotations.R"

rule benchmark_biccn:
    input:
        "/scratch/tim/biccn/downsampling/50000.rds",
        "data/biccn/annotations.rds"
    output:
        "data/biccn/downsampling/svd_50000.rds"
    message: "Running BICCN benchmarks"
    threads: 8
    shell: "bash code/biccn_downsampling/benchmark.sh"
    
rule benchmark_pbmc_atac:
    input:
        "/scratch/tim/pbmc_atac/downsampling/25000.rds",
        "data/pbmc_atac/annotations.rds"
    output:
        "data/pbmc_atac/downsampling/svd_25000.rds"
    message: "Running PBMC scATAC benchmarks"
    threads: 8
    shell: "bash code/pbmc_atac_downsampling/benchmark.sh"

rule archr_pbmc_atac:
    input: "archr_pbmc/1000/Save-ArchR-Project.rds"
    output: "data/pbmc_atac/archr_lsi_runtime_1000_1.txt"
    message: "Running PBMC scATAC ArchR benchmarks"
    threads: 8
    shell: "bash code/pbmc_atac_downsampling/benchmark_archr.sh"

rule archr_biccn:
    input: "archr_biccn/700000/Save-ArchR-Project.rds"
    output: "data/biccn/archr_lsi_runtime_700000_1.txt"
    message: "Running BICCN scATAC ArchR benchmarks"
    threads: 8
    shell: "bash code/biccn_downsampling/benchmark_archr.sh"

rule gather_benchmark_timings:
    input:
        "data/biccn/downsampling/svd_50000.rds",
        "data/pbmc_atac/downsampling/svd_25000.rds",
    output:
        "data/biccn/timings.tsv",
        "data/pbmc_atac/timings.tsv"
    message: "Collating benchmark data"
    shell:
        """
        Rscript code/biccn_downsampling/collate_timings.R
        Rscript code/pbmc_atac_downsampling/collate_timings.R
        """

# ------- Count Downsampling -------

rule tfidf_downsampling:
    input: "objects/pbmc.rds"
    output: "data/pbmc/downsamples/1.rds"
    message: "Downsampling count matrix"
    threads: 1
    shell: "Rscript code/pbmc_downsampling/run_pbmc_downsample.R"

rule downsample_lsi:
    input:
        "data/pbmc/downsamples/1.rds",
        "data/chen/scATAC-benchmarking-master/Synthetic_Data/BoneMarrow_cov250/input/bonemarrow_cov250.rds"
    output:
        "data/pbmc/downsamples/lsi_1_1.rds",
        "data/chen/embeddings/lsi_250_1.rds"
    message: "Running LSI on downsampled counts"
    threads: 1
    shell: "Rscript code/pbmc_downsampling/run_lsi.R"

rule downsample_cistopic:
    input:
        "data/pbmc/downsamples/1.rds",
        "data/chen/scATAC-benchmarking-master/Synthetic_Data/BoneMarrow_cov250/input/bonemarrow_cov250.rds"
    output:
        "data/pbmc/downsamples/cistopic_warp_1.rds",
        "data/pbmc/downsamples/cistopic_cgs_1.rds",
        "data/chen/embeddings/cistopic_warp_250.rds",
        "data/chen/embeddings/cistopic_cgs_250.rds"
    message: "Running cisTopic on downsampled counts"
    threads: 1
    shell: "Rscript code/pbmc_downsampling/run_cistopic.R"
    
rule downsample_snapatac:
    input:
        "data/pbmc/downsamples/1.rds",
        "data/chen/scATAC-benchmarking-master/Synthetic_Data/BoneMarrow_cov250/input/bonemarrow_cov250.rds"
    output:
        "data/pbmc/downsamples/snapatac_1.rds",
        "data/chen/embeddings/snapatac_250.rds"
    message: "Running SnapATAC on downsampled counts"
    threads: 1
    shell: "Rscript code/pbmc_downsampling/run_snapatac.R"

rule downsample_scale:
    input:
        "data/pbmc/downsamples/1.rds",
        "data/chen/scATAC-benchmarking-master/Synthetic_Data/BoneMarrow_cov250/input/bonemarrow_cov250.rds"
    output:
        "data/pbmc/downsamples/scale_1/adata.h5ad",
        "data/chen/embeddings/scale_250/adata.h5ad"
    message: "Running SCALE on downsampled counts"
    threads: 6
    shell: "Rscript code/pbmc_downsampling/run_scale.R"

rule evaluate_dimreduc:
    input:
        "data/pbmc/downsamples/lsi_1_1.rds",
        "data/pbmc/downsamples/cistopic_warp_1.rds",
        "data/pbmc/downsamples/snapatac_1.rds",
        "data/pbmc/downsamples/scale_1/adata.h5ad",
        "data/chen/embeddings/lsi_250_1.rds",
        "data/chen/embeddings/cistopic_warp_250.rds",
        "data/chen/embeddings/snapatac_250.rds",
        "data/chen/embeddings/scale_250/adata.h5ad"
    output: "figures/dimreduc_pbmc.png", "figures/dimreduc_chen.png"
    message: "Evaluating dimension reduction methods"
    threads: 1
    shell: "Rscript code/pbmc_downsampling/evaluate_dimreducs.R"
    
# ------- Total runtime -------

rule get_total_runtime:
    input: "objects/biccn.rds", "objects/pbmc_atac.rds"
    output:
        "data/biccn/archr_total_runtime.txt",
        "data/biccn/signac_total_runtime.txt",
        "data/pbmc_atac/archr_total_runtime.txt",
        "data/pbmc_atac/signac_total_runtime.txt"
    threads: 8
    message: "Computing total runtime for ArchR and Signac"
    shell:
        """
        Rscript code/create_biccn_signac.R
        Rscript code/create_biccn_archr.R
        Rscript code/create_pbmc_atac_signac.R
        Rscript code/create_pbmc_atac_archr.R
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

rule multimodal_transfer:
    input:
        "objects/pbmc.rds"
    output:
        "objects/multimodal_label_transfer.rds"
    message: "Running multimodal label transfer"
    threads: 1
    shell:
        """
        Rscript code/multimodal_label_transfer.R
        """

rule analyze_pbmc:
    input:
        "objects/pbmc.rds",
        "objects/pbmc_links.rds",
        "objects/multimodal_label_transfer.rds",
        "eqtl.done"
    output:
        "figures/tss_enrichment.rds",
        "figures/multimodal_label_transfer.png"
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

rule fig_s1:
    input:
        "objects/pbmc.rds"
    output:
        "figures/peakcalls.png"
    shell:
        """
        Rscript code/peak_calling.R
        """

rule fig_s2:
    input:
        "objects/pbmc.rds"
    output:
        "figures/cluster_param_sweep.png"
    shell:
        """
        Rscript code/clustering.R
        """
