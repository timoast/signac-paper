#! /bin/bash
declare -a ncell=("1000" "3000" "5000" "7000" "9000" "11000" "13000" "15000" "17000" "19000" "21000" "23000" "25000")
declare -a ncore=("1" "2" "4" "8")

[ -d data/pbmc_atac/benchmarks ] || mkdir data/pbmc_atac/benchmarks

# feature matrix
for i in "${ncell[@]}"; do
  # need to run cores separately since it doesn't seem to obey the set number of threads
  taskset --cpu-list 1 Rscript --vanilla code/downsampling_code/run_archr_peakmatrix.R \
      1 \
      archr_pbmc/$i \
      data/pbmc_atac/peaks.bed \
      3 \
      data/pbmc_atac/benchmarks/archr_featmat_runtime_${i}_1.txt \
      "hg19"
  
  taskset --cpu-list 1,2 Rscript --vanilla code/downsampling_code/run_archr_peakmatrix.R \
      2 \
      archr_pbmc/$i \
      data/pbmc_atac/peaks.bed \
      3 \
      data/pbmc_atac/benchmarks/archr_featmat_runtime_${i}_2.txt \
      "hg19"
    
  taskset --cpu-list 1,2,3,4 Rscript --vanilla code/downsampling_code/run_archr_peakmatrix.R \
      4 \
      archr_pbmc/$i \
      data/pbmc_atac/peaks.bed \
      3 \
      data/pbmc_atac/benchmarks/archr_featmat_runtime_${i}_4.txt \
      "hg19"

  taskset --cpu-list 1,2,3,4,5,6,7,8 Rscript --vanilla code/downsampling_code/run_archr_peakmatrix.R \
      8 \
      archr_pbmc/$i \
      data/pbmc_atac/peaks.bed \
      3 \
      data/pbmc_atac/benchmarks/archr_featmat_runtime_${i}_8.txt \
      "hg19"
done


# gene activity
for i in "${ncell[@]}"; do
  # need to run cores separately since it doesn't seem to obey the set number of threads
  taskset --cpu-list 1 Rscript --vanilla code/downsampling_code/run_archr_geneactivity.R \
      1 \
      archr_pbmc/$i \
      3 \
      data/pbmc_atac/benchmarks/archr_geneactivity_runtime_${i}_1.txt \
      "hg19"
      
  taskset --cpu-list 1,2 Rscript --vanilla code/downsampling_code/run_archr_geneactivity.R \
      2 \
      archr_pbmc/$i \
      3 \
      data/pbmc_atac/benchmarks/archr_geneactivity_runtime_${i}_2.txt \
      "hg19"
      
  taskset --cpu-list 1,2,3,4 Rscript --vanilla code/downsampling_code/run_archr_geneactivity.R \
      4 \
      archr_pbmc/$i \
      3 \
      data/pbmc_atac/benchmarks/archr_geneactivity_runtime_${i}_4.txt \
      "hg19"
      
  taskset --cpu-list 1,2,3,4,5,6,7,8 Rscript --vanilla code/downsampling_code/run_archr_geneactivity.R \
      8 \
      archr_pbmc/$i \
      3 \
      data/pbmc_atac/benchmarks/archr_geneactivity_runtime_${i}_8.txt \
      "hg19"
done

# lsi
for i in "${ncell[@]}"; do
  taskset --cpu-list 1 Rscript --vanilla code/downsampling_code/run_archr_lsi.R \
    archr_pbmc/$i \
    data/pbmc_atac/peaks.bed \
    3 \
    data/pbmc_atac/benchmarks/archr_lsi_runtime_${i}.txt \
    "hg19"
done

# estimated lsi
for i in "${ncell[@]}"; do
  taskset --cpu-list 1 Rscript --vanilla code/downsampling_code/run_archr_estimated_lsi.R \
    archr_pbmc/$i \
    data/pbmc_atac/peaks.bed \
    3 \
    data/pbmc_atac/benchmarks/archr_est_lsi_runtime_${i}.txt \
    "hg19"
done