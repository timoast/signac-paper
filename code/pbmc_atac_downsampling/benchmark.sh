#! /bin/bash

declare -a ncell=("1000" "3000" "5000" "7000" "9000" "11000" "13000" "15000" "17000" "19000" "21000" "23000" "25000")
declare -a ncore=("1" "2" "4" "8")

[ -d data/pbmc_atac/benchmarks ] || mkdir data/pbmc_atac/benchmarks
[ -d data/pbmc_atac/downsampling ] || mkdir data/pbmc_atac/downsampling

# run each step with different numbers of cores, profile max memory usage

# feature matrix
for i in "${ncell[@]}"; do
  for j in "${ncore[@]}"; do
    /usr/bin/time -o data/pbmc_atac/benchmarks/featmat_mem_${i}_${j}.txt \
      -v Rscript --vanilla code/downsampling_code/run_featurematrix.R \
      $j \
      /scratch/tim/pbmc_atac/downsampling/$i.rds \
      data/pbmc_atac/peaks.bed \
      3 \
      data/pbmc_atac/benchmarks/featmat_runtime_${i}_${j}.txt \
      data/pbmc_atac/downsampling/counts_${i}.rds
  done
done

# nucleosome signal
for i in "${ncell[@]}"; do
  /usr/bin/time -o data/pbmc_atac/benchmarks/nucleosome_mem_${i}.txt \
    -v Rscript --vanilla code/downsampling_code/run_nucleosome.R \
    data/pbmc_atac/downsampling/counts_${i}.rds \
    /scratch/tim/pbmc_atac/downsampling/$i.rds \
    data/pbmc_atac/annotations.rds \
    3 \
    data/pbmc_atac/benchmarks/nucleosome_runtime_${i}.txt \
    data/pbmc_atac/downsampling/nucleosome_${i}.rds
done

# tss enrichment
for i in "${ncell[@]}"; do
  for j in "${ncore[@]}"; do
    /usr/bin/time -o data/biccn/benchmarks/tss_mem_${i}_${j}.txt \
      -v Rscript --vanilla code/downsampling_code/run_tss.R \
      data/pbmc_atac/downsampling/nucleosome_${i}.rds \
      3 \
      data/pbmc_atac/benchmarks/tss_runtime_${i}_${j}.txt \
      data/pbmc_atac/downsampling/tss_${i}.rds \
      $j
  done
done

# gene activity matrix
for i in "${ncell[@]}"; do
  for j in "${ncore[@]}"; do
    /usr/bin/time -o data/pbmc_atac/benchmarks/ga_mem_${i}_${j}.txt \
      -v Rscript --vanilla code/downsampling_code/run_gene_activity.R \
      data/pbmc_atac/downsampling/tss_${i}.rds \
      $j \
      3 \
      data/pbmc_atac/benchmarks/ga_runtime_${i}_${j}.txt
  done
done

# tf-idf
for i in "${ncell[@]}"; do
  /usr/bin/time -o data/pbmc_atac/benchmarks/tfidf_mem_${i}.txt \
    -v Rscript --vanilla code/downsampling_code/run_tfidf.R \
    data/pbmc_atac/downsampling/tss_${i}.rds \
    3 \
    data/pbmc_atac/benchmarks/tfidf_runtime_${i}.txt \
    data/pbmc_atac/downsampling/tfidf_${i}.rds
done

# svd
for i in "${ncell[@]}"; do
  /usr/bin/time -o data/pbmc_atac/benchmarks/svd_mem_${i}.txt \
    -v Rscript --vanilla code/downsampling_code/run_svd.R \
    data/pbmc_atac/downsampling/tfidf_${i}.rds \
    3 \
    data/pbmc_atac/benchmarks/svd_runtime_${i}.txt \
    data/pbmc_atac/downsampling/svd_${i}.rds
done
