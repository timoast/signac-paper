#! /bin/bash

declare -a ncell=("50000" "100000" "200000" "300000" "400000" "500000" "600000" "700000")
declare -a ncore=("1" "2" "4" "8")

[ -d data/biccn/benchmarks ] || mkdir data/biccn/benchmarks
[ -d data/biccn/downsampling ] || mkdir data/biccn/downsampling

# run each step with different numbers of cores, profile max memory usage

# feature matrix
for i in "${ncell[@]}"; do
  for j in "${ncore[@]}"; do
    /usr/bin/time -o data/biccn/benchmarks/featmat_mem_${i}_${j}.txt \
      -v Rscript --vanilla code/downsampling_code/run_featurematrix.R \
      $j \
      /scratch/tim/biccn/downsampling/$i.rds \
      data/biccn/unified_peaks.bed \
      3 \
      data/biccn/benchmarks/featmat_runtime_${i}_${j}.txt \
      data/biccn/downsampling/counts_${i}.rds
  done
done

# nucleosome signal
for i in "${ncell[@]}"; do
  /usr/bin/time -o data/biccn/benchmarks/nucleosome_mem_${i}.txt \
    -v Rscript --vanilla code/downsampling_code/run_nucleosome.R \
    data/biccn/downsampling/counts_${i}.rds \
    /scratch/tim/biccn/downsampling/$i.rds \
    data/biccn/annotations.rds \
    3 \
    data/biccn/benchmarks/nucleosome_runtime_${i}.txt \
    data/biccn/downsampling/nucleosome_${i}.rds
done


# tss enrichment
for i in "${ncell[@]}"; do
  for j in "${ncore[@]}"; do
    /usr/bin/time -o data/biccn/benchmarks/tss_mem_${i}_${j}.txt \
      -v Rscript --vanilla code/downsampling_code/run_tss.R \
      data/biccn/downsampling/nucleosome_${i}.rds \
      3 \
      data/biccn/benchmarks/tss_runtime_${i}_${j}.txt \
      data/biccn/downsampling/tss_${i}.rds \
      $j
  done
done

# gene activity matrix
for i in "${ncell[@]}"; do
  for j in "${ncore[@]}"; do
    /usr/bin/time -o data/biccn/benchmarks/ga_mem_${i}_${j}.txt \
      -v Rscript --vanilla code/downsampling_code/run_gene_activity.R \
      data/biccn/downsampling/tss_${i}.rds \
      $j \
      3 \
      data/biccn/benchmarks/ga_runtime_${i}_${j}.txt
  done
done

# tf-idf
for i in "${ncell[@]}"; do
  /usr/bin/time -o data/biccn/benchmarks/tfidf_mem_${i}.txt \
    -v Rscript --vanilla code/downsampling_code/run_tfidf.R \
    data/biccn/downsampling/tss_${i}.rds \
    3 \
    data/biccn/benchmarks/tfidf_runtime_${i}.txt \
    data/biccn/downsampling/tfidf_${i}.rds
done

# svd
for i in "${ncell[@]}"; do
  /usr/bin/time -o data/biccn/benchmarks/svd_mem_${i}.txt \
    -v Rscript --vanilla code/downsampling_code/run_svd.R \
    data/biccn/downsampling/tfidf_${i}.rds \
    3 \
    data/biccn/benchmarks/svd_runtime_${i}.txt \
    data/biccn/downsampling/svd_${i}.rds
done
