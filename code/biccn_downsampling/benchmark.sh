#! /bin/bash

declare -a ncell=("50000" "100000" "200000" "300000" "400000" "500000" "600000" "700000")
declare -a ncore=("1" "2" "4" "8")

[ -d data/biccn/benchmarks ] || mkdir data/biccn/benchmarks

# run each step with different numbers of cores, profile max memory usage

# feature matrix
for i in "${ncell[@]}"; do
  for j in "${ncore[@]}"; do
    /usr/bin/time -o data/biccn/benchmarks/featmat_mem_${i}_${j}.txt \
      -v Rscript --vanilla code/biccn_downsampling/run_featurematrix.R \
      $j \
      data/biccn/downsampling/$i.rds \
      data/biccn/unified_peaks.bed \
      3 \
      data/biccn/benchmarks/featmat_runtime_${i}_${j}.txt \
      data/biccn/downsampling/counts_${i}.rds
  done
done

# nucleosome signal
for i in "${ncell[@]}"; do
  /usr/bin/time -o data/biccn/benchmarks/nucleosome_mem_${i}.txt \
    -v Rscript --vanilla code/biccn_downsampling/run_nucleosome.R \
    data/biccn/downsampling/counts_${i}.rds \
    data/biccn/downsampling/$i.rds \
    data/biccn/annotations.rds \
    3 \
    data/biccn/benchmarks/nucleosome_runtime_${i}.txt \
    data/biccn/downsampling/nucleosome_${i}.rds
done

# tss enrichment
for i in "${ncell[@]}"; do
  /usr/bin/time -o data/biccn/benchmarks/tss_mem_${i}.txt \
    -v Rscript --vanilla code/biccn_downsampling/run_tss.R \
    data/biccn/downsampling/nucleosome_${i}.rds \
    3 \
    data/biccn/benchmarks/tss_runtime_${i}.txt \
    data/biccn/downsampling/tss_${i}.rds
done

# tf-idf
for i in "${ncell[@]}"; do
  /usr/bin/time -o data/biccn/benchmarks/tfidf_mem_${i}.txt \
    -v Rscript --vanilla code/biccn_downsampling/run_tfidf.R \
    data/biccn/downsampling/tss_${i}.rds \
    3 \
    data/biccn/benchmarks/tfidf_runtime_${i}.txt \
    data/biccn/downsampling/tfidf_${i}.rds
done

# svd
for i in "${ncell[@]}"; do
  /usr/bin/time -o data/biccn/benchmarks/svd_mem_${i}.txt \
    -v Rscript --vanilla code/biccn_downsampling/run_svd.R \
    data/biccn/downsampling/tfidf_${i}.rds \
    3 \
    data/biccn/benchmarks/svd_runtime_${i}.txt \
    data/biccn/downsampling/svd_${i}.rds
done
