#! /bin/bash

declare -a ncell=("50000" "100000" "200000" "300000" "400000" "500000" "600000" "700000")
declare -a ncore=("1" "2" "4" "8")

[ -d data/biccn/benchmarks ] || mkdir data/biccn/benchmarks

# feature matrix
for i in "${ncell[@]}"; do
  # need to run cores separately since it doesn't seem to obey the set number of threads
  taskset --cpu-list 1 Rscript --vanilla code/downsampling_code/run_archr_peakmatrix.R \
      1 \
      archr_biccn/$i \
      data/biccn/unified_peaks.bed \
      3 \
      data/biccn/benchmarks/archr_featmat_runtime_${i}_1.txt \
      "mm10"

  taskset --cpu-list 1,2 Rscript --vanilla code/downsampling_code/run_archr_peakmatrix.R \
      2 \
      archr_biccn/$i \
      data/biccn/unified_peaks.bed \
      3 \
      data/biccn/benchmarks/archr_featmat_runtime_${i}_2.txt \
      "mm10"

  taskset --cpu-list 1,2,3,4 Rscript --vanilla code/downsampling_code/run_archr_peakmatrix.R \
      4 \
      archr_biccn/$i \
      data/biccn/unified_peaks.bed \
      3 \
      data/biccn/benchmarks/archr_featmat_runtime_${i}_4.txt \
      "mm10"

  taskset --cpu-list 1,2,3,4,5,6,7,8 Rscript --vanilla code/downsampling_code/run_archr_peakmatrix.R \
      8 \
      archr_biccn/$i \
      data/biccn/unified_peaks.bed \
      3 \
      data/biccn/benchmarks/archr_featmat_runtime_${i}_8.txt \
      "mm10"
done

for i in "${ncell[@]}"; do
  # need to run cores separately since it doesn't seem to obey the set number of threads
  taskset --cpu-list 1 Rscript --vanilla code/downsampling_code/run_archr_geneactivity.R \
      1 \
      archr_biccn/$i \
      3 \
      data/biccn/benchmarks/archr_geneactivity_runtime_${i}_1.txt \
      "mm10"

  taskset --cpu-list 1,2 Rscript --vanilla code/downsampling_code/run_archr_geneactivity.R \
      2 \
      archr_biccn/$i \
      3 \
      data/biccn/benchmarks/archr_geneactivity_runtime_${i}_2.txt \
      "mm10"

  taskset --cpu-list 1,2,3,4 Rscript --vanilla code/downsampling_code/run_archr_geneactivity.R \
      4 \
      archr_biccn/$i \
      3 \
      data/biccn/benchmarks/archr_geneactivity_runtime_${i}_4.txt \
      "mm10"

  taskset --cpu-list 1,2,3,4,5,6,7,8 Rscript --vanilla code/downsampling_code/run_archr_geneactivity.R \
      8 \
      archr_biccn/$i \
      3 \
      data/biccn/benchmarks/archr_geneactivity_runtime_${i}_8.txt \
      "mm10"
done

# lsi
for i in "${ncell[@]}"; do
    taskset --cpu-list 1 Rscript --vanilla code/downsampling_code/run_archr_lsi.R \
    archr_biccn/$i \
    data/biccn/unified_peaks.bed \
    1 \
    data/biccn/benchmarks/archr_lsi_runtime_${i}.txt \
    "mm10"
done
