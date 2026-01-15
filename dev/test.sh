#!/bin/bash

cd "mip-poolfreq"

# test pipelines
Rscript dev/test_estimate_genotype_frequencies_pipeline.R
Rscript dev/test_frequency_accuracy_pipeline.R
Rscript dev/test_depth_per_marker_pipeline.R
Rscript dev/test_MAF_hist_pipeline.R


# test CLI
(
  cd Gfreq
  Rscript -e "devtools::install('.')"
)

input_dir="Gfreq/example_data/input_files"

mkdir -p dev/outputs/estimate_genotype_frequencies
Gfreq estimate_genotype_frequencies \
  -v ${input_dir}/comp_genotypes.vcf \
  -a ${input_dir}/mix_ref_all_freqs.tsv \
  -d ${input_dir}/mix_read_depths.tsv \
  -t 40 \
  -l ${input_dir}/comp_libnames_corresp.tsv \
  -o dev/outputs/estimate_genotype_frequencies/genotype_frequencies.tsv

mkdir -p dev/outputs/frequency_accuracy
Gfreq evaluate_frequency_accuracy \
  -v ${input_dir}/comp_genotypes.vcf \
  --allele_freqs_mix ${input_dir}/mix_ref_all_freqs.tsv \
  --depths_mix ${input_dir}/mix_read_depths.tsv \
  --exp_freqs_mix ${input_dir}/mix_exp_geno_freqs.tsv \
  -t 40 \
  --allele_freqs_comp ${input_dir}/comp_ref_all_freqs.tsv \
  --depths_comp ${input_dir}/comp_read_depths.tsv \
  --exp_freqs_comp ${input_dir}/comp_exp_geno_freqs.tsv \
  -l ${input_dir}/comp_libnames_corresp.tsv \
  -s 5 \
  -r 20 \
  -o dev/outputs/frequency_accuracy

mkdir -p dev/outputs/MAF_hist
Gfreq plot_MAF_hist \
  -v ${input_dir}/comp_genotypes.vcf \
  -l ${input_dir}/comp_libnames_corresp.tsv \
  -o dev/outputs/MAF_hist

mkdir -p dev/outputs/depth_per_marker
Gfreq plot_depth_per_marker \
  -d dev/depth_files.list \
  -l 50 \
  -o dev/outputs/depth_per_marker
