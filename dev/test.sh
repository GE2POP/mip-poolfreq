#!/bin/bash
# Note: works in WSL

cd "mip-poolfreq"

# test build
(
  cd Gfreq
  rm -f Gfreq_*.tar.gz
  Rscript -e "devtools::document()"
  R CMD build .
  R CMD check --no-manual Gfreq_*.tar.gz
  rm Gfreq_*.tar.gz
)

# test pipelines
Rscript dev/test_estimate_pipeline.R
Rscript dev/test_eval_pipeline.R
Rscript dev/test_depth_pipeline.R
Rscript dev/test_maf_pipeline.R


# test CLI
(
  cd Gfreq
  Rscript -e "devtools::install('.')"
)

input_dir="Gfreq/example_data/input_files"

mkdir -p dev/outputs/estimate
Gfreq estimate \
  -v ${input_dir}/comp_genotypes.vcf \
  -a ${input_dir}/mix_ref_all_freqs.tsv \
  -d ${input_dir}/mix_read_depths.tsv \
  -t 40 \
  -l ${input_dir}/comp_libnames_corresp.tsv \
  -o dev/outputs/estimate/genotype_frequencies.tsv

mkdir -p dev/outputs/eval
Gfreq eval \
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
  --sampling_seed 123 \
  -o dev/outputs/eval

mkdir -p dev/outputs/maf
Gfreq maf \
  -v ${input_dir}/comp_genotypes.vcf \
  -l ${input_dir}/comp_libnames_corresp.tsv \
  -o dev/outputs/maf

mkdir -p dev/outputs/depth
Gfreq depth \
  -d dev/depth_files.list \
  -l 50 \
  -o dev/outputs/depth

# test Docker image
version=0.1.0
docker build -t gfreq:${version} .
docker run --rm gfreq:${version} --help

input_dir="Gfreq/example_data/input_files"

mkdir -p dev/outputs/estimate
docker run --rm \
  -v "$(pwd)/${input_dir}:/in:ro" \
  -v "$(pwd)/dev/outputs/estimate:/out" \
  gfreq:${version} estimate \
  -v /in/comp_genotypes.vcf \
  -a /in/mix_ref_all_freqs.tsv \
  -d /in/mix_read_depths.tsv \
  -t 40 \
  -l /in/comp_libnames_corresp.tsv \
  -o /out/genotype_frequencies.tsv

mkdir -p dev/outputs/eval
docker run --rm \
  -v "$(pwd)/${input_dir}:/in:ro" \
  -v "$(pwd)/dev/outputs/eval:/out" \
  gfreq:${version} eval \
  -v /in/comp_genotypes.vcf \
  --allele_freqs_mix /in/mix_ref_all_freqs.tsv \
  --depths_mix /in/mix_read_depths.tsv \
  --exp_freqs_mix /in/mix_exp_geno_freqs.tsv \
  -t 40 \
  --allele_freqs_comp /in/comp_ref_all_freqs.tsv \
  --depths_comp /in/comp_read_depths.tsv \
  --exp_freqs_comp /in/comp_exp_geno_freqs.tsv \
  -l /in/comp_libnames_corresp.tsv \
  -s 5 \
  -r 20 \
  --sampling_seed 123 \
  -o /out

mkdir -p dev/outputs/maf
docker run --rm \
  -v "$(pwd)/${input_dir}:/in:ro" \
  -v "$(pwd)/dev/outputs/maf:/out" \
  gfreq:${version} maf \
  -v /in/comp_genotypes.vcf \
  -l /in/comp_libnames_corresp.tsv \
  -o /out

mkdir -p dev/outputs/depth
docker run --rm \
  -v "$(pwd)/dev:/in:ro" \
  -v "$(pwd)/dev/outputs/depth:/out" \
  gfreq:${version} depth \
  -d /in/depth_files.list \
  -l 50 \
  -o /out
