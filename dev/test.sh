#!/usr/bin/env bash
set -euo pipefail


#### Test Gfreq R package build
(
  cd Gfreq
  rm -f Gfreq_*.tar.gz
  Rscript -e "devtools::document()"
  R CMD build .
  R CMD check --no-manual Gfreq_*.tar.gz
  rm Gfreq_*.tar.gz
)

#### Tests in Docker
if pwd -W >/dev/null 2>&1; then
  HOST_PWD="$(pwd -W)"
  DOCKER_RUN_PREFIX=(env MSYS_NO_PATHCONV=1 docker run --rm -v "${HOST_PWD}:/work" -w /work)
else
  HOST_PWD="$(pwd)"
  DOCKER_RUN_PREFIX=(docker run --rm -v "${HOST_PWD}:/work" -w /work)
fi

### Build Docker image
echo -e "\n### Building Docker image"
docker build -t mip-poolfreq .

echo -e "\n### Checking container help"
docker run --rm mip-poolfreq --help
docker run --rm mip-poolfreq Afreq --help
docker run --rm mip-poolfreq Afreq estimate --help
docker run --rm mip-poolfreq Gfreq --help
docker run --rm mip-poolfreq Gfreq estimate --help

### Test Afreq in Docker image
echo -e "\n### Testing Afreq in Docker image"
"${DOCKER_RUN_PREFIX[@]}" mip-poolfreq bash ./dev/Afreq/test.sh

### Test Gfreq in Docker image
echo -e "\n### Testing Gfreq in Docker image"
"${DOCKER_RUN_PREFIX[@]}" mip-poolfreq bash ./dev/Gfreq/test.sh

#### Tests in Apptainer
REPO_ROOT=$(pwd)

### Build Apptainer image
docker save mip-poolfreq:latest -o mip-poolfreq.tar
apptainer build --force mip-poolfreq.sif docker-archive:mip-poolfreq.tar
SIF=$(realpath mip-poolfreq.sif)

### Test Afreq in Apptainer image
echo -e "\n### Testing Afreq in Apptainer image"
apptainer exec --bind "${REPO_ROOT}:${REPO_ROOT}" "$SIF" \
  bash "${REPO_ROOT}/dev/Afreq/test.sh"

### Test Gfreq in Apptainer image
echo -e "\n### Testing Gfreq in Apptainer image"
apptainer exec --bind "${REPO_ROOT}:${REPO_ROOT}" "$SIF" \
  bash "${REPO_ROOT}/dev/Gfreq/test.sh"

exit 0




#### Gfreq extra hand testing
# Note: works in WSL

# test pipelines
Rscript dev/Gfreq/test_estimate_pipeline.R
Rscript dev/Gfreq/test_eval_pipeline.R
Rscript dev/Gfreq/test_depth_pipeline.R
Rscript dev/Gfreq/test_maf_pipeline.R


# test CLI
(
  cd Gfreq
  Rscript -e "devtools::install('.')"
)

input_dir="Gfreq/example_data/input_files"

mkdir -p dev/Gfreq/outputs/estimate
Gfreq estimate \
  -v ${input_dir}/comp_genotypes.vcf \
  -a ${input_dir}/mix_ref_all_freqs.tsv \
  -d ${input_dir}/mix_read_depths.tsv \
  -t 40 \
  -l ${input_dir}/comp_libnames_corresp.tsv \
  -o dev/Gfreq/outputs/estimate/genotype_frequencies.tsv

mkdir -p dev/Gfreq/outputs/eval
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
  -o dev/Gfreq/outputs/eval

mkdir -p dev/Gfreq/outputs/maf
Gfreq maf \
  -v ${input_dir}/comp_genotypes.vcf \
  -l ${input_dir}/comp_libnames_corresp.tsv \
  -o dev/Gfreq/outputs/maf

mkdir -p dev/Gfreq/outputs/depth
Gfreq depth \
  -d dev/Gfreq/data/depth_files.list \
  -l 50 \
  -o dev/Gfreq/outputs/depth
