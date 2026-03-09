#!/usr/bin/env bash
set -euo pipefail

TEST_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$TEST_DIR/../.." && pwd)"

INPUT_DIR="$REPO_ROOT/Gfreq/example_data/input_files"
EXPECTED_DIR="$REPO_ROOT/Gfreq/example_data/output_files"

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

OUT_ESTIMATE="$TMP/estimate"
OUT_MAF="$TMP/maf"
OUT_DEPTH="$TMP/depth"
OUT_EVAL="$TMP/eval"

mkdir -p "$OUT_ESTIMATE" "$OUT_EVAL" "$OUT_MAF" "$OUT_DEPTH"

step() {
  echo
  echo "## [$1/7] $2"
}

ok() {
  echo -e "\nOK: $1"
}

step 1 "Running Gfreq estimate"
Gfreq estimate \
  -v "$INPUT_DIR/comp_genotypes.vcf" \
  -a "$INPUT_DIR/mix_ref_all_freqs.tsv" \
  -d "$INPUT_DIR/mix_read_depths.tsv" \
  -t 40 \
  -l "$INPUT_DIR/comp_libnames_corresp.tsv" \
  -o "$OUT_ESTIMATE/genotype_frequencies.tsv"
ok "Gfreq estimate completed"

step 2 "Running Gfreq eval"
Gfreq eval \
  -v "$INPUT_DIR/comp_genotypes.vcf" \
  --allele_freqs_mix "$INPUT_DIR/mix_ref_all_freqs.tsv" \
  --depths_mix "$INPUT_DIR/mix_read_depths.tsv" \
  --exp_freqs_mix "$INPUT_DIR/mix_exp_geno_freqs.tsv" \
  -t 40 \
  --allele_freqs_comp "$INPUT_DIR/comp_ref_all_freqs.tsv" \
  --depths_comp "$INPUT_DIR/comp_read_depths.tsv" \
  --exp_freqs_comp "$INPUT_DIR/comp_exp_geno_freqs.tsv" \
  -l "$INPUT_DIR/comp_libnames_corresp.tsv" \
  -s 50 \
  -r 2 \
  --sampling_seed 123 \
  -o "$OUT_EVAL"
ok "Gfreq eval completed"

step 3 "Running Gfreq maf"
Gfreq maf \
  -v "$INPUT_DIR/comp_genotypes.vcf" \
  -l "$INPUT_DIR/comp_libnames_corresp.tsv" \
  -o "$OUT_MAF"
ok "Gfreq maf completed"

step 4 "Running Gfreq depth"
Gfreq depth \
  -d "$TEST_DIR/data/depth_files.list" \
  -l 50 \
  -o "$OUT_DEPTH"
ok "Gfreq depth completed"

step 5 "Checking estimate output"
diff -u \
  "$EXPECTED_DIR/estimate/genotype_frequencies.tsv" \
  "$OUT_ESTIMATE/genotype_frequencies.tsv"
ok "Estimate output matches expected file"

step 6 "Checking maf output"
diff -u \
  "$EXPECTED_DIR/maf/MAFs.tsv" \
  "$OUT_MAF/MAFs.tsv"
ok "maf output matches expected file"

step 7 "Checking depth output"
diff -u \
  "$EXPECTED_DIR/depth/mean_depth_per_marker.tsv" \
  "$OUT_DEPTH/mean_depth_per_marker.tsv"
ok "depth output matches expected file"

echo
echo "All Gfreq tests passed."
