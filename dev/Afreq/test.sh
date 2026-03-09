#!/usr/bin/env bash
set -euo pipefail


TEST_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="$TEST_DIR/data/test_cases.vcf"
EXPECTED="$TEST_DIR/expected"

TMP=$(mktemp -d)
trap "rm -rf $TMP" EXIT

OUT_PREFIX="$TMP/out"

echo "Running Afreq estimate"
if ! Afreq estimate \
    --vcf "$DATA" \
    --out-prefix "$OUT_PREFIX" \
    >"$TMP/stdout.txt" \
    2>"$TMP/stderr.txt"
then
    echo "Command failed:"
    cat "$TMP/stderr.txt"
    exit 1
fi

echo "Checking allele frequencies"
diff -u \
  "$EXPECTED/expected_ref_allelic_freqs.tsv" \
  "$OUT_PREFIX"_ref_allelic_freqs.tsv

echo "Checking depths"
diff -u \
  "$EXPECTED/expected_total_depths.tsv" \
  "$OUT_PREFIX"_total_depths.tsv

echo "Checking warnings"
diff -u \
  "$EXPECTED/expected_warnings.txt" \
  "$TMP/stderr.txt"

echo
echo "All Afreq tests passed."
