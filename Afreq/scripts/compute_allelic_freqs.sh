#!/usr/bin/env bash
set -euo pipefail

usage() {
    local prog="${1:-compute_allelic_freqs.sh}"

    cat <<EOF
Usage:
  ${prog} --vcf <file.vcf[.gz]> --out-prefix <prefix>

Description:
  Extract allele depths (AD) from a VCF using bcftools query and compute:
    - reference allele frequencies
    - total depths

  Only biallelic sites are processed.
  Multiallelic sites are ignored with a warning.

Outputs:
  <prefix>_ref_allelic_freqs.tsv
  <prefix>_total_depths.tsv

Options:
  --vcf PATH          Input VCF file
  --out-prefix STR    Prefix for output files
  -h, --help          Show this help message
EOF
}

VCF_PATH=""
OUT_PREFIX=""
USAGE_NAME="compute_allelic_freqs.sh"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --vcf)
            VCF_PATH="${2:-}"
            shift 2
            ;;
        --out-prefix)
            OUT_PREFIX="${2:-}"
            shift 2
            ;;
        --usage-name)
            USAGE_NAME="${2:-}"
            shift 2
            ;;
        -h|--help)
            usage "$USAGE_NAME"
            exit 0
            ;;
        *)
            echo "Error: unknown argument '$1'" >&2
            usage "$USAGE_NAME" >&2
            exit 1
            ;;
    esac
done

if [[ -z "$VCF_PATH" || -z "$OUT_PREFIX" ]]; then
    echo "Error: --vcf and --out-prefix are required." >&2
    usage "$USAGE_NAME" >&2
    exit 1
fi

if [[ ! -f "$VCF_PATH" ]]; then
    echo "Error: input VCF not found: $VCF_PATH" >&2
    exit 1
fi

if ! command -v bcftools >/dev/null 2>&1; then
    echo "Error: bcftools is not available in PATH." >&2
    exit 1
fi

output_freqs="${OUT_PREFIX}_ref_allelic_freqs.tsv"
output_depths="${OUT_PREFIX}_total_depths.tsv"

sample_header="$(
    bcftools query -l "$VCF_PATH" | paste -sd $'\t' -
)"

if [[ -z "$sample_header" ]]; then
    echo "Error: no sample names found in VCF: $VCF_PATH" >&2
    exit 1
fi

printf 'CHROM\tPOS\t%s\n' "$sample_header" > "$output_freqs"
printf 'CHROM\tPOS\t%s\n' "$sample_header" > "$output_depths"

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' "$VCF_PATH" | \
while IFS=$'\t' read -r chr pos ref alt rest; do
    if [[ "$alt" == *,* ]]; then
        echo "Warning: ${chr} ${pos}: ignored because it does not look like a biallelic site" >&2
        continue
    fi

    freqs=()
    depths=()

    IFS=$'\t' read -r -a samples <<< "$rest"

    for sample in "${samples[@]}"; do
        if [[ -z "$sample" || "$sample" == "." || "$sample" == ".,." ]]; then
            freqs+=("NA")
            depths+=("NA")
            continue
        fi

        IFS=',' read -r -a al_depths <<< "$sample"

        ref_depth="${al_depths[0]}"
        alt_depth="${al_depths[1]}"

        sum=$((ref_depth + alt_depth))

        if (( sum == 0 )); then
            freqs+=("NA")
            depths+=("NA")
        else
            freq="$(awk -v ref="$ref_depth" -v sum="$sum" 'BEGIN { printf "%.2f", ref / sum }')"
            freqs+=("$freq")
            depths+=("$sum")
        fi
    done

    printf '%s\t%s\t%s\n' "$chr" "$pos" "$(IFS=$'\t'; echo "${freqs[*]}")" >> "$output_freqs"
    printf '%s\t%s\t%s\n' "$chr" "$pos" "$(IFS=$'\t'; echo "${depths[*]}")" >> "$output_depths"
done
