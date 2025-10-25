#!/bin/bash

# Check parameters
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <VCF_PATH> <SUFFIX>"
    exit 1
fi

VCF_PATH="$1"
SUFFIX="$2"

btquery_output="$SUFFIX".btquery
output_freqs=ref_allelic_freqs_"$SUFFIX".tsv
output_depths=total_depths_"$SUFFIX".tsv

# Extract allele depths (AD) with bcftools
bcftools query -f '%CHROM\t%POS[\t%AD]\n' "$VCF_PATH" > "$btquery_output"

# Retrieve the header (sample names)
grep '#' "$VCF_PATH" | tail -1 | sed 's/ID.*FORMAT\t//' | sed 's/^#//' > "$output_freqs"
grep '#' "$VCF_PATH" | tail -1 | sed 's/ID.*FORMAT\t//' | sed 's/^#//' > "$output_depths"

# Compute allele frequencies
while IFS=$'\t' read -r chr pos rest; do
  freqs=""
  depths=""

  IFS=$'\t' read -ra samples <<< "$rest"

  for sample in "${samples[@]}"; do
    IFS=',' read -ra al_depths <<< "$sample"

    if (( ${#al_depths[@]} != 2 )); then
      [[ -n "$freqs" ]] && freqs="$freqs\terror" || freqs="error"
      [[ -n "$depths" ]] && depths="$depths\terror" || depths="error"
      echo "Warning: more than 2 alleles for $chr $pos"
    else
      sum=$((al_depths[0] + al_depths[1]))

      if (( sum == 0 )); then
        [[ -n "$freqs" ]] && freqs="$freqs\tNA" || freqs="NA"
        [[ -n "$depths" ]] && depths="$depths\tNA" || depths="NA"
      else
        freq=$(echo "${al_depths[0]} / $sum" | bc -l | awk '{printf "%.2f", $0}')
        [[ -n "$freqs" ]] && freqs="$freqs\t$freq" || freqs="$freq"
        [[ -n "$depths" ]] && depths="$depths\t$sum" || depths="$sum"
      fi
    fi
  done

  echo -e "$chr\t$pos\t$freqs" >> "$output_freqs"
  echo -e "$chr\t$pos\t$depths" >> "$output_depths"
done < "$btquery_output"
