# Afreq

*A lightweight command-line tool to extract allele frequencies and sequencing depths from a VCF file.*

Afreq extracts allele depths (AD) from a multi-sample VCF file using bcftools and computes:

- **reference allele frequencies**

- **total sequencing depths**

These tables can then be used as inputs for **Gfreq** to estimate genotype frequencies in mixtures.

Only **biallelic SNPs** are processed. Multiallelic sites are ignored with a warning. 


## Table of Contents
- [Installation](#installation)
- [Command-line usage](#command-line-usage)
- [Output files](#output-files)

## Installation
### Using the container (recommended)
See the container instructions in the [main README](../README.md#container-image).

Run help:
```bash
apptainer exec mip-poolfreq.sif Afreq --help
```

## Command-line usage

Afreq currently provides a single command:

`Afreq estimate`

It reads a VCF file containing allele depth information (AD field) and produces two tables:

- reference allele frequencies per SNP

- total sequencing depth per SNP

These outputs can be used directly as inputs for **Gfreq**.


### Arguments
```
--vcf PATH          Input VCF file
--out-prefix STR    Prefix for output files
-h, --help          Show help message
```

### Example
```
SIF=$(realpath mip-poolfreq.sif)
apptainer exec $SIF Afreq estimate \
  --vcf input.vcf.gz \
  --out-prefix mixture
```

## Output files

Afreq produces two TSV files:

- **`*_ref_allelic_freqs.tsv`**

Reference allele frequencies per SNP and per sample.

```text
CHROM  POS  Sample1  Sample2  Sample3
chr1   123  0.50     0.75     0.12
chr1   456  0.33     0.66     0.00
...
```


Values correspond to:
<p align="center">
<code>REF_depth / (REF_depth + ALT_depth)</code>
</p>

- **`*_total_depths.tsv`**

Total sequencing depth per SNP and per sample.
```text
CHROM  POS  Sample1  Sample2  Sample3
chr1   123  42       58       35
chr1   456  61       73       22
...
```

Values correspond to:
<p align="center">
<code>REF_depth + ALT_depth</code>
</p>
