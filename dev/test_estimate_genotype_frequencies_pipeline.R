#!/usr/bin/env Rscript

## ---- Load packages in dev mode ----------------------------------------------

devtools::load_all("Gfreq")
library(reshape2)
library(vcfR)
library(zeallot)

## ---- Paths -----------------------------------------------------------------

input_dir <- "Gfreq/example_data/input_files"
out_dir <- "dev/outputs/estimate_genotype_frequencies"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_file <- file.path(out_dir, "genotype_frequencies.tsv")

## ---- Run pipeline -----------------------------------------------------------

res <- estimate_genotype_frequencies_pipeline(
  genotyping_vcf_path = file.path(input_dir, "comp_genotypes.vcf"),
  allele_freqs_path = file.path(input_dir, "mix_ref_all_freqs.tsv"),
  snp_depths_path = file.path(input_dir, "mix_read_depths.tsv"),
  lib_names_corresp_path = file.path(input_dir, "comp_libnames_corresp.tsv"),
  min_depth = 40,
  output_file_name = out_file
)
