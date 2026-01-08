#!/usr/bin/env Rscript

## ---- Load packages in dev mode ----------------------------------------------

devtools::load_all("Gfreq")
library(ggplot2)
library(glue)
library(multcompView)
library(reshape2)
library(scales)
library(vcfR)
library(zeallot)

## ---- Paths -----------------------------------------------------------------

input_dir <- "Gfreq/example_data/input_files"
out_dir <- "dev/outputs/frequency_accuracy"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---- Run pipeline -----------------------------------------------------------

res <- frequency_accuracy_pipeline(
  genotyping_vcf_path = file.path(input_dir, "comp_genotypes.vcf"),
  allele_freqs_mixtures_path = file.path(input_dir, "mix_ref_all_freqs.tsv"),
  snp_depths_mixtures_path = file.path(input_dir, "mix_read_depths.tsv"),
  expected_freqs_mixtures_path = file.path(input_dir, "mix_exp_geno_freqs.tsv"),
  allele_freqs_components_path = file.path(input_dir, "comp_ref_all_freqs.tsv"),
  snp_depths_components_path = file.path(input_dir, "comp_read_depths.tsv"),
  expected_freqs_components_path = file.path(input_dir, "comp_exp_geno_freqs.tsv"),
  lib_names_corresp_path = file.path(input_dir, "comp_libnames_corresp.tsv"),
  min_depth = 40,
  subsampling_step = 30,
  out_dir = out_dir
)
