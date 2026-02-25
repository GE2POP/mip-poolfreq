#!/usr/bin/env Rscript

## ---- Load packages in dev mode ----------------------------------------------

devtools::load_all("Gfreq")
library(glue)
library(UpSetR)
library(vcfR)

## ---- Paths -----------------------------------------------------------------

input_dir <- "Gfreq/example_data/input_files"
out_dir <- "dev/outputs/maf"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---- Run pipeline -----------------------------------------------------------

res <- maf_pipeline(
  genotyping_vcf_path = file.path(input_dir, "comp_genotypes.vcf"),
  lib_names_corresp_path = file.path(input_dir, "comp_libnames_corresp.tsv"),
  out_dir = out_dir
)
