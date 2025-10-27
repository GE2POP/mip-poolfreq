#!/usr/bin/env Rscript

# Print full error traceback when running from command line
options(error = function() {
  traceback(2)
  q(status = 1)
})

# v ... to replace with CLI params ... v
## Input parameters
# genotypeFreqEstimation_path="~/GitHub/mip-poolfreq/scripts/genotypeFreqEstimation/example_data/input_files"
# genotyping_vcf_path=paste0(genotypeFreqEstimation_path, "/comp_genotypes.vcf")
# allele_freqs_path=paste0(genotypeFreqEstimation_path, "/mix_ref_all_freqs.tsv")
# snp_depths_path=paste0(genotypeFreqEstimation_path, "/mix_read_depths.tsv")
# min_depth=40
# lib_names_corresp_path=paste0(genotypeFreqEstimation_path, "/comp_libnames_corresp.tsv")
# output_file_name="~/GitHub/mip-poolfreq/scripts/genotypeFreqEstimation/example_data/output_files/estimate_genotype_frequencies/genotype_frequencies.tsv"
# ^ .................................. ^

## Import dependencies
user_lib <- Sys.getenv("R_LIBS_USER")
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE)
.libPaths(user_lib)

pkgs<-c(
  "optparse",
  "devtools",
  "this.path",
  "reshape2",
  "vcfR",
  "zeallot"
)

for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing missing package: ", pkg)
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

user_dir <- getwd()
script_dir <- dirname(this.path::this.path())
project_root <- normalizePath(file.path(script_dir, ".."))
setwd(project_root)
load_all(".")
setwd(user_dir)

## Parse arguments
option_list <- list(
  make_option(c("--vcf", "-v"), type = "character", help = "Path to genotyping VCF"),
  make_option(c("--allele_freqs", "-a"), type = "character", help = "Path to allele frequencies file (with a 'CHROM' and a 'POS' columns)"),
  make_option(c("--depths", "-d"), type = "character", help = "Path to SNP depths file (with a 'CHROM' and a 'POS' columns)"),
  make_option(c("--min_depth", "-t"), type = "numeric", default = 0, help = "Minimum read depth threshold. SNPs with at least one genotype < threshold we be removed."),
  make_option(c("--libs", "-l"), type = "character", default = NULL, help = "Path to library name correspondence file"),
  make_option(c("--output_file", "-o"), type = "character", help = "Path to output file")
)

opt <- parse_args(OptionParser(option_list = option_list))

genotyping_vcf_path <- opt$vcf
allele_freqs_path <- opt$allele_freqs
snp_depths_path <- opt$depths
min_depth <- opt$min_depth
lib_names_corresp_path <- opt$libs
output_file_name <- opt$output_file

required_files <- list(
  vcf = genotyping_vcf_path,
  allele_freqs = allele_freqs_path,
  depths = snp_depths_path
)

optional_files <- list(
  libs = lib_names_corresp_path
)

check_missing_args(args = c(required_files, output_file_name = output_file_name, min_depth = min_depth))

check_input_files(
  required_files = required_files,
  optional_files = optional_files
)



## Import input files
inputs <- load_inputs(
  genotyping_vcf_path = genotyping_vcf_path,
  lib_names_corresp_path = lib_names_corresp_path,
  allele_freqs_mixtures_path = allele_freqs_path,
  snp_depths_mixtures_path = snp_depths_path,
  min_depth = min_depth
)

## Analysis
genotype_frequencies<-estimate_genotype_freqs(
  genotyping_matrix = inputs$genotyping_matrix, 
  allele_freqs = inputs$allele_freqs_mixtures, 
  snp_depths = inputs$snp_depths_mixtures,
  out_file = output_file_name
)
