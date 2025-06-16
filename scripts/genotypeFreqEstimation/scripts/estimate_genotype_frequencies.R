#!/usr/bin/env Rscript

# Print full error traceback when running from command line
options(error = function() {
  traceback(2)
  q(status = 1)
})

# v ... to replace with CLI params ... v
## Input parameters
# genotypeFreqEstimation_path="~/MIPs/PUBLI/03_scripts/genotypeFreqEstimation"
# genotyping_vcf_path=paste0(genotypeFreqEstimation_path, "/tests/data/GENOTYPING_MATRIX/04__Genotype_Locus1_Sample_Locus2_Filtered.vcf")
# allele_freqs_path=paste0(genotypeFreqEstimation_path, "/tests/data/ALL_FREQ/27MIXTURES/ref_allelic_freqs_27mixtures_246SNPs.tsv")
# snp_depths_path=paste0(genotypeFreqEstimation_path, "/tests/data/ALL_FREQ/27MIXTURES/total_depths_27mixtures_246SNPs.tsv")
# lib_names_corresp_path=paste0(genotypeFreqEstimation_path, "/tests/data/infos_files/corresp_comp_genotypes_libnames.tsv")
# output_file_name=paste0(genotypeFreqEstimation_path, "/tests/results/estimate_genotype_frequencies/genotype_frequencies.tsv")
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
  "vcfR"
)

for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing missing package: ", pkg)
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

script_dir <- dirname(this.path::this.path())
project_root <- normalizePath(file.path(script_dir, ".."))
setwd(project_root)

load_all(".")

## Parse arguments
option_list <- list(
  make_option(c("--vcf", "-v"), type = "character", help = "Path to genotyping VCF"),
  make_option(c("--allele_freqs", "-a"), type = "character", help = "Path to allele frequencies file (with a 'CHROM' and a 'POS' columns)"),
  make_option(c("--depths", "-d"), type = "character", help = "Path to SNP depths file (with a 'CHROM' and a 'POS' columns)"),
  make_option(c("--libs", "-l"), type = "character", default = NULL, help = "Path to library name correspondence file"),
  make_option(c("--output_file", "-o"), type = "character", help = "Path to output file")
)

opt <- parse_args(OptionParser(option_list = option_list))

genotyping_vcf_path <- opt$vcf
allele_freqs_path <- opt$allele_freqs
snp_depths_path <- opt$depths
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

check_missing_args(args = c(required_files, output_file_name = output_file_name))

check_input_files(
  required_files = required_files,
  optional_files = optional_files
)



## Import input files
genotyping_matrix<-vcf_to_numeric_matrix(
  vcf_path = genotyping_vcf_path,
  lib_names_corresp_path = lib_names_corresp_path
)

allele_freqs<-read.table(allele_freqs_path, header=T)
snp_depths<-read.table(snp_depths_path, header=T)

allele_freqs<-set_rownames_from_chrom_pos(allele_freqs)
snp_depths<-set_rownames_from_chrom_pos(snp_depths)


## Analysis
genotype_frequencies<-estimate_genotype_freqs(
  genotyping_matrix = genotyping_matrix, 
  allele_freqs = allele_freqs, 
  snp_depths = snp_depths,
  out_file = output_file_name
)
