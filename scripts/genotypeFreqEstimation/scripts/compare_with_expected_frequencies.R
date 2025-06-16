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
# allele_freqs_mixtures_path=paste0(genotypeFreqEstimation_path, "/tests/data/ALL_FREQ/27MIXTURES/ref_allelic_freqs_27mixtures_246SNPs.tsv")
# snp_depths_mixtures_path=paste0(genotypeFreqEstimation_path, "/tests/data/ALL_FREQ/27MIXTURES/total_depths_27mixtures_246SNPs.tsv")
# expected_freqs_mixtures_path=paste0(genotypeFreqEstimation_path, "/tests/data/infos_files/expectedGenoFreqs_mixtures.tsv")
# allele_freqs_components_path=paste0(genotypeFreqEstimation_path, "/tests/data/ALL_FREQ/4COMPONENTS/ref_allelic_freqs_4components_246SNPs.tsv")
# snp_depths_components_path=paste0(genotypeFreqEstimation_path, "/tests/data/ALL_FREQ/4COMPONENTS/total_depths_4components_246SNPs.tsv")
# expected_freqs_components_path=paste0(genotypeFreqEstimation_path, "/tests/data/infos_files/expectedGenoFreqs_components.tsv")
# lib_names_corresp_path=paste0(genotypeFreqEstimation_path, "/tests/data/infos_files/corresp_comp_genotypes_libnames.tsv")
# out_dir=paste0(genotypeFreqEstimation_path, "/tests/results/compare_with_expected_frequencies")

# ^ .................................. ^

## Import dependencies
user_lib <- Sys.getenv("R_LIBS_USER")
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE)
.libPaths(user_lib)

pkgs<-c(
  "optparse",
  "devtools",
  "this.path",
  "glue",
  "vcfR",
  "reshape2",
  "scales",
  "ggplot2",
  "multcompView"
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
  make_option(c("--allele_freqs_mix"), type = "character", help = "Path to allele frequencies file of mixtures (with a 'CHROM' and a 'POS' columns)"),
  make_option(c("--depths_mix"), type = "character", help = "Path to SNP depths file of mixtures (with a 'CHROM' and a 'POS' columns)"),
  make_option(c("--exp_freqs_mix"), type = "character", help = "Path to expected genotype frequencies file of mixtures"),
  make_option(c("--allele_freqs_comp"), type = "character", default = NULL, help = "Path to allele frequencies file of components (with a 'CHROM' and a 'POS' columns)"),
  make_option(c("--depths_comp"), type = "character", default = NULL, help = "Path to SNP depths file of components (with a 'CHROM' and a 'POS' columns)"),
  make_option(c("--exp_freqs_comp"), type = "character", default = NULL, help = "Path to expected genotype frequencies file of components"),
  make_option(c("--libs", "-l"), type = "character", default = NULL, help = "Path to library name correspondence file"),
  make_option(c("--out_dir", "-o"), type = "character", help = "Path to output directory")
)

opt <- parse_args(OptionParser(option_list = option_list))

genotyping_vcf_path <- opt$vcf
allele_freqs_mixtures_path <- opt$allele_freqs_mix
snp_depths_mixtures_path <- opt$depths_mix
expected_freqs_mixtures_path<- opt$exp_freqs_mix
allele_freqs_components_path <- opt$allele_freqs_comp
snp_depths_components_path <- opt$depths_comp
expected_freqs_components_path<- opt$exp_freqs_comp
lib_names_corresp_path <- opt$libs
out_dir <- opt$out_dir

required_files <- list(
  vcf = genotyping_vcf_path,
  allele_freqs_mix = allele_freqs_mixtures_path,
  depths_mix = snp_depths_mixtures_path,
  exp_freqs_mix = expected_freqs_mixtures_path,
  out_dir = out_dir
)

optional_files <- list(
  libs = lib_names_corresp_path,
  allele_freqs_comp = allele_freqs_components_path,
  depths_comp = snp_depths_components_path,
  exp_freqs_comp = expected_freqs_components_path
)

check_missing_args(args = required_files)

check_input_files(
  required_files = required_files,
  optional_files = optional_files
)



## Import input files
genotyping_matrix<-vcf_to_numeric_matrix(
  vcf_path = genotyping_vcf_path,
  lib_names_corresp_path = lib_names_corresp_path
)

allele_freqs_mixtures<-read.table(allele_freqs_mixtures_path, header=T)
snp_depths_mixtures<-read.table(snp_depths_mixtures_path, header=T)
expected_freqs_mixtures<-read.table(expected_freqs_mixtures_path, header=T)

allele_freqs_mixtures<-set_rownames_from_chrom_pos(allele_freqs_mixtures)
snp_depths_mixtures<-set_rownames_from_chrom_pos(snp_depths_mixtures)
expected_freqs_mixtures_melt<-melt_genotype_freqs(expected_freqs_mixtures, "ExpFreq")


if (!is.null(allele_freqs_components_path)){
  allele_freqs_components<-read.table(allele_freqs_components_path, header=T)
  snp_depths_components<-read.table(snp_depths_components_path, header=T)
  expected_freqs_components<-read.table(expected_freqs_components_path, header=T)
  
  allele_freqs_components<-set_rownames_from_chrom_pos(allele_freqs_components)
  snp_depths_components<-set_rownames_from_chrom_pos(snp_depths_components)
  expected_freqs_components_melt<-melt_genotype_freqs(expected_freqs_components, "ExpFreq")
  
}



## Analysis

# Create output sub-directories
general_subdir<-glue("{out_dir}/general")
weights_subdir<-glue("{out_dir}/weight_vector_effect")
subsampling_subdir<-glue("{out_dir}/snp_subsampling_effect")

for (subdir in c(general_subdir, weights_subdir, subsampling_subdir)){
  dir.create(subdir)
}

# Estimate genotype frequencies, compare to expected frequencies and compute biases
writeLines("\n\nEstimate genotype frequencies, compare to expected frequencies and compute biases:\n")
genotype_frequencies_mixtures<-estimate_genotype_freqs(
  genotyping_matrix = genotyping_matrix, 
  allele_freqs = allele_freqs_mixtures, 
  snp_depths = snp_depths_mixtures,
  expected_freqs_melt = expected_freqs_mixtures_melt,
  out_file = glue("{general_subdir}/est_geno_freqs_mixtures.tsv")
)

compute_stats_per_exp_freq(
  freqs_df = genotype_frequencies_mixtures,
  out_dir = general_subdir
)

if (! is.null(allele_freqs_components)){
  genotype_frequencies_components<-estimate_genotype_freqs(
    genotyping_matrix = genotyping_matrix, 
    allele_freqs = allele_freqs_components, 
    snp_depths = snp_depths_components,
    expected_freqs_melt = expected_freqs_components_melt,
    out_file = glue("{general_subdir}/est_geno_freqs_components.tsv")
  )
} else {
  genotype_frequencies_components<-NULL
}

plot_correlation_with_exp_freq(
  freqs_df = genotype_frequencies_mixtures,
  extra_freqs_df = genotype_frequencies_components,
  out_file = glue("{general_subdir}/correlation_plot.png")
)

compute_bias(
  freqs_df = genotype_frequencies_mixtures,
  out_dir = general_subdir
)

# Evaluate the effect of the weight vector
compare_weight_vector_effect(
  genotyping_matrix = genotyping_matrix,
  allele_freqs = allele_freqs_mixtures,
  snp_depths = snp_depths_mixtures,
  expected_freqs_melt = expected_freqs_mixtures_melt,
  out_dir = weights_subdir
)

# Evaluate the effect of SNP sub-sampling
compare_subsampling_effect(
  genotyping_matrix = genotyping_matrix,
  allele_freqs = allele_freqs_mixtures,
  snp_depths = snp_depths_mixtures,
  expected_freqs_melt = expected_freqs_mixtures_melt,
  step_size = 50,
  out_dir = subsampling_subdir
)

writeLines(c(
  "",
  "Output files saved in:",
  glue("- {general_subdir}"),
  glue("- {weights_subdir}"),
  glue("- {subsampling_subdir}"),
  ""
))
