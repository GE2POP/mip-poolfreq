#!/usr/bin/env Rscript

# Print full error traceback when running from command line
options(error = function() {
  traceback(2)
  q(status = 1)
})

# v ... to replace with CLI params ... v
## Input parameters
# genotypeFreqEstimation_path="~/GitHub/mip-poolfreq/scripts/genotypeFreqEstimation"
# genotyping_vcf_path=paste0(genotypeFreqEstimation_path, "/data_test/data/GENOTYPING_MATRIX/04__Genotype_Locus1_Sample_Locus2_Filtered.vcf")
# lib_names_corresp_path=paste0(genotypeFreqEstimation_path, "/data_test/data/infos_files/corresp_comp_genotypes_libnames.tsv")
# out_dir=paste0(genotypeFreqEstimation_path, "/data_test/results/plot_MAF_hist")

# ^ .................................. ^


## Import dependencies
user_lib <- Sys.getenv("R_LIBS_USER")
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE)
.libPaths(user_lib)

pkgs<-c(
  "optparse",
  "devtools",
  "this.path",
  "vcfR",
  "glue",
  "UpSetR"
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
  make_option(c("--libs", "-l"), type = "character", default = NULL, help = "Path to library name correspondence file"),
  make_option(c("--out_dir", "-o"), type = "character", help = "Path to output directory")
)

opt <- parse_args(OptionParser(option_list = option_list))

genotyping_vcf_path <- opt$vcf
lib_names_corresp_path <- opt$libs
out_dir <- opt$out_dir

required_files <- list(
  vcf = genotyping_vcf_path,
  out_dir = out_dir
)

optional_files <- list(
  libs = lib_names_corresp_path
)

check_missing_args(args = required_files)

check_input_files(
  required_files = required_files,
  optional_files = optional_files
)


## Import input file
vcf <- read.vcfR(genotyping_vcf_path)
majmin_gt_matrix<-vcf_to_majmin_num_gt_matrix(
  vcf_path = genotyping_vcf_path,
  lib_names_corresp_path = lib_names_corresp_path
)

## Analysis
plot_MAF_hist(
  vcf = vcf, 
  x_lim_values = c(0, 0.5),
  out_file = glue("{out_dir}/MAF_hist.png")
)

#majmin_gt_matrix <- majmin_gt_matrix[rowSums(majmin_gt_matrix == 0.5, na.rm = TRUE) == 0, , drop = FALSE]

plot_marker_set_intersections(
  numeric_matrix = majmin_gt_matrix,
  out_file = glue("{out_dir}/market_set_upsetplot.png")
)


