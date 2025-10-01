#!/usr/bin/env Rscript

# Print full error traceback when running from command line
options(error = function() {
  traceback(2)
  q(status = 1)
})

# v ... to replace with CLI params ... v
## Input parameters
genotypeFreqEstimation_path="~/GitHub/mip-poolfreq/scripts/genotypeFreqEstimation"
genotyping_vcf_path=paste0(genotypeFreqEstimation_path, "/data_test/data/GENOTYPING_MATRIX/04__Genotype_Locus1_Sample_Locus2_Filtered.vcf")
lib_names_corresp_path=paste0(genotypeFreqEstimation_path, "/data_test/data/infos_files/corresp_comp_genotypes_libnames.tsv")
output_file_name=paste0(genotypeFreqEstimation_path, "/data_test/results/plot_MAF_hist/MAF_hist.png")
out_dir=paste0(genotypeFreqEstimation_path, "/data_test/results/plot_MAF_hist")

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
  make_option(c("--output_file", "-o"), type = "character", help = "Path to output file")
)

opt <- parse_args(OptionParser(option_list = option_list))

genotyping_vcf_path <- opt$vcf
output_file_name <- opt$output_file

required_files <- list(
  vcf = genotyping_vcf_path
)

check_missing_args(args = c(required_files, output_file_name = output_file_name))

check_input_files(
  required_files = required_files
)

## Import input file
vcf <- read.vcfR(genotyping_vcf_path)
# vcf_num <- vcf_to_numeric_matrix(
#   vcf_path = genotyping_vcf_path,
#   lib_names_corresp_path = lib_names_corresp_path
# )

## Analysis
plot_MAF_hist(
  vcf = vcf, 
  x_lim_values = c(0, 0.5), 
  out_file = output_file_name
)


num_matrix_to_major_minor <- function(numeric_matrix) {
  f_ref <- rowMeans(numeric_matrix, na.rm = TRUE)
  inv <- !is.na(f_ref) & (f_ref < 0.5)
  numeric_matrix[inv, ] <- 1 - numeric_matrix[inv, ]
  numeric_matrix
}

vcf_to_majmin_numeric_matrix <- function(vcf_path, lib_names_corresp_path = NULL) {
  numeric_matrix <- vcf_to_numeric_matrix(
    vcf_path = vcf_path,
    lib_names_corresp_path = lib_names_corresp_path
  )
  majmin_num_matrix<-num_matrix_to_major_minor(numeric_matrix)
  return(majmin_num_matrix)
}
vcf_num<-vcf_to_majmin_numeric_matrix(
  vcf_path = genotyping_vcf_path,
  lib_names_corresp_path = lib_names_corresp_path
)
#vcf_num <- vcf_num[rowSums(vcf_num == 0.5, na.rm = TRUE) == 0, , drop = FALSE]

plot_marker_set_intersections(vcf_num)


