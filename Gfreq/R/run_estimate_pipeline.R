#' Estimate genotype frequencies pipeline
#'
#' Run the complete pipeline to estimate genotype frequencies from a genotyping
#' VCF, allele frequencies and SNP depths.
#'
#' @param genotyping_vcf_path Path to genotyping VCF
#' @param allele_freqs_path Path to allele frequencies file
#' @param snp_depths_path Path to SNP depths file
#' @param lib_names_corresp_path Optional path to library name correspondence file
#' @param min_depth Minimum read depth threshold
#' @param output_file_name Output file path
#'
#' @return Invisibly returns the output file path
#' @export
estimate_pipeline <- function(
  genotyping_vcf_path,
  allele_freqs_path,
  snp_depths_path,
  lib_names_corresp_path = NULL,
  min_depth = 0,
  output_file_name
) {
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


  inputs <- load_inputs(
    genotyping_vcf_path = genotyping_vcf_path,
    lib_names_corresp_path = lib_names_corresp_path,
    allele_freqs_mixtures_path = allele_freqs_path,
    snp_depths_mixtures_path = snp_depths_path,
    min_depth = min_depth
  )

  genotype_frequencies<-estimate_genotype_freqs(
    genotyping_matrix = inputs$genotyping_matrix,
    allele_freqs = inputs$allele_freqs_mixtures,
    snp_depths = inputs$snp_depths_mixtures,
    out_file = output_file_name
  )

  invisible(output_file_name)
}
