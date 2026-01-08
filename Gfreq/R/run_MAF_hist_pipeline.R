#' Minor allele frequency analysis pipeline
#'
#' Run the complete pipeline to compute and plot minor allele frequency (MAF)
#' distributions and marker set intersections from a genotyping VCF.
#'
#' @param genotyping_vcf_path Path to genotyping VCF
#' @param lib_names_corresp_path Optional path to library name correspondence file
#' @param out_dir Output directory
#'
#' @return Invisibly returns the output file path
#' @export
MAF_hist_pipeline <- function(
  genotyping_vcf_path,
  lib_names_corresp_path = NULL,
  out_dir
) {
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

  vcf <- read.vcfR(genotyping_vcf_path)
  majmin_gt_matrix<-vcf_to_majmin_num_gt_matrix(
    vcf_path = genotyping_vcf_path,
    lib_names_corresp_path = lib_names_corresp_path
  )

  plot_MAF_hist(
    vcf = vcf,
    x_lim_values = c(0, 0.5),
    out_dir = out_dir
  )

  upset_file<-glue("{out_dir}/market_set_upsetplot.png")
  plot_marker_set_intersections(
    numeric_matrix = majmin_gt_matrix,
    out_file = upset_file
  )

  invisible(list(
    upset_plot = upset_file
  ))
}
