#' Depth per marker analysis pipeline
#'
#' This function merges SNP depth files, computes mean depth per
#' marker and generates depth distribution plots.
#'
#' @param depth_list Path to a file listing depth files
#' @param hline Optional numeric value for a horizontal reference line
#' @param out_dir Output directory
#'
#' @return Invisibly returns a list of output file paths
#' @export
depth_pipeline <- function(
  depth_list,
  hline = NULL,
  out_dir
) {
  required_files <- list(
    depth_list = depth_list,
    out_dir = out_dir
  )

  optional_files <- list()

  check_missing_args(args = required_files)

  check_input_files(
    required_files = required_files,
    optional_files = optional_files
  )

  depths<-merge_depth_files(
    depth_list = depth_list
  )

  compute_mean_depth_per_marker(
    depths=depths,
    out_file=glue("{out_dir}/mean_depth_per_marker.tsv")
  )

  plot_depth_per_marker(
    depths=depths,
    hline=hline,
    out_file=glue("{out_dir}/depth_per_marker_boxplots.png")
  )


  invisible(list(
    mean_depth = depths
  ))
}
