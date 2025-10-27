#!/usr/bin/env Rscript

# Print full error traceback when running from command line
options(error = function() {
  traceback(2)
  q(status = 1)
})

# v ... to replace with CLI params ... v
## Input parameters
# genotypeFreqEstimation_path="~/GitHub/mip-poolfreq/scripts/genotypeFreqEstimation"
# depth_list=paste0(genotypeFreqEstimation_path, "/data_test/data/ALL_FREQ/depth_files.list")
# hline=50
# out_dir=paste0(genotypeFreqEstimation_path, "/data_test/results/plot_depth_per_marker")

# ^ .................................. ^


## Import dependencies
user_lib <- Sys.getenv("R_LIBS_USER")
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE)
.libPaths(user_lib)

pkgs<-c(
  "optparse",
  "devtools",
  "glue",
  "ggplot2",
  "reshape2"
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
  make_option(c("--depth_files_list", "-d"), type = "character", help = "Path to a list of depth files"),
  make_option(c("--hline", "-l"), type = "numeric", default = NULL, help = "Numeric value for a horizontal reference line (read depth value)"),
  make_option(c("--out_dir", "-o"), type = "character", help = "Path to output directory")
)

opt <- parse_args(OptionParser(option_list = option_list))

depth_list <- opt$depth_files_list
hline <- opt$hline
out_dir <- opt$out_dir


required_files <- list(
  depth_list = depth_list,
  out_dir = out_dir
)

optional_files <- list(
  
)

check_missing_args(args = required_files)

check_input_files(
  required_files = required_files,
  optional_files = optional_files
)

## Plot depths per marker
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

