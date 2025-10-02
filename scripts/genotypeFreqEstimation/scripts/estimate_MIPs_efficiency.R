#!/usr/bin/env Rscript

# Print full error traceback when running from command line
options(error = function() {
  traceback(2)
  q(status = 1)
})

# v ... to replace with CLI params ... v
## Input parameters
genotypeFreqEstimation_path="~/GitHub/mip-poolfreq/scripts/genotypeFreqEstimation"
depth_list=paste0(genotypeFreqEstimation_path, "/data_test/data/ALL_FREQ/depth_files.list")
hline=50
# ^ .................................. ^


## Import dependencies
user_lib <- Sys.getenv("R_LIBS_USER")
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE)
.libPaths(user_lib)

pkgs<-c(
  "optparse",
  "devtools",
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

script_dir <- dirname(this.path::this.path())
project_root <- normalizePath(file.path(script_dir, ".."))
setwd(project_root)

load_all(".")

## Parse arguments
option_list <- list(
  make_option(c("--depth_files_list", "-d"), type = "character", help = "Path to a list of depth files"),
  make_option(c("--hline", "-h"), type = "numeric", default = NULL, help = "Numeric value for a horizontal reference line (read depth value)"),
  make_option(c("--output_file", "-o"), type = "character", help = "Path to output file (PNG)")
)

opt <- parse_args(OptionParser(option_list = option_list))

depth_list <- opt$depth_files_list
hline <- opt$hline
output_file_name <- opt$output_file


required_files <- list(
  depth_list = depth_list
)

optional_files <- list(
  
)

check_missing_args(args = c(required_files, output_file_name = output_file_name))

check_input_files(
  required_files = required_files,
  optional_files = optional_files
)

## Plot depths per marker
depths<-merge_depth_files(
  depth_list = depth_list
)

plot_mip_effect_on_depth(
  depths=depths,
  hline=hline,
  out_file=output_file_name
)


