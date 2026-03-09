#!/usr/bin/env Rscript

## ---- Load packages in dev mode ----------------------------------------------

devtools::load_all("Gfreq")
library(ggplot2)
library(glue)
library(reshape2)


## ---- Paths -----------------------------------------------------------------

out_dir <- "dev/outputs/depth"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

depth_list <- "dev/depth_files.list"

## ---- Run pipeline -----------------------------------------------------------

res <- depth_pipeline(
  depth_list = depth_list,
  hline = 50,
  out_dir = out_dir
)
