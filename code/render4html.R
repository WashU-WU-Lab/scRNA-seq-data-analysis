# Coding: utf-8
# Details: rendering .html files from R script

rm(list=ls())
library(rmarkdown)


# Parameters
workdir <- "C:/Users/Federico/Desktop/WU_Lab/scRNA-seq-data-analysis/code/"
output_dir <- "../results/html/"
outputname <- "pipeline4bt_BHT6&7_mt&cc_seurat_npc_20.html"

# Set work dir
setwd(workdir)

# Render and save results
rmarkdown::render(input = "pipeline4bt.R", 
                  output_format = "html_document",  # or "pdf_document", "all" for all formats
                  output_file = outputname,
                  output_dir = output_dir)
