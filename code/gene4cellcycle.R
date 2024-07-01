# Coding: utf-8
# Details: generating cell cycle genes (for homo sapiens or mus musculus)

rm(list=ls())
library(AnnotationHub)
library(ensembldb)
library(dplyr)


# Set work dir
workdir <- "C:/Users/Federico/Desktop/WU_Lab/scRNA-seq-data-analysis/code/"
setwd(workdir)

# Read gene data downloaded from https://github.com/hbc/tinyatlas/tree/master/cell_cycle
cell_cycle_genes <- read.csv("../data/cell_cycle/Homo_sapiens.csv")

# Connect to AnnotationHub
ah <- AnnotationHub()

# Access the Ensembl database for organism
ahDb <- query(ah, 
              pattern = c("Homo sapiens", "EnsDb"), 
              ignore.case = TRUE)

# Acquire the latest annotation files
id <- ahDb %>%
  mcols() %>%
  rownames() %>%
  tail(n = 1)

# Download the appropriate Ensembldb database
edb <- ah[[id]]

# Extract gene-level information from database
annotations <- genes(edb, 
                     return.type = "data.frame")

# Select annotations of interest
annotations <- annotations %>%
  dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)

# Get gene names for Ensembl IDs for each gene
cell_cycle_markers <- dplyr::left_join(cell_cycle_genes, annotations, by = c("geneID" = "gene_id"))

# Acquire the S phase genes
s_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "S") %>%
  pull("gene_name")

# Acquire the G2M phase genes        
g2m_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "G2/M") %>%
  pull("gene_name")

# Save data
save(s_genes, g2m_genes, file = "../data/cell_cycle/Homo_sapiens.rda")
