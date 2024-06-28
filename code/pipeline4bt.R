# Coding: utf-8
# Details: Processing and analyzing sc_bt.rds

rm(list=ls())
library(Seurat)
library(dplyr)
library(ggplot2)
library(Nebulosa)
library(harmony)


# Parameters
workdir <- ""
datadir <- "../data/sc_bt.rds"
selected <- c("BHT6", "BHT7")  # For cell selection
nFeature_RNA_min <- 250  # For cell filtering
nCount_RNA_min <- 500  # For cell filtering
mt_ratio_max <- 20  # For cell filtering
GenesPerUMI_min <- 0.75  # For cell filtering
ncell_min <- 10  # For gene filtering
n_pc <- 30  # PC component num used in KNN & UMAP
integration_method <- "Harmony"  # "Harmony" for harmony, 
                                 # "Seurat" for integration
                                 # "" for no integration
resolution <- 0.8  # For clustering
gene_list <- c("Thy1", "Cd34")  # Check gene expression
compute_marker <- FALSE  # Whether computing marker genes for clusters
save_prefix <- "../processed/sc_bt_BHT6&7_npc_30_rsln_0.3"
seed <- 2024  # Random seed


# Read data & pre-computation
set.seed(seed)
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  # Set work dir
setwd <- workdir
data <- readRDS(datadir)
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt-")
data$log10GenesPerUMI <- log10(data$nFeature_RNA) / log10(data$nCount_RNA) # Novelty score


# For BHT 6 & 7
sub_data_67 <- subset(data, subset = hash.ID %in% selected)
## Filtering
### Cell-level filtering
VlnPlot(sub_data_67, features = c("nFeature_RNA", "nCount_RNA", 
                                  "percent.mt", "log10GenesPerUMI"), ncol = 4)
metadata <- sub_data_67@meta.data
metadata %>% 
  ggplot(aes(x=nFeature_RNA, olor=hash.ID, fill= hash.ID)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = nFeature_RNA_min)
metadata %>% 
  ggplot(aes(x=nCount_RNA, color=hash.ID, fill= hash.ID)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = nCount_RNA_min)
metadata %>% 
  ggplot(aes(x=percent.mt, color=hash.ID, fill=hash.ID)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = mt_ratio_max)
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = hash.ID, fill = hash.ID)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = GenesPerUMI_min)
print(paste("Original cell num:", ncol(sub_data_67)))
sub_data_67 <- subset(sub_data_67, subset = nFeature_RNA >= nFeature_RNA_min & 
                                            nCount_RNA >= nCount_RNA_min &
                                            percent.mt <= mt_ratio_max &
                                            log10GenesPerUMI >= GenesPerUMI_min)
print(paste("Cell num after filtering:", ncol(sub_data_67)))
### Gene-level filtering
counts <- GetAssayData(object = sub_data_67, layer = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= ncell_min
print(paste("Original gene num:", nrow(sub_data_67)))
sub_data_67 <- CreateSeuratObject(counts[keep_genes, ],
                                  meta.data = sub_data_67@meta.data)
print(paste("Cell num after filtering:", nrow(sub_data_67)))
## Preprocessing
sub_data_67 <- NormalizeData(sub_data_67, normalization.method = "LogNormalize", 
                             scale.factor = 10000, layer = "counts")
sub_data_67 <- FindVariableFeatures(sub_data_67, selection.method = "vst", 
                                    nfeatures = 2000, layer = "counts")
sub_data_67 <- ScaleData(sub_data_67, features = VariableFeatures(sub_data_67))
sub_data_67 <- SCTransform(sub_data_67, vars.to.regress = "percent.mt")
sub_data_67 <- RunPCA(sub_data_67, assay = "SCT", npcs = 50)
DefaultAssay(sub_data_67) <- "SCT"
ElbowPlot(sub_data_67, ndims = 50)
sub_data_67 <- RunUMAP(sub_data_67, dims = 1:n_pc, reduction = "pca")
Idents(sub_data_67) <- "hash.ID"
DimPlot(sub_data_67)
# Integration
if(integration_method == ""){
  reduction <- "pca"
}else if(integration_method == "Harmony"){
  sub_data_67 <- RunHarmony(sub_data_67, group.by.vars = "hash.ID", 
                            reduction = "pca", assay.use = "SCT", 
                            reduction.save = "harmony")
  reduction <- "harmony"
  sub_data_67 <- RunUMAP(sub_data_67, reduction = reduction, dims = 1:n_pc)
  Idents(sub_data_67) <- "hash.ID"
  DimPlot(sub_data_67)
}else if(integration_method == "Seurat"){
  seurat_obj_list <- c()
  for(id in selected){
    seurat_obj_list <- c(seurat_obj_list, 
                         subset(sub_data_67, subset = hash.ID == id))
  }
  integ_features <- SelectIntegrationFeatures(object.list = seurat_obj_list, 
                                              nfeatures = 3000) 
  seurat_obj_list <- PrepSCTIntegration(object.list = seurat_obj_list, 
                                        anchor.features = integ_features)
  integ_anchors <- FindIntegrationAnchors(object.list = seurat_obj_list, 
                                          normalization.method = "SCT", 
                                          anchor.features = integ_features)
  sub_data_67 <- IntegrateData(anchorset = integ_anchors, 
                               normalization.method = "SCT")
  reduction <- "pca"
  sub_data_67 <- RunPCA(object = sub_data_67)
  sub_data_67 <- RunUMAP(sub_data_67, reduction = reduction, dims = 1:n_pc)
  Idents(sub_data_67) <- "hash.ID"
  DimPlot(sub_data_67)
}
sub_data_67 <- FindNeighbors(sub_data_67, reduction = reduction, dims = 1:n_pc)
sub_data_67 <- FindClusters(sub_data_67, resolution = resolution)
DimPlot(sub_data_67, label = TRUE)
DefaultAssay(sub_data_67) <- "RNA"
for(gene in gene_list){
  FeaturePlot(sub_data_67, gene)
  VlnPlot(sub_data_67, gene, layer="data")
  plot_density(sub_data_67, gene, reduction = "umap")
}
if(length(gene_list) > 1){
  plot_density(sub_data_67, features = gene_list, joint = TRUE, 
               reduction = "umap", combine = FALSE)[length(gene_list) + 1]
}
if(compute_marker){
  markers <- FindAllMarkers(sub_data_67, only.pos = TRUE, min.pct = 0.25, 
                            logfc.threshold = 0.25)
  write.csv(markers, paste(save_prefix, ".csv", sep=""))
}
saveRDS(sub_data_67, paste(save_prefix, ".rds", sep=""))
