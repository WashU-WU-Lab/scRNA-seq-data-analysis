# Coding: utf-8
# Details: Processing and analyzing sc_bt.rds

rm(list=ls())
library(Seurat)
library(dplyr)
library(ggplot2)
library(Nebulosa)
library(harmony)


# Parameters
workdir <- "C:/Users/Federico/Desktop/WU_Lab/scRNA-seq-data-analysis/code/"  # current work directory
datadir <- "../data/sc_bt.rds"  # .rds data path
cell_cycle_genes <- "../data/cell_cycle/Mus_musculus.rda"  # cell-cycle gene path
selected <- c("BHT6", "BHT7")  # For cell selection
nFeature_RNA_min <- 250  # For cell filtering
nCount_RNA_min <- 500  # For cell filtering
mt_ratio_max <- 20  # For cell filtering
GenesPerUMI_min <- 0.75  # For cell filtering
ncell_min <- 3  # For gene filtering
n_pc <- 20  # PC component num used in KNN & UMAP
regress_out <- c("percent.mt",  # Unwanted variations for regressing out
                 "S.Score",     # S.Score & G2M.Score are about cell cycle
                 "G2M.Score")  
integration_method <- "Seurat"  # "Harmony" for harmony, 
                                # "Seurat" for integration
                                # "" for no integration
resolution <- c(0.4, 0.6, 0.8, 1.0, 1.4)  # For clustering
gene_list <- c("Thy1", "Cd34")  # Check gene expression
compute_marker <- FALSE  # Whether to compute marker genes for clusters
save_prefix <- "sc_bt_BHT6&7_mt&cc_seurat_npc_20"  # save file prefix (like [prefix].rds)
seed <- 2024  # Random seed


# Read data & pre-computation
set.seed(seed)
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  # Set work dir
setwd(workdir)
data <- readRDS(datadir)
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt-")
data$log10GenesPerUMI <- log10(data$nFeature_RNA) / log10(data$nCount_RNA) # Novelty score


# Selecting samples
sub_data <- subset(data, subset = hash.ID %in% selected)
## Filtering
### Cell-level filtering
VlnPlot(sub_data, features = c("nFeature_RNA", "nCount_RNA", 
                                  "percent.mt", "log10GenesPerUMI"), ncol = 4)
metadata <- sub_data@meta.data
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
print(paste("Original cell num:", ncol(sub_data)))
sub_data <- subset(sub_data, subset = nFeature_RNA >= nFeature_RNA_min & 
                                            nCount_RNA >= nCount_RNA_min &
                                            percent.mt <= mt_ratio_max &
                                            log10GenesPerUMI >= GenesPerUMI_min)
print(paste("Cell num after filtering:", ncol(sub_data)))
### Gene-level filtering
counts <- GetAssayData(object = sub_data, layer = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= ncell_min
print(paste("Original gene num:", nrow(sub_data)))
sub_data <- CreateSeuratObject(counts[keep_genes, ],
                                  meta.data = sub_data@meta.data)
print(paste("Cell num after filtering:", nrow(sub_data)))
## Preprocessing
sub_data <- NormalizeData(sub_data, normalization.method = "LogNormalize", 
                             scale.factor = 10000, layer = "counts")
sub_data <- FindVariableFeatures(sub_data, selection.method = "vst", 
                                    nfeatures = 2000, layer = "counts")
sub_data <- ScaleData(sub_data)

### Check whether regress out the cell cycle
load(cell_cycle_genes)
sub_data <- CellCycleScoring(sub_data, 
                             g2m.features = g2m_genes, 
                             s.features = s_genes)
sub_data <- RunPCA(sub_data, seed.use = NULL)
DimPlot(sub_data,
        reduction = "pca",
        group.by= "Phase")
DimPlot(sub_data,
        reduction = "pca",
        group.by= "Phase",
        split.by= "Phase")

### Check whether regress out the mt gene ratio
sub_data@meta.data$cls.mt <- cut(sub_data@meta.data$percent.mt, 
                                 breaks=c(-Inf, mt_ratio_max * 0.25, 
                                          mt_ratio_max * 0.5, 
                                          mt_ratio_max * 0.75, Inf), 
                                 labels=c("Low","Medium","Medium high", "High"))
DimPlot(sub_data,
        reduction = "pca",
        group.by = "cls.mt")
DimPlot(sub_data,
        reduction = "pca",
        group.by = "cls.mt",
        split.by = "cls.mt")

sub_data <- SCTransform(sub_data, vars.to.regress = regress_out)
sub_data <- RunPCA(sub_data, assay = "SCT", npcs = 50, seed.use = NULL)
DefaultAssay(sub_data) <- "SCT"
ElbowPlot(sub_data, ndims = 50)
sub_data <- RunUMAP(sub_data, dims = 1:n_pc, reduction = "pca", seed.use = NULL)
Idents(sub_data) <- "hash.ID"
DimPlot(sub_data)
# Integration
if(integration_method == ""){
  reduction <- "pca"
}else if(integration_method == "Harmony"){
  sub_data <- RunHarmony(sub_data, group.by.vars = "hash.ID", 
                            reduction = "pca", assay.use = "SCT", 
                            reduction.save = "harmony")
  reduction <- "harmony"
  sub_data <- RunUMAP(sub_data, reduction = reduction, dims = 1:n_pc, 
                      seed.use = NULL)
  Idents(sub_data) <- "hash.ID"
  DimPlot(sub_data)
}else if(integration_method == "Seurat"){
  seurat_obj_list <- c()
  for(id in selected){
    temp_data <- subset(sub_data, subset = hash.ID == id)
    temp_data <- NormalizeData(temp_data, normalization.method = "LogNormalize", 
                              scale.factor = 10000, layer = "counts")
    temp_data <- FindVariableFeatures(temp_data, selection.method = "vst", 
                                     nfeatures = 2000, layer = "counts")
    temp_data <- SCTransform(temp_data, vars.to.regress = regress_out)
    seurat_obj_list <- c(seurat_obj_list, temp_data)
  }
  integ_features <- SelectIntegrationFeatures(object.list = seurat_obj_list, 
                                              nfeatures = 3000) 
  seurat_obj_list <- PrepSCTIntegration(object.list = seurat_obj_list, 
                                        anchor.features = integ_features)
  integ_anchors <- FindIntegrationAnchors(object.list = seurat_obj_list, 
                                          normalization.method = "SCT", 
                                          anchor.features = integ_features)
  sub_data <- IntegrateData(anchorset = integ_anchors, 
                               normalization.method = "SCT")
  sub_data <- ScaleData(sub_data)
  sub_data <- RunPCA(sub_data, npcs = 50, seed.use = NULL)
  reduction <- "pca"
  sub_data <- RunUMAP(sub_data, reduction = reduction, dims = 1:n_pc, 
                      seed.use = NULL)
  Idents(sub_data) <- "hash.ID"
  DimPlot(sub_data)
}
sub_data <- FindNeighbors(sub_data, reduction = reduction, dims = 1:n_pc)
for(r in resolution){
  sub_data <- FindClusters(sub_data, resolution = r)
  print(DimPlot(sub_data, label = TRUE))
  if(compute_marker){
    markers <- FindAllMarkers(sub_data, only.pos = TRUE, min.pct = 0.25, 
                              logfc.threshold = 0.25)
    write.csv(markers, paste("../results/markers/", save_prefix, 
                             "_rsln_", r, ".csv", sep=""))
  }
}
DefaultAssay(sub_data) <- "RNA"
for(gene in gene_list){
  print(FeaturePlot(sub_data, gene))
  print(VlnPlot(sub_data, gene, layer="data"))
  print(plot_density(sub_data, gene, reduction = "umap"))
}
if(length(gene_list) > 1){
  plot_density(sub_data, features = gene_list, joint = TRUE, 
               reduction = "umap", combine = FALSE)[length(gene_list) + 1]
}
saveRDS(sub_data, paste("../results/rds/", save_prefix, ".rds", sep=""))
