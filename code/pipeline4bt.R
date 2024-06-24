# Coding: utf-8
# Details: Processing and analyzing sc_bt.rds

rm(list=ls())
library(Seurat)
library(dplyr)
library(ggplot2)
library(Nebulosa)
library(harmony)

set.seed(2024)
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  # Set work dir
data <- readRDS("../data/sc_bt.rds")
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt-")
data$log10GenesPerUMI <- log10(data$nFeature_RNA) / log10(data$nCount_RNA) # Novelty score


# For BHT 6 & 7
sub_data_67 <- subset(data, subset = hash.ID %in% c("BHT6", "BHT7"))
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
  geom_vline(xintercept = 250)
metadata %>% 
  ggplot(aes(x=nCount_RNA, color=hash.ID, fill= hash.ID)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)
metadata %>% 
  ggplot(aes(x=percent.mt, color=hash.ID, fill=hash.ID)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 20)
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = hash.ID, fill = hash.ID)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.75)
print(paste("Original cell num:", ncol(sub_data_67)))
sub_data_67 <- subset(sub_data_67, subset = nFeature_RNA >= 250 & 
                                            nCount_RNA >= 500 &
                                            percent.mt <= 20 &
                                            log10GenesPerUMI >= 0.75)
print(paste("Cell num after filtering:", ncol(sub_data_67)))
### Gene-level filtering
counts <- GetAssayData(object = sub_data_67, layer = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
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
sub_data_67 <- RunUMAP(sub_data_67, dims = 1:30, reduction = "pca")
Idents(sub_data_67) <- "hash.ID"
DimPlot(sub_data_67)
# Integration
sub_data_67 <- RunHarmony(sub_data_67, group.by.vars = "hash.ID", 
                          reduction = "pca", assay.use = "SCT", 
                          reduction.save = "harmony")
sub_data_67 <- RunUMAP(sub_data_67, reduction = "harmony", dims = 1:30)
Idents(sub_data_67) <- "hash.ID"
DimPlot(sub_data_67)
sub_data_67 <- FindNeighbors(sub_data_67, reduction = "harmony", dims = 1:30)
sub_data_67 <- FindClusters(sub_data_67, resolution = 0.3)
DimPlot(sub_data_67, label = TRUE)
DefaultAssay(sub_data_67) <- "RNA"
FeaturePlot(sub_data_67, "Thy1", slot = "data")
FeaturePlot(sub_data_67, "Cd34", slot = "data")
VlnPlot(sub_data_67, "Thy1", slot="data")
VlnPlot(sub_data_67, "Cd34", slot="data")
plot_density(sub_data_67, "Thy1", reduction = "umap")
plot_density(sub_data_67, "Cd34", reduction = "umap")
plot_density(sub_data_67, features = c("Thy1", "Cd34"), joint = TRUE, 
             reduction = "umap", combine = FALSE)[3]
markers <- FindAllMarkers(sub_data_67, only.pos = TRUE, min.pct = 0.25, 
                          logfc.threshold = 0.25)
write.csv(markers, "../processed/sc_bt_BHT6&7_Xinyang_20240624.csv")
saveRDS(sub_data_67, "../processed/sc_bt_BHT6&7_Xinyang_20240624.Rds")

