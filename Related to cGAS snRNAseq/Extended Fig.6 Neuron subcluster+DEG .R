#install.packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(reshape2)
library(MAST)

setwd()

Cluster_EN <- readRDS(file = 'EN_cGAS_integrated_Annotated.rds')
DefaultAssay(Cluster_EN) <- 'integrated'
all.genes <- rownames(Cluster_EN)
Cluster_EN <- ScaleData(Cluster_EN, features = all.genes)
Cluster_EN <- FindVariableFeatures(object = Cluster_EN)
Cluster_EN <- RunPCA(Cluster_EN, features = VariableFeatures(object = Cluster_EN))

DimPlot(Cluster_EN, reduction = "pca")
ElbowPlot(Cluster_EN)

Cluster_EN <- FindNeighbors(Cluster_EN, dims = 1:10)
Cluster_EN <- FindClusters(Cluster_EN, resolution = 0.1)
Cluster_EN <- RunUMAP(Cluster_EN, dims = 1:10)

Idents(Cluster_EN) <- 'seurat_clusters'
DimPlot(Cluster_EN, reduction = "umap", label = T)

Cluster_EN <- RenameIdents(Cluster_EN, '0' = "1", '1' = "2", '2' = "3",
                           '3' = "4", '4' = "5", '5' = "6", '6' = "7", '7' = "8",'8' = "9",'9' = "10",'10' = "11", '11' = "12")
pdf("UMAP_ExcitatoryNeuron.pdf", width=4, height=4)
DimPlot(Cluster_EN, reduction = "umap", label = F)
dev.off()

saveRDS(Cluster_EN, 'EN_cGAS_subclustered.rds')

#Inhibitory Neuron
Cluster_IN <- readRDS(file = 'IN_cGAS_integrated_Annotated_new.rds')
DefaultAssay(Cluster_IN) <- 'RNA'
all.genes <- rownames(Cluster_IN)
Cluster_IN <- NormalizeData(Cluster_IN, normalization.method = 'LogNormalize', scale.factor = 10000)
Cluster_IN <- ScaleData(Cluster_IN, features = all.genes)
Cluster_IN <- FindVariableFeatures(object = Cluster_IN)
Cluster_IN <- RunPCA(Cluster_IN, features = VariableFeatures(object = Cluster_IN))

top10 <- head(VariableFeatures(Cluster_IN), 10)
plot1 <- VariableFeaturePlot(Cluster_IN)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

DimPlot(Cluster_IN, reduction = "pca")
ElbowPlot(Cluster_IN)

Cluster_IN <- FindNeighbors(Cluster_IN, dims = 1:15)
Cluster_IN <- FindClusters(Cluster_IN, resolution = 0.1)
Cluster_IN <- RunUMAP(Cluster_IN, dims = 1:15)

Idents(Cluster_IN) <- 'seurat_clusters'
Cluster_IN <- RenameIdents(Cluster_IN, '0' = "1", '1' = "2", '2' = "3", '3' = "4", '4' = "5", '5' = "6", '6' = "7", '7' = "8",'8' = "9",'9' = "10",'10' = "11")
pdf("UMAP_InhibitoryNeuron.pdf", width=4, height=4)
DimPlot(Cluster_IN, reduction = "umap", label = F)
dev.off()

DimPlot(Cluster_IN, reduction = "umap", label = T)
saveRDS(Cluster_IN, 'IN_cGAS_reclustered.rds')
DimPlot(Cluster_IN, reduction = "umap", label = T)
Cluster_IN_sub <- subset(Cluster_IN, subset= Genotype == 'P301S tau: +; cGAS: +/+'| Genotype =='P301S tau: +; cGAS: +/-'| Genotype == 'P301S tau: +; cGAS: -/-' | Genotype =='P301S tau: -; cGAS: +/+' )

#DEG analysis
Idents(Cluster_EN) <- 'Genotype'
PS19KOvsPS19WT_DE <- FindMarkers(Cluster_EN, ident.1 = "P301S tau: +; cGAS: -/-", ident.2 = "P301S tau: +; cGAS: +/+", logfc.threshold = 0,
                                 test.use = "MAST", assay ='RNA')
write.csv(PS19KOvsPS19WT_DE, file = 'PS19KOvsPS19WT_DE_MAST_RNAassay_ExcitatoryNeuron.csv')
# refer to supplementary table S7

Idents(Cluster_IN) <- 'Genotype'
PS19KOvsPS19WT_DE <- FindMarkers(Cluster_IN, ident.1 = "P301S tau: +; cGAS: -/-", ident.2 = "P301S tau: +; cGAS: +/+", logfc.threshold = 0,
                                 test.use = "MAST", assay ='RNA')
write.csv(PS19KOvsPS19WT_DE, file = 'PS19KOvsPS19WT_DE_MAST_RNAassay_InhibitoryNeuron.csv')
# refer to supplementary table S7


r = DotPlot(Cluster_EN, features = c('Pdzrn3','Dpp10', 'Efna5','Mef2c', 'Satb2', 
                                     'Satb1', 'Nrg1','Pcdh15','Lrrc4c'), assay = 'RNA') + RotatedAxis()

r+scale_colour_gradient2(low = "darkblue", mid = "white", high = "darkred")
dev.off()

