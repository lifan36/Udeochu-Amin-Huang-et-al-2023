#install.packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(reshape2)

#load in data from Cell Ranger or other counts data ====

#for loading Cell Ranger counts:
setwd("/athena/ganlab/scratch/lif4001/Joe_cGas/data_analysis/DF_2ndRound")
cGAS_25 <- readRDS(file = "cGAS_25_singlets.rds")
cGAS_23 <- readRDS(file = "cGAS_23_singlets.rds")
cGAS_18 <- readRDS(file = "cGAS_18_singlets.rds")
cGAS_20 <- readRDS(file = "cGAS_20_singlets.rds")
cGAS_55 <- readRDS(file = "cGAS_55_singlets.rds")
cGAS_61 <- readRDS(file = "cGAS_61_singlets.rds")
cGAS_21 <- readRDS(file = "cGAS_21_singlets.rds")
cGAS_22 <- readRDS(file = "cGAS_22_singlets.rds")
cGAS_19 <- readRDS(file = "cGAS_19_singlets.rds")
cGAS_30 <- readRDS(file = "cGAS_30_singlets.rds")
cGAS_29 <- readRDS(file = "cGAS_29_singlets.rds")
cGAS_63 <- readRDS(file = "cGAS_63_singlets.rds")
cGAS_36 <- readRDS(file = "cGAS_36_singlets.rds")
cGAS_43 <- readRDS(file = "cGAS_43_singlets.rds")
cGAS_33 <- readRDS(file = "cGAS_33_singlets.rds")
cGAS_59 <- readRDS(file = "cGAS_59_singlets.rds")
cGAS_42 <- readRDS(file = "cGAS_42_singlets.rds")
cGAS_74 <- readRDS(file = "cGAS_74_singlets.rds")
cGAS_45 <- readRDS(file = "cGAS_45_singlets.rds")
cGAS_46 <- readRDS(file = "cGAS_46_singlets.rds")

setwd("/athena/ganlab/scratch/jcu4001/JCUPs19cGas/TestRunJU1/Joe_cGas_3rd")
cGAS_27 = readRDS(file = "cGAS_27_singlets_PCA.rds")
cGAS_41 = readRDS(file = "cGAS_41_singlets_PCA.rds")
cGAS_57 = readRDS(file = "cGAS_57_singlets_PCA.rds")
cGAS_62 = readRDS(file = "cGAS_62_singlets_PCA.rds")
cGAS_65 = readRDS(file = "cGAS_65_singlets_PCA.rds")
cGAS_66 = readRDS(file = "cGAS_66_singlets_PCA.rds")
cGAS_77 = readRDS(file = "cGAS_77_singlets_PCA.rds")

setwd("/athena/ganlab/scratch/lif4001/Joe_cGas/test")

Batch1 <- c(cGAS_25, cGAS_23, cGAS_18, cGAS_20, cGAS_55, cGAS_61)
anchors_Batch1 <- FindIntegrationAnchors(object.list = Batch1, dims = 1:30)
Batch1_integrated <- IntegrateData(anchorset = anchors_Batch1, dims = 1:30)
rm(ccGAS_25, cGAS_23, cGAS_18, cGAS_20, cGAS_55, cGAS_61, Batch1)

Batch2 <- c(cGAS_21, cGAS_22, cGAS_19, cGAS_30, cGAS_29, cGAS_63)
anchors_Batch2 <- FindIntegrationAnchors(object.list = Batch2, dims = 1:30)
Batch2_integrated <- IntegrateData(anchorset = anchors_Batch2, dims = 1:30)
rm(cGAS_21, cGAS_22, cGAS_19, cGAS_30, cGAS_29, cGAS_6, Batch2)

Batch3 <- c(cGAS_36, cGAS_43, cGAS_33, cGAS_59, cGAS_42, cGAS_74, cGAS_45, cGAS_46)
anchors_Batch3 <- FindIntegrationAnchors(object.list = Batch3, dims = 1:30)
Batch3_integrated <- IntegrateData(anchorset = anchors_Batch3, dims = 1:30)
rm(cGAS_36, cGAS_43, cGAS_33, cGAS_59, cGAS_42, cGAS_74, cGAS_45, cGAS_46, Batch3)

Batch4 <- c(cGAS_27, cGAS_41, cGAS_57, cGAS_62, cGAS_65, cGAS_66, cGAS_77)
anchors_Batch4 <- FindIntegrationAnchors(object.list = Batch4, dims = 1:30)
Batch4_integrated <- IntegrateData(anchorset = anchors_Batch4, dims = 1:30)
rm(cGAS_27, cGAS_41, cGAS_57, cGAS_62, cGAS_65, cGAS_66, cGAS_77, Batch4)


cGAS <- c(Batch1_integrated, Batch2_integrated, Batch3_integrated, Batch4_integrated)
anchors_cGAS <- FindIntegrationAnchors(object.list = cGAS, dims = 1:30)
cGAS_integrated <- IntegrateData(anchorset = anchors_cGAS, dims = 1:30)
rm(Batch1_integrated, Batch2_integrated, Batch3_integrated, Batch4_integrated)

pdf("cGAS_QC_byBatch.pdf", width=15, height=4)
Idents(cGAS_integrated) <- "Sample_Name"
VlnPlot(object = cGAS_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
Idents(cGAS_integrated) <- "orig.ident"
VlnPlot(object = cGAS_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
Idents(cGAS_integrated) <- "Condition_2"
VlnPlot(object = cGAS_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()

#saveRDS(cGAS_integrated, file = "cGAS_integrated_byBatch.rds")

DefaultAssay(cGAS_integrated) <- 'integrated'

# cGAS_integrated <- NormalizeData(cGAS_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
# cGAS_integrated <- FindVariableFeatures(cGAS_integrated, selection.method = "vst", nfeatures = 3000)

cGAS_integrated <- ScaleData(cGAS_integrated, verbose = FALSE)
cGAS_integrated <- RunPCA(cGAS_integrated, features = VariableFeatures(object = cGAS_integrated), verbose = FALSE)

cGAS_integrated <- FindNeighbors(cGAS_integrated, dims = 1:20)
cGAS_integrated <- FindClusters(cGAS_integrated, resolution = 0.1)
cGAS_integrated <- RunUMAP(cGAS_integrated, dims = 1: 20)

str(cGAS_integrated)

DefaultAssay(cGAS_integrated) <- 'RNA'
cGAS_integrated <- NormalizeData(cGAS_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
cGAS_integrated <- ScaleData(cGAS_integrated, features = rownames(cGAS_integrated))

pdf("cGAS_integrated_umap_byBatch.pdf", width=5, height=4)
DimPlot(cGAS_integrated, reduction = 'umap', label = T)
dev.off()
pdf("cGAS_integrated_umap_split_individual_byBatch.pdf", width=12, height=15)
DimPlot(cGAS_integrated, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 4)
dev.off()
pdf("cGAS_integrated_umap_split_Condition_byBatch.pdf", width=15, height=10)
DimPlot(cGAS_integrated, reduction = "umap", split.by = "Condition_2", label = T, ncol = 3)
dev.off()

saveRDS(cGAS_integrated, file = 'cGAS_integrated_PCA_0.1_byBatch.rds')

DefaultAssay(cGAS_integrated) <- 'RNA'

cGAS_markers <- FindAllMarkers(cGAS_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1, test.use = "MAST")
write.csv(cGAS_markers, "cGAS_markers_byBatch.csv")

cGAS_markers <- read.csv(file = "cGAS_markers_byBatch.csv", header=T,row.names =1)
top5 <- cGAS_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
top5$gene <- as.character(top5$gene)
pdf("cGAS_HeatMapTop5_0.1_byBatch.pdf", width=24, height=16)
DoHeatmap(cGAS_integrated, features = top5$gene) + NoLegend()
dev.off()

DefaultAssay(cGAS_integrated) <- 'RNA'
pdf("cGAS_umap_test_byBatch.pdf", width=8, height=6)
DimPlot(cGAS_integrated, reduction = 'umap', label = T)
dev.off()

#Add marker genes

sig_EN<-c("Plp1", "Mbp", "Mobp","Slc17a7", "Nrgn", "Gad1", "Gad2", "Clu", "Aldoc", 
          "Pla2g7", "Cx3cr1", "P2ry12", "Csf1r","Scrg1", "Pdgfra", "Vtn", "Igfbp7",
          "Bnc2", "Slc47a1", "Ttr","Prox1","Mpped1","Map3k15","Cdh24")
markers.to.plot <- as.matrix(sig_EN)
pdf("cGAS_annotation_combine_byBatch.pdf", width=10, height=5)
DotPlot(object = cGAS_integrated, features = rev(x = markers.to.plot)) + RotatedAxis()
dev.off()

clusterSummary = as.data.frame(table(Idents(cGAS_integrated), cGAS_integrated$Genotype))
write.csv(clusterSummary, 'Idents(cGAS_integrated_ClusterSummary_byGenotype.csv')

clusterSummary = as.data.frame(table(Idents(cGAS_integrated), cGAS_integrated$orig.ident))
write.csv(clusterSummary, 'Idents(cGAS_integrated_ClusterSummary_byOrigIdent.csv')
