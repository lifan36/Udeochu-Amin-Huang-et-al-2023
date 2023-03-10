
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(dplyr)
library(cowplot)
library(reshape2)
library(MAST)
setwd("/athena/ganlab/scratch/lif4001/IFNARKO/integration")
Cluster_EN <- readRDS("IFNARKO_EN_subset.rds")
Cluster_IN <- readRDS("IFNARKO_IN_subset.rds")
Cluster_MG <- readRDS("IFNARKO_MG_subset.rds")
Cluster_AST <- readRDS("IFNARKO_AST_subset.rds")
Cluster_OL <- readRDS("IFNARKO_OL_subset.rds")
Cluster_OPC <- readRDS("IFNARKO_OPC_subset.rds")
Cluster_VC <- readRDS("IFNARKO_VC_subset.rds")

Idents(Cluster_EN) <- "Condition"
Idents(Cluster_IN) <- "Condition"
Idents(Cluster_MG) <- "Condition"
Idents(Cluster_AST) <- "Condition"
Idents(Cluster_OL) <- "Condition"
Idents(Cluster_OPC) <- "Condition"
Idents(Cluster_VC) <- "Condition"

Cluster_EN <- subset(Cluster_EN, idents=c("Ctrl","DMXAA"))
Cluster_IN <- subset(Cluster_IN, idents=c("Ctrl","DMXAA"))
Cluster_MG <- subset(Cluster_MG, idents=c("Ctrl","DMXAA"))
Cluster_AST <- subset(Cluster_AST, idents=c("Ctrl","DMXAA"))
Cluster_OL <- subset(Cluster_OL, idents=c("Ctrl","DMXAA"))
Cluster_OPC <- subset(Cluster_OPC, idents=c("Ctrl","DMXAA"))
Cluster_VC <- subset(Cluster_VC, idents=c("Ctrl","DMXAA"))

#######################################################################
setwd("/athena/ganlab/scratch/lif4001/IFNARKO/integration/subclustering_DMXAA")
MG <- Cluster_MG
DefaultAssay(MG) <- 'integrated'
MG <- ScaleData(MG, verbose = FALSE)
MG <- RunPCA(MG, features = VariableFeatures(object = MG), verbose = FALSE)
ElbowPlot(MG)
MG <- FindNeighbors(MG, dims = 1:15)
MG <- FindClusters(MG, resolution = 0.15)
MG <- RunUMAP(MG, dims = 1: 15)
# rename cluster
n <- dim(table(MG@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
MG@active.ident <- plyr::mapvalues(x = MG@active.ident, from = current.cluster.ids, to = new.cluster.ids)
MG@active.ident <- factor(MG@active.ident, levels=1:n)
saveRDS(MG, file = 'IFNARKO_MG_reclusted_res0.15.rds')
pdf("IFNARKO_MG_umap_res0.15.pdf", width=3.3, height=2.7)
DimPlot(MG, reduction = 'umap', label = T)
dev.off()
pdf("IFNARKO_MG_umap_Condition_res0.15.pdf", width=5.5, height=2.7)
DimPlot(MG, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("IFNARKO_MG_umap_Sample_res0.15.pdf", width=8, height=4)
DimPlot(MG, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 4)
dev.off()
DefaultAssay(MG) <- 'RNA'
IFNARKO_MG_markers <- FindAllMarkers(MG, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(IFNARKO_MG_markers, "IFNARKO_MG_markers_res0.15.csv")
write.csv(table(MG$seurat_clusters, MG$Sample_Name), "MG_subcluster_cell_counts_res0.15.csv")
#######################################################################
setwd("/athena/ganlab/scratch/lif4001/IFNARKO/integration/subclustering_DMXAA")
MG <- Cluster_MG
DefaultAssay(MG) <- 'integrated'
MG <- ScaleData(MG, verbose = FALSE)
MG <- RunPCA(MG, features = VariableFeatures(object = MG), verbose = FALSE)
ElbowPlot(MG)
MG <- FindNeighbors(MG, dims = 1:15)
MG <- FindClusters(MG, resolution = 0.2)
MG <- RunUMAP(MG, dims = 1: 15)
# rename cluster
n <- dim(table(MG@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
MG@active.ident <- plyr::mapvalues(x = MG@active.ident, from = current.cluster.ids, to = new.cluster.ids)
MG@active.ident <- factor(MG@active.ident, levels=1:n)
saveRDS(MG, file = 'IFNARKO_MG_reclusted_res0.2.rds')
pdf("IFNARKO_MG_umap_res0.2.pdf", width=3.3, height=2.7)
DimPlot(MG, reduction = 'umap', label = T)
dev.off()
pdf("IFNARKO_MG_umap_Condition_res0.2.pdf", width=5.5, height=2.7)
DimPlot(MG, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("IFNARKO_MG_umap_Sample_res0.2.pdf", width=8, height=4)
DimPlot(MG, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 4)
dev.off()
DefaultAssay(MG) <- 'RNA'
IFNARKO_MG_markers <- FindAllMarkers(MG, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(IFNARKO_MG_markers, "IFNARKO_MG_markers_res0.2.csv")
write.csv(table(MG$seurat_clusters, MG$Sample_Name), "MG_subcluster_cell_counts_res0.2.csv")
#######################################################################
setwd("/athena/ganlab/scratch/lif4001/IFNARKO/integration/subclustering_DMXAA")
MG <- Cluster_MG
DefaultAssay(MG) <- 'integrated'
MG <- ScaleData(MG, verbose = FALSE)
MG <- RunPCA(MG, features = VariableFeatures(object = MG), verbose = FALSE)
ElbowPlot(MG)
MG <- FindNeighbors(MG, dims = 1:15)
MG <- FindClusters(MG, resolution = 0.3)
MG <- RunUMAP(MG, dims = 1: 15)
# rename cluster
n <- dim(table(MG@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
MG@active.ident <- plyr::mapvalues(x = MG@active.ident, from = current.cluster.ids, to = new.cluster.ids)
MG@active.ident <- factor(MG@active.ident, levels=1:n)
saveRDS(MG, file = 'IFNARKO_MG_reclusted_res0.3.rds')
pdf("IFNARKO_MG_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(MG, reduction = 'umap', label = T)
dev.off()
pdf("IFNARKO_MG_umap_Condition_res0.3.pdf", width=5.5, height=2.7)
DimPlot(MG, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("IFNARKO_MG_umap_Sample_res0.3.pdf", width=8, height=4)
DimPlot(MG, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 4)
dev.off()
DefaultAssay(MG) <- 'RNA'
IFNARKO_MG_markers <- FindAllMarkers(MG, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(IFNARKO_MG_markers, "IFNARKO_MG_markers_res0.3.csv")
write.csv(table(MG$seurat_clusters, MG$Sample_Name), "MG_subcluster_cell_counts_res0.3.csv")
#######################################################################
AST <- Cluster_AST
DefaultAssay(AST) <- 'integrated'
AST <- ScaleData(AST, verbose = FALSE)
AST <- RunPCA(AST, features = VariableFeatures(object = AST), verbose = FALSE)
ElbowPlot(AST)
AST <- FindNeighbors(AST, dims = 1:15)
AST <- FindClusters(AST, resolution = 0.15)
AST <- RunUMAP(AST, dims = 1: 15)
# rename cluster
n <- dim(table(AST@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
AST@active.ident <- plyr::mapvalues(x = AST@active.ident, from = current.cluster.ids, to = new.cluster.ids)
AST@active.ident <- factor(AST@active.ident, levels=1:n)
saveRDS(AST, file = 'IFNARKO_AST_reclusted_res0.15.rds')
pdf("IFNARKO_AST_umap.pdf", width=3.3, height=2.7)
DimPlot(AST, reduction = 'umap', label = T)
dev.off()
pdf("IFNARKO_AST_umap_Condition.pdf", width=5.5, height=2.7)
DimPlot(AST, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("IFNARKO_AST_umap_Sample.pdf", width=8, height=4)
DimPlot(AST, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 4)
dev.off()
DefaultAssay(AST) <- 'RNA'
IFNARKO_AST_markers <- FindAllMarkers(AST, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(IFNARKO_AST_markers, "IFNARKO_AST_markers.csv")
write.csv(table(AST$seurat_clusters, AST$Sample_Name), "AST_subcluster_cell_counts.csv")
setwd("/athena/ganlab/scratch/lif4001/IFNARKO/integration/subclustering_DMXAA")
AST <- Cluster_AST
DefaultAssay(AST) <- 'integrated'
AST <- ScaleData(AST, verbose = FALSE)
AST <- RunPCA(AST, features = VariableFeatures(object = AST), verbose = FALSE)
ElbowPlot(AST)
AST <- FindNeighbors(AST, dims = 1:15)
AST <- FindClusters(AST, resolution = 0.2)
AST <- RunUMAP(AST, dims = 1: 15)
# rename cluster
n <- dim(table(AST@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
AST@active.ident <- plyr::mapvalues(x = AST@active.ident, from = current.cluster.ids, to = new.cluster.ids)
AST@active.ident <- factor(AST@active.ident, levels=1:n)
saveRDS(AST, file = 'IFNARKO_AST_reclusted_res0.2.rds')
pdf("IFNARKO_AST_umap_res0.2.pdf", width=3.3, height=2.7)
DimPlot(AST, reduction = 'umap', label = T)
dev.off()
pdf("IFNARKO_AST_umap_Condition_res0.2.pdf", width=5.5, height=2.7)
DimPlot(AST, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("IFNARKO_AST_umap_Sample_res0.2.pdf", width=8, height=4)
DimPlot(AST, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 4)
dev.off()
DefaultAssay(AST) <- 'RNA'
IFNARKO_AST_markers <- FindAllMarkers(AST, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(IFNARKO_AST_markers, "IFNARKO_AST_markers_res0.2.csv")
write.csv(table(AST$seurat_clusters, AST$Sample_Name), "AST_subcluster_cell_counts_res0.2.csv")
#######################################################################
setwd("/athena/ganlab/scratch/lif4001/IFNARKO/integration/subclustering_DMXAA")
AST <- Cluster_AST
DefaultAssay(AST) <- 'integrated'
AST <- ScaleData(AST, verbose = FALSE)
AST <- RunPCA(AST, features = VariableFeatures(object = AST), verbose = FALSE)
ElbowPlot(AST)
AST <- FindNeighbors(AST, dims = 1:15)
AST <- FindClusters(AST, resolution = 0.3)
AST <- RunUMAP(AST, dims = 1: 15)
# rename cluster
n <- dim(table(AST@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
AST@active.ident <- plyr::mapvalues(x = AST@active.ident, from = current.cluster.ids, to = new.cluster.ids)
AST@active.ident <- factor(AST@active.ident, levels=1:n)
saveRDS(AST, file = 'IFNARKO_AST_reclusted_res0.3.rds')
pdf("IFNARKO_AST_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(AST, reduction = 'umap', label = T)
dev.off()
pdf("IFNARKO_AST_umap_Condition_res0.3.pdf", width=5.5, height=2.7)
DimPlot(AST, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("IFNARKO_AST_umap_Sample_res0.3.pdf", width=8, height=4)
DimPlot(AST, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 4)
dev.off()
DefaultAssay(AST) <- 'RNA'
IFNARKO_AST_markers <- FindAllMarkers(AST, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(IFNARKO_AST_markers, "IFNARKO_AST_markers_res0.3.csv")
write.csv(table(AST$seurat_clusters, AST$Sample_Name), "AST_subcluster_cell_counts_res0.3.csv")
#######################################################################
#######################################################################
OL <- Cluster_OL
DefaultAssay(OL) <- 'integrated'
OL <- ScaleData(OL, verbose = FALSE)
OL <- RunPCA(OL, features = VariableFeatures(object = OL), verbose = FALSE)
ElbowPlot(OL)
OL <- FindNeighbors(OL, dims = 1:15)
OL <- FindClusters(OL, resolution = 0.15)
OL <- RunUMAP(OL, dims = 1: 15)
# rename cluster
n <- dim(table(OL@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
OL@active.ident <- plyr::mapvalues(x = OL@active.ident, from = current.cluster.ids, to = new.cluster.ids)
OL@active.ident <- factor(OL@active.ident, levels=1:n)
saveRDS(OL, file = 'IFNARKO_OL_reclusted_res0.15.rds')
pdf("IFNARKO_OL_umap.pdf", width=3.3, height=2.7)
DimPlot(OL, reduction = 'umap', label = T)
dev.off()
pdf("IFNARKO_OL_umap_Condition.pdf", width=5.5, height=2.7)
DimPlot(OL, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("IFNARKO_OL_umap_Sample.pdf", width=8, height=4)
DimPlot(OL, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 4)
dev.off()
DefaultAssay(OL) <- 'RNA'
IFNARKO_OL_markers <- FindAllMarkers(OL, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(IFNARKO_OL_markers, "IFNARKO_OL_markers.csv")
write.csv(table(OL$seurat_clusters, OL$Sample_Name), "OL_subcluster_cell_counts.csv")
#######################################################################
setwd("/athena/ganlab/scratch/lif4001/IFNARKO/integration/subclustering_DMXAA")
OL <- Cluster_OL
DefaultAssay(OL) <- 'integrated'
OL <- ScaleData(OL, verbose = FALSE)
OL <- RunPCA(OL, features = VariableFeatures(object = OL), verbose = FALSE)
ElbowPlot(OL)
OL <- FindNeighbors(OL, dims = 1:15)
OL <- FindClusters(OL, resolution = 0.2)
OL <- RunUMAP(OL, dims = 1: 15)
# rename cluster
n <- dim(table(OL@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
OL@active.ident <- plyr::mapvalues(x = OL@active.ident, from = current.cluster.ids, to = new.cluster.ids)
OL@active.ident <- factor(OL@active.ident, levels=1:n)
saveRDS(OL, file = 'IFNARKO_OL_reclusted_res0.2.rds')
pdf("IFNARKO_OL_umap_res0.2.pdf", width=3.3, height=2.7)
DimPlot(OL, reduction = 'umap', label = T)
dev.off()
pdf("IFNARKO_OL_umap_Condition_res0.2.pdf", width=5.5, height=2.7)
DimPlot(OL, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("IFNARKO_OL_umap_Sample_res0.2.pdf", width=8, height=4)
DimPlot(OL, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 4)
dev.off()
DefaultAssay(OL) <- 'RNA'
IFNARKO_OL_markers <- FindAllMarkers(OL, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(IFNARKO_OL_markers, "IFNARKO_OL_markers_res0.2.csv")
write.csv(table(OL$seurat_clusters, OL$Sample_Name), "OL_subcluster_cell_counts_res0.2.csv")
#######################################################################
setwd("/athena/ganlab/scratch/lif4001/IFNARKO/integration/subclustering_DMXAA")
OL <- Cluster_OL
DefaultAssay(OL) <- 'integrated'
OL <- ScaleData(OL, verbose = FALSE)
OL <- RunPCA(OL, features = VariableFeatures(object = OL), verbose = FALSE)
ElbowPlot(OL)
OL <- FindNeighbors(OL, dims = 1:15)
OL <- FindClusters(OL, resolution = 0.3)
OL <- RunUMAP(OL, dims = 1: 15)
# rename cluster
n <- dim(table(OL@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
OL@active.ident <- plyr::mapvalues(x = OL@active.ident, from = current.cluster.ids, to = new.cluster.ids)
OL@active.ident <- factor(OL@active.ident, levels=1:n)
saveRDS(OL, file = 'IFNARKO_OL_reclusted_res0.3.rds')
pdf("IFNARKO_OL_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(OL, reduction = 'umap', label = T)
dev.off()
pdf("IFNARKO_OL_umap_Condition_res0.3.pdf", width=5.5, height=2.7)
DimPlot(OL, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("IFNARKO_OL_umap_Sample_res0.3.pdf", width=8, height=4)
DimPlot(OL, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 4)
dev.off()
DefaultAssay(OL) <- 'RNA'
IFNARKO_OL_markers <- FindAllMarkers(OL, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(IFNARKO_OL_markers, "IFNARKO_OL_markers_res0.3.csv")
write.csv(table(OL$seurat_clusters, OL$Sample_Name), "OL_subcluster_cell_counts_res0.3.csv")
#######################################################################
OPC <- Cluster_OPC
DefaultAssay(OPC) <- 'integrated'
OPC <- ScaleData(OPC, verbose = FALSE)
OPC <- RunPCA(OPC, features = VariableFeatures(object = OPC), verbose = FALSE)
ElbowPlot(OPC)
OPC <- FindNeighbors(OPC, dims = 1:15)
OPC <- FindClusters(OPC, resolution = 0.15)
OPC <- RunUMAP(OPC, dims = 1: 15)
# rename cluster
n <- dim(table(OPC@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
OPC@active.ident <- plyr::mapvalues(x = OPC@active.ident, from = current.cluster.ids, to = new.cluster.ids)
OPC@active.ident <- factor(OPC@active.ident, levels=1:n)
saveRDS(OPC, file = 'IFNARKO_OPC_reclusted_res0.15.rds')
pdf("IFNARKO_OPC_umap.pdf", width=3.3, height=2.7)
DimPlot(OPC, reduction = 'umap', label = T)
dev.off()
pdf("IFNARKO_OPC_umap_Condition.pdf", width=5.5, height=2.7)
DimPlot(OPC, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("IFNARKO_OPC_umap_Sample.pdf", width=8, height=4)
DimPlot(OPC, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 4)
dev.off()
DefaultAssay(OPC) <- 'RNA'
IFNARKO_OPC_markers <- FindAllMarkers(OPC, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(IFNARKO_OPC_markers, "IFNARKO_OPC_markers.csv")
write.csv(table(OPC$seurat_clusters, OPC$Sample_Name), "OPC_subcluster_cell_counts.csv")
#######################################################################
setwd("/athena/ganlab/scratch/lif4001/IFNARKO/integration/subclustering_DMXAA")
OPC <- Cluster_OPC
DefaultAssay(OPC) <- 'integrated'
OPC <- ScaleData(OPC, verbose = FALSE)
OPC <- RunPCA(OPC, features = VariableFeatures(object = OPC), verbose = FALSE)
ElbowPlot(OPC)
OPC <- FindNeighbors(OPC, dims = 1:15)
OPC <- FindClusters(OPC, resolution = 0.2)
OPC <- RunUMAP(OPC, dims = 1: 15)
# rename cluster
n <- dim(table(OPC@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
OPC@active.ident <- plyr::mapvalues(x = OPC@active.ident, from = current.cluster.ids, to = new.cluster.ids)
OPC@active.ident <- factor(OPC@active.ident, levels=1:n)
saveRDS(OPC, file = 'IFNARKO_OPC_reclusted_res0.2.rds')
pdf("IFNARKO_OPC_umap_res0.2.pdf", width=3.3, height=2.7)
DimPlot(OPC, reduction = 'umap', label = T)
dev.off()
pdf("IFNARKO_OPC_umap_Condition_res0.2.pdf", width=5.5, height=2.7)
DimPlot(OPC, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("IFNARKO_OPC_umap_Sample_res0.2.pdf", width=8, height=4)
DimPlot(OPC, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 4)
dev.off()
DefaultAssay(OPC) <- 'RNA'
IFNARKO_OPC_markers <- FindAllMarkers(OPC, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(IFNARKO_OPC_markers, "IFNARKO_OPC_markers_res0.2.csv")
write.csv(table(OPC$seurat_clusters, OPC$Sample_Name), "OPC_subcluster_cell_counts_res0.2.csv")
#######################################################################
setwd("/athena/ganlab/scratch/lif4001/IFNARKO/integration/subclustering_DMXAA")
OPC <- Cluster_OPC
DefaultAssay(OPC) <- 'integrated'
OPC <- ScaleData(OPC, verbose = FALSE)
OPC <- RunPCA(OPC, features = VariableFeatures(object = OPC), verbose = FALSE)
ElbowPlot(OPC)
OPC <- FindNeighbors(OPC, dims = 1:15)
OPC <- FindClusters(OPC, resolution = 0.3)
OPC <- RunUMAP(OPC, dims = 1: 15)
# rename cluster
n <- dim(table(OPC@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
OPC@active.ident <- plyr::mapvalues(x = OPC@active.ident, from = current.cluster.ids, to = new.cluster.ids)
OPC@active.ident <- factor(OPC@active.ident, levels=1:n)
saveRDS(OPC, file = 'IFNARKO_OPC_reclusted_res0.3.rds')
pdf("IFNARKO_OPC_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(OPC, reduction = 'umap', label = T)
dev.off()
pdf("IFNARKO_OPC_umap_Condition_res0.3.pdf", width=5.5, height=2.7)
DimPlot(OPC, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("IFNARKO_OPC_umap_Sample_res0.3.pdf", width=8, height=4)
DimPlot(OPC, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 4)
dev.off()
DefaultAssay(OPC) <- 'RNA'
IFNARKO_OPC_markers <- FindAllMarkers(OPC, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(IFNARKO_OPC_markers, "IFNARKO_OPC_markers_res0.3.csv")
write.csv(table(OPC$seurat_clusters, OPC$Sample_Name), "OPC_subcluster_cell_counts_res0.3.csv")
#######################################################################
setwd("/athena/ganlab/scratch/lif4001/IFNARKO/integration/subclustering_DMXAA")
IN <- Cluster_IN
DefaultAssay(IN) <- 'integrated'
IN <- ScaleData(IN, verbose = FALSE)
IN <- RunPCA(IN, features = VariableFeatures(object = IN), verbose = FALSE)
ElbowPlot(IN)
IN <- FindNeighbors(IN, dims = 1:15)
IN <- FindClusters(IN, resolution = 0.1)
IN <- RunUMAP(IN, dims = 1: 15)
# rename cluster
n <- dim(table(IN@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
IN@active.ident <- plyr::mapvalues(x = IN@active.ident, from = current.cluster.ids, to = new.cluster.ids)
IN@active.ident <- factor(IN@active.ident, levels=1:n)
saveRDS(IN, file = 'IFNARKO_IN_reclusted_res0.1.rds')
pdf("IFNARKO_IN_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(IN, reduction = 'umap', label = T)
dev.off()
pdf("IFNARKO_IN_umap_Condition_res0.3.pdf", width=5.5, height=2.7)
DimPlot(IN, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("IFNARKO_IN_umap_Sample_res0.3.pdf", width=8, height=4)
DimPlot(IN, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 4)
dev.off()
DefaultAssay(IN) <- 'RNA'
IFNARKO_IN_markers <- FindAllMarkers(IN, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(IFNARKO_IN_markers, "IFNARKO_IN_markers_res0.3.csv")
write.csv(table(IN$seurat_clusters, IN$Sample_Name), "IN_subcluster_cell_counts_res0.3.csv")
#######################################################################
setwd("/athena/ganlab/scratch/lif4001/IFNARKO/integration/subclustering_DMXAA")
EN <- Cluster_EN
DefaultAssay(EN) <- 'integrated'
EN <- ScaleData(EN, verbose = FALSE)
EN <- RunPCA(EN, features = VariableFeatures(object = EN), verbose = FALSE)
ElbowPlot(EN)
EN <- FindNeighbors(EN, dims = 1:15)
EN <- FindClusters(EN, resolution = 0.1)
EN <- RunUMAP(EN, dims = 1: 15)
# rename cluster
n <- dim(table(EN@active.ident))
current.cluster.ids <- c(0, seq(1:(n-1)))
new.cluster.ids <- current.cluster.ids + 1
EN@active.ident <- plyr::mapvalues(x = EN@active.ident, from = current.cluster.ids, to = new.cluster.ids)
EN@active.ident <- factor(EN@active.ident, levels=1:n)
saveRDS(EN, file = 'IFNARKO_EN_reclusted_res0.1.rds')
pdf("IFNARKO_EN_umap_res0.3.pdf", width=3.3, height=2.7)
DimPlot(EN, reduction = 'umap', label = T)
dev.off()
pdf("IFNARKO_EN_umap_Condition_res0.3.pdf", width=5.5, height=2.7)
DimPlot(EN, reduction = "umap", split.by = "Condition", label = T, ncol = 2)
dev.off()
pdf("IFNARKO_EN_umap_Sample_res0.3.pdf", width=8, height=4)
DimPlot(EN, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 4)
dev.off()
DefaultAssay(EN) <- 'RNA'
IFNARKO_EN_markers <- FindAllMarkers(EN, logfc.threshold = 0.1, test.use = "MAST",min.pct = 0.25, only.pos = T)
write.csv(IFNARKO_EN_markers, "IFNARKO_EN_markers_res0.3.csv")
write.csv(table(EN$seurat_clusters, EN$Sample_Name), "EN_subcluster_cell_counts_res0.3.csv")
#######################################################################
