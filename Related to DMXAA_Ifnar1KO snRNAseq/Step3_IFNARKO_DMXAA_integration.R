
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(dplyr)
library(cowplot)
library(reshape2)
library(MAST)

#load in data from Cell Ranger or other counts data ====

#for loading Cell Ranger counts:
setwd("/athena/ganlab/scratch/lif4001/DMXAA/DF_2ndRound")
Ctrl_1 <- readRDS(file = "Ctrl_1_singlets_PCA.rds")
Ctrl_2 <- readRDS(file = "Ctrl_2_singlets_PCA.rds")
Ctrl_3 <- readRDS(file = "Ctrl_3_singlets_PCA.rds")
Ctrl_4 <- readRDS(file = "Ctrl_4_singlets_PCA.rds")

DMXAA_1 <- readRDS(file = "DMXAA_1_singlets_PCA.rds")
DMXAA_2 <- readRDS(file = "DMXAA_2_singlets_PCA.rds")
#DMXAA_3 <- readRDS(file = "DMXAA_3_singlets_PCA.rds")
DMXAA_4 <- readRDS(file = "DMXAA_4_singlets_PCA.rds")

setwd("/athena/ganlab/scratch/lif4001/IFNARKO/DF_2ndRound")
Vehicle_1 <- readRDS(file = "Vehicle_1_singlets_PCA.rds")
Vehicle_2 <- readRDS(file = "Vehicle_2_singlets_PCA.rds")
Vehicle_3 <- readRDS(file = "Vehicle_3_singlets_PCA.rds")
Vehicle_4 <- readRDS(file = "Vehicle_4_singlets_PCA.rds")

IFNARKO_1 <- readRDS(file = "IFNARKO_1_singlets_PCA.rds")
IFNARKO_2 <- readRDS(file = "IFNARKO_2_singlets_PCA.rds")
IFNARKO_3 <- readRDS(file = "IFNARKO_3_singlets_PCA.rds")
IFNARKO_4 <- readRDS(file = "IFNARKO_4_singlets_PCA.rds")

setwd("/athena/ganlab/scratch/lif4001/IFNARKO/integration")
Ctrl <- c(Ctrl_1, Ctrl_2, Ctrl_3, Ctrl_4)
anchors_Ctrl <- FindIntegrationAnchors(object.list = Ctrl, dims = 1:30)
Ctrl_integrated <- IntegrateData(anchorset = anchors_Ctrl, dims = 1:30)
rm(Ctrl_1, Ctrl_2, Ctrl_3, Ctrl_4, Ctrl)

DMX <- c(DMXAA_1, DMXAA_2, DMXAA_4)
anchors_DMX <- FindIntegrationAnchors(object.list = DMX, dims = 1:30)
DMX_integrated <- IntegrateData(anchorset = anchors_DMX, dims = 1:30)
rm(DMXAA_1, DMXAA_2, DMXAA_4, DMX)

Vehicle <- c(Vehicle_1, Vehicle_2, Vehicle_3, Vehicle_4)
anchors_Vehicle <- FindIntegrationAnchors(object.list = Vehicle, dims = 1:30)
Vehicle_integrated <- IntegrateData(anchorset = anchors_Vehicle, dims = 1:30)
rm(Vehicle_1, Vehicle_2, Vehicle_3, Vehicle_4, Vehicle)

IFNAR <- c(IFNARKO_1, IFNARKO_2, IFNARKO_3, IFNARKO_4)
anchors_IFNAR <- FindIntegrationAnchors(object.list = IFNAR, dims = 1:30)
IFNAR_integrated <- IntegrateData(anchorset = anchors_IFNAR, dims = 1:30)
rm(IFNARKO_1, IFNARKO_2, IFNARKO_3, IFNARKO_4, IFNAR)

IFNARKO <- c(Ctrl_integrated, DMX_integrated, Vehicle_integrated, IFNAR_integrated)
anchors_IFNARKO <- FindIntegrationAnchors(object.list = IFNARKO, dims = 1:30)
IFNARKO_integrated <- IntegrateData(anchorset = anchors_IFNARKO, dims = 1:30)
rm(Ctrl_integrated, DMX_integrated, Vehicle_integrated, IFNAR_integrated, IFNARKO)


#saveRDS(IFNARKO_integrated, file = "IFNARKO_integrated.rds")

DefaultAssay(IFNARKO_integrated) <- 'integrated'

# IFNARKO_integrated <- NormalizeData(IFNARKO_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
# IFNARKO_integrated <- FindVariableFeatures(IFNARKO_integrated, selection.method = "vst", nfeatures = 3000)

IFNARKO_integrated <- ScaleData(IFNARKO_integrated, verbose = FALSE)
IFNARKO_integrated <- RunPCA(IFNARKO_integrated, features = VariableFeatures(object = IFNARKO_integrated), verbose = FALSE)

IFNARKO_integrated <- FindNeighbors(IFNARKO_integrated, dims = 1:15)
IFNARKO_integrated <- FindClusters(IFNARKO_integrated, resolution = 0.1)
IFNARKO_integrated <- RunUMAP(IFNARKO_integrated, dims = 1: 15)

DefaultAssay(IFNARKO_integrated) <- 'RNA'
IFNARKO_integrated <- NormalizeData(IFNARKO_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
IFNARKO_integrated <- ScaleData(IFNARKO_integrated, features = rownames(IFNARKO_integrated))

#saveRDS(IFNARKO_integrated, file = 'IFNARKO_integrated_PCA_0.1.rds')
#IFNARKO_integrated <- readRDS(file = "IFNARKO_integrated_PCA_0.1.rds")

IFNARKO_integrated$Condition <- factor(x = IFNARKO_integrated$Condition, levels = c("Ctrl","DMXAA","IFNARKO_Vehicle","IFNARKO_DMXAA"))
IFNARKO_integrated$Sample_Name <- factor(x = IFNARKO_integrated$Sample_Name, levels = c("Ctrl_1","Ctrl_2","Ctrl_3","Ctrl_4",
                                                                                    "DMXAA_1","DMXAA_2","DMXAA_4",
                                                                                    "IFNARKO_Vehicle_1","IFNARKO_Vehicle_2","IFNARKO_Vehicle_3","IFNARKO_Vehicle_4",
                                                                                    "IFNARKO_DMXAA_1","IFNARKO_DMXAA_2","IFNARKO_DMXAA_3","IFNARKO_DMXAA_4"))

pdf("IFNARKO_QC.pdf", width=9, height=4)
Idents(IFNARKO_integrated) <- "Condition"
VlnPlot(object = IFNARKO_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()


Idents(IFNARKO_integrated) <- "Sample_Name"
pdf("IFNARKO_QC_Sample.pdf", width=12, height=4)

VlnPlot(object = IFNARKO_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()

Idents(IFNARKO_integrated) <- "seurat_clusters"
pdf("IFNARKO_integrated_umap.pdf", width=5, height=4)
DimPlot(IFNARKO_integrated, reduction = 'umap', label = T)
dev.off()
pdf("IFNARKO_integrated_umap_split_individual.pdf", width=11, height=10)
DimPlot(IFNARKO_integrated, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 4)
dev.off()
pdf("IFNARKO_integrated_umap_split_Condition.pdf", width=13, height=3)
DimPlot(IFNARKO_integrated, reduction = "umap", split.by = "Condition", label = T, ncol = 4)
dev.off()

write.csv(table(IFNARKO_integrated$seurat_clusters, IFNARKO_integrated$Sample_Name), "IFNARKO_cell_counts_cluster_by_sample.csv")

DefaultAssay(IFNARKO_integrated) <- 'RNA'

IFNARKO_markers <- FindAllMarkers(IFNARKO_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1, test.use = "MAST")
write.csv(IFNARKO_markers, "IFNARKO_markers.csv")


saveRDS(IFNARKO_integrated, file = 'IFNARKO_integrated_PCA_0.1.rds')

IFNARKO_markers <- read.csv(file = "IFNARKO_markers.csv", header=T,row.names =1)
top5 <- IFNARKO_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5$gene <- as.character(top5$gene)
pdf("IFNARKO_HeatMapTop5_0.1_new.pdf", width=24, height=16)
DoHeatmap(IFNARKO_integrated, features = top5$gene) + NoLegend()
dev.off()

#Add marker genes

sig_EN<-c("Snap25","Slc17a7", "Nrgn","Gad1", "Gad2","Plp1", "Mbp", "Mobp","Sntn","Aqp4", "Clu", "Aldoc", "Pla2g7","Cx3cr1", "P2ry12", "Csf1r",
          "Pdgfra", "Vcan", "Flt1","Vtn", "Igfbp7")
markers.to.plot <- as.matrix(sig_EN)
pdf("IFNARKO_annotation_combine.pdf", width=10, height=5)
DotPlot(object = IFNARKO_integrated, features = rev(x = markers.to.plot)) + RotatedAxis()
dev.off()

