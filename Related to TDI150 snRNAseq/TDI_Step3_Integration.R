#setwd("/athena/ganlab/scratch/jcu4001/JCUPs19cGas/TestRunJU1/TDI_sNucSeq150")
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(RColorBrewer)

C90 <- readRDS(file = "singlets_PCA_C90.rds")
D45<- readRDS(file = "singlets_PCA_D45.rds")
D46<- readRDS(file = "singlets_PCA_D46.rds")
D27<- readRDS(file = "singlets_PCA_D27.rds")
C95 <- readRDS(file = "singlets_PCA_C95.rds")
D18<- readRDS(file = "singlets_PCA_D18.rds")
D17<- readRDS(file = "singlets_PCA_D17.rds")
D15<- readRDS(file = "singlets_PCA_D15.rds")
C96 <- readRDS(file = "singlets_PCA_C96.rds")
D13<- readRDS(file = "singlets_PCA_D13.rds")
D22<- readRDS(file = "singlets_PCA_D22.rds")
D50<- readRDS(file = "singlets_PCA_D50.rds")
C92 <- readRDS(file = "singlets_PCA_C92.rds")
C94 <- readRDS(file = "singlets_PCA_C94.rds")
D49<- readRDS(file = "singlets_PCA_D49.rds")
D28<- readRDS(file = "singlets_PCA_D28.rds")


Batch1 <- c(C90, D45, D46, D27, C95, D18, D17, D15)
anchors_Batch1 <- FindIntegrationAnchors(object.list = Batch1, dims = 1:30)
Batch1_integrated <- IntegrateData(anchorset = anchors_Batch1, dims = 1:30)
rm(C90, D45, D46, D27, C95, D18, D17, D15)

Batch2 <- c(C96, D13, D22, D50, C92, C94, D49, D28)
anchors_Batch2 <- FindIntegrationAnchors(object.list = Batch2, dims = 1:30)
Batch2_integrated <- IntegrateData(anchorset = anchors_Batch2, dims = 1:30)
rm(C96, D13, D22, D50, C92, C94, D49, D28)

TDI150 <- c(Batch1_integrated, Batch2_integrated)
anchors_TDI150 <- FindIntegrationAnchors(object.list = TDI150, dims = 1:30)
TDI150_integrated <- IntegrateData(anchorset = anchors_TDI150, dims = 1:30)
rm(Batch1_integrated, Batch2_integrated, Batch1, Batch2)


Idents(TDI150_integrated) <- "Sample_Name"
pdf("TDI150_bySampleName_QC.pdf", width=12, height=12)
VlnPlot(object = TDI150_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 4, pt.size=0, idents=NULL)
dev.off()

Idents(TDI150_integrated) <- "Genotype"
pdf("TDI150_byGenotype_QC.pdf", width=8, height=6)
VlnPlot(object = TDI150_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 2, pt.size=0, idents=NULL)
dev.off()

saveRDS(TDI150_integrated, file = "TDI150integrated.rds")

DefaultAssay(TDI150_integrated) <- 'integrated'

TDI150_integrated <- NormalizeData(TDI150_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
TDI150_integrated <- FindVariableFeatures(TDI150_integrated, selection.method = "vst", nfeatures = 3000)

TDI150_integrated <- ScaleData(TDI150_integrated, verbose = FALSE)
TDI150_integrated <- RunPCA(TDI150_integrated, features = VariableFeatures(object = TDI150_integrated), verbose = FALSE)

pdf("ElbowPlot_TDI150integrated.pdf", width=8, height=6)
ElbowPlot(TDI150_integrated)
dev.off()

TDI150_integrated <- FindNeighbors(TDI150_integrated, dims = 1:20)
TDI150_integrated <- FindClusters(TDI150_integrated, resolution = 0.5)
TDI150_integrated <- RunUMAP(TDI150_integrated, dims = 1: 20)
saveRDS(TDI150_integrated, file = 'PCA_0.5_TDI150integrated.rds')

str(TDI150_integrated)

#DefaultAssay(TDI150_integrated) <- 'RNA'
#TDI150_integrated <- NormalizeData(TDI150_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
#TDI150_integrated <- ScaleData(TDI150_integrated, features = rownames(TDI150_integrated))

Idents(TDI150_integrated) <- 'seurat_clusters'

pdf("UMAP_TDI150integrated.pdf", width=11, height=9)
DimPlot(TDI150_integrated, reduction = 'umap', label = F)
dev.off()

pdf("UMAP_byBatch_TDI150integrated.pdf", width=11, height=9)
DimPlot(TDI150_integrated, reduction = "umap", split.by = "Batch", label = F, ncol = 2)
dev.off()

pdf("UMAP_bySampleNam_TDI150integrated.pdf", width=12, height=15)
DimPlot(TDI150_integrated, reduction = "umap", split.by = "Sample_Name", label = F, ncol = 4)
dev.off()


DefaultAssay(TDI150_integrated) <- 'RNA'

TDI150_markers <- FindAllMarkers(TDI150_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1, test.use = "MAST")
write.csv(TDI150_markers, "ClusterMarkers_TDI150integrated.csv")


#Add marker genes
CellTypeMarkers <-c("Plp1", "Mbp", "Mobp","Slc17a7", "Nrgn", "Gad1", "Gad2", "Clu", "Aldoc", 
          "Pla2g7", "Cx3cr1", "P2ry12", "Csf1r","Scrg1", "Pdgfra", "Vtn", "Igfbp7",
          "Bnc2", "Slc47a1", "Ttr","Prox1","Mpped1","Map3k15","Cdh24")
markers.to.plot <- as.matrix(CellTypeMarkers)
p = DotPlot(object = TDI150_integrated, features = rev(x = markers.to.plot)) + RotatedAxis()

pdf("TDI150_annotation_combine_byBatch.pdf", width=10, height=10)
p+scale_colour_gradient2(low = "darkblue", mid = "white", high = "darkred")
dev.off()

clusterSummary = as.data.frame(table(Idents(TDI150_integrated), TDI150_integrated$Genotype))
write.csv(clusterSummary, 'Idents(ClusterSummary_byGenotype_TDI150integrated.csv')

clusterSummary = as.data.frame(table(Idents(TDI150_integrated), TDI150_integrated$Sample_Name))
write.csv(clusterSummary, 'Idents(ClusterSummary_bySampleName_TDI150integrated.csv')

