setwd("/athena/ganlab/scratch/jcu4001/JCUPs19cGas/TestRunJU1/TDI_sNucSeq150")
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)

###############
D28 <- readRDS(file = 'D28_Step1.rds')
annotations <- D28@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.046*10216) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
D28<- doubletFinder_v3(D28, PCs=1:15, pN=0.25, pK=0.01, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
D28<- doubletFinder_v3(D28, PCs = 1:15, pN = 0.25, pK = 0.01, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.01_470", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Elbow_D28.pdf", width=8, height=6)
ElbowPlot(D28,ndims=50)
dev.off()
D28<- FindNeighbors(object = D28, dims = 1:15)
D28<- FindClusters(object = D28, resolution = 0.5)
D28<- RunUMAP(D28, dims = 1:15)
Idents(object = D28) <- "DF.classifications_0.25_0.01_470"
pdf("UMAP_D28.pdf", width=8, height=6)
DimPlot(object = D28, reduction = 'umap', label = T)
dev.off()

saveRDS(D28,"after_doublet_detection_D28.rds")
#processing singlets ====
#remove doublets
singlets <- subset(D28, idents=c("Singlet"))
rm(D28)
saveRDS(singlets,"singlets_D28.rds")

Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture

singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.5)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("UMAP_singlets_after_processing_D28.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"singlets_PCA_D28.rds")

##############