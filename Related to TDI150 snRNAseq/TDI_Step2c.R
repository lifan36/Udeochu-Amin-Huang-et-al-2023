setwd("/athena/ganlab/scratch/jcu4001/JCUPs19cGas/TestRunJU1/TDI_sNucSeq150")
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)

####################
#homotypic doublet proportion estimate
C90 <- readRDS(file = 'C90_Step1.rds')
annotations <- C90@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.031*6601) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
C90<- doubletFinder_v3(C90, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
#doublet finder with different classification stringencies

C90<- doubletFinder_v3(C90, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_205", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Elbow_C90.pdf", width=8, height=6)
ElbowPlot(C90,ndims=50)
dev.off()
C90<- FindNeighbors(object = C90, dims = 1:15)
C90<- FindClusters(object = C90, resolution = 0.5)
C90<- RunUMAP(C90, dims = 1:15)
#visualizing the singlet vs doublet cells
Idents(object = C90) <- "DF.classifications_0.25_0.005_205" 
pdf("UMAP_C90.pdf", width=8, height=6)
DimPlot(object = C90, reduction = 'umap', label = T)
dev.off()

saveRDS(C90,"after_doublet_detection_C90.rds")
#processing singlets ====
#remove doublets

singlets <- subset(C90, idents=c("Singlet"))
rm(C90)
saveRDS(singlets,"singlets_C90.rds")

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
pdf("UMAP_singlets_after_processing_C90.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"singlets_PCA_C90.rds")

######################
D45 <- readRDS(file = 'D45_Step1.rds')
annotations <- D45@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.039*8178) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
D45<- doubletFinder_v3(D45, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
D45<- doubletFinder_v3(D45, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_319", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Elbow_D45.pdf", width=8, height=6)
ElbowPlot(D45,ndims=50)
dev.off()
D45<- FindNeighbors(object = D45, dims = 1:15)
D45<- FindClusters(object = D45, resolution = 0.5)
D45<- RunUMAP(D45, dims = 1:15)
#visualizing the singlet vs doublet cells
Idents(oject = D45) <- "DF.classifications_0.25_0.005_319" 
pdf("UMAP_D45.pdf", width=8, height=6)
DimPlot(object = D45, reduction = 'umap', label = T)
dev.off()

saveRDS(D45,"after_doublet_detection_D45.rds")
#processing singlets ====
#remove doublets
singlets <- subset(D45, idents=c("Singlet"))
rm(D45)
saveRDS(singlets,"singlets_D45.rds")

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
pdf("UMAP_singlets_after_processing_D45.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"singlets_PCA_D45.rds")

##################
D46 <- readRDS(file = 'D46_Step1.rds')
annotations <- D46@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.031*7138) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
D46<- doubletFinder_v3(D46, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
D46<- doubletFinder_v3(D46, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_221", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Elbow_D46.pdf", width=8, height=6)
ElbowPlot(D46,ndims=50)
dev.off()
D46<- FindNeighbors(object = D46, dims = 1:15)
D46<- FindClusters(object = D46, resolution = 0.5)
D46<- RunUMAP(D46, dims = 1:15)
Idents(object = D46) <- "DF.classifications_0.25_0.005_221"
pdf("UMAP_D46.pdf", width=8, height=6)
DimPlot(object = D46, reduction = 'umap', label = T)
dev.off()

saveRDS(D46,"after_doublet_detection_D46.rds")
#processing singlets ====
#remove doublets
singlets <- subset(D46, idents=c("Singlet"))
rm(D46)
saveRDS(singlets,"singlets_D46.rds")

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
pdf("UMAP_singlets_after_processing_D46.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"singlets_PCA_D46.rds")

#################
D27 <- readRDS(file = 'D27_Step1.rds')
annotations <- D27@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.046*9400) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
D27<- doubletFinder_v3(D27, PCs=1:15, pN=0.25, pK=0.01, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
D27<- doubletFinder_v3(D27, PCs = 1:15, pN = 0.25, pK = 0.01, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.01_432", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Elbow_D27.pdf", width=8, height=6)
ElbowPlot(D27,ndims=50)
dev.off()
D27<- FindNeighbors(object = D27, dims = 1:15)
D27<- FindClusters(object = D27, resolution = 0.5)
D27<- RunUMAP(D27, dims = 1:15)
Idents(object = D27) <- "DF.classifications_0.25_0.01_432"
pdf("UMAP_D27.pdf", width=8, height=6)
DimPlot(object = D27, reduction = 'umap', label = T)
dev.off()

saveRDS(D27,"after_doublet_detection_D27.rds")
#processing singlets ====
#remove doublets
singlets <- subset(D27, idents=c("Singlet"))
rm(D27)
saveRDS(singlets,"singlets_D27.rds")

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
pdf("UMAP_singlets_after_processing_D27.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"singlets_PCA_D27.rds")

#############
C95 <- readRDS(file = 'C95_Step1.rds')
annotations <- C95@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.039*8587) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
C95<- doubletFinder_v3(C95, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
C95<- doubletFinder_v3(C95, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_335", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Elbow_C95.pdf", width=8, height=6)
ElbowPlot(C95,ndims=50)
dev.off()
C95<- FindNeighbors(object = C95, dims = 1:15)
C95<- FindClusters(object = C95, resolution = 0.5)
C95<- RunUMAP(C95, dims = 1:15)
Idents(object = C95) <- "DF.classifications_0.25_0.005_335"
pdf("UMAP_C95.pdf", width=8, height=6)
DimPlot(object = C95, reduction = 'umap', label = T)
dev.off()

saveRDS(C95,"after_doublet_detection_C95.rds")
#processing singlets ====
#remove doublets
singlets <- subset(C95, idents=c("Singlet"))
rm(C95)
saveRDS(singlets,"singlets_C95.rds")

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
pdf("UMAP_singlets_after_processing_C95.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"singlets_PCA_C95.rds")

#############
D18 <- readRDS(file = 'D18_Step1.rds')
annotations <- D18@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.031*6337) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
D18<- doubletFinder_v3(D18, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
D18 <- readRDS(file = 'D18_Step2b.rds')
D18<- doubletFinder_v3(D18, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_196", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Elbow_D18.pdf", width=8, height=6)
ElbowPlot(D18,ndims=50)
dev.off()
D18<- FindNeighbors(object = D18, dims = 1:15)
D18<- FindClusters(object = D18, resolution = 0.5)
D18<- RunUMAP(D18, dims = 1:15)
Idents(object = D18) <- "DF.classifications_0.25_0.005_196"
pdf("UMAP_D18.pdf", width=8, height=6)
DimPlot(object = D18, reduction = 'umap', label = T)
dev.off()

saveRDS(D18,"after_doublet_detection_D18.rds")
#processing singlets ====
#remove doublets
singlets <- subset(D18, idents=c("Singlet"))
rm(D18)
saveRDS(singlets,"singlets_D18.rds")

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
pdf("UMAP_singlets_after_processing_D18.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"singlets_PCA_D18.rds")

##############
D17 <- readRDS(file = 'D17_Step1.rds')
annotations <- D17@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.039*7631) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
D17<- doubletFinder_v3(D17, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
D17<- doubletFinder_v3(D17, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_298", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Elbow_D17.pdf", width=8, height=6)
ElbowPlot(D17,ndims=50)
dev.off()
D17<- FindNeighbors(object = D17, dims = 1:15)
D17<- FindClusters(object = D17, resolution = 0.5)
D17<- RunUMAP(D17, dims = 1:15)
Idents(object = D17) <- "DF.classifications_0.25_0.005_298"
pdf("UMAP_D17.pdf", width=8, height=6)
DimPlot(object = D17, reduction = 'umap', label = T)
dev.off()

saveRDS(D17,"after_doublet_detection_D17.rds")
#processing singlets ====
#remove doublets
singlets <- subset(D17, idents=c("Singlet"))
rm(D17)
saveRDS(singlets,"singlets_D17.rds")

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
pdf("UMAP_singlets_after_processing_D17.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"singlets_PCA_D17.rds")

################
D15 <- readRDS(file = 'D15_Step1.rds')
annotations <- D15@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.046*9428) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
D15<- doubletFinder_v3(D15, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
D15<- doubletFinder_v3(D15, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_434", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Elbow_D15.pdf", width=8, height=6)
ElbowPlot(D15,ndims=50)
dev.off()
D15<- FindNeighbors(object = D15, dims = 1:15)
D15<- FindClusters(object = D15, resolution = 0.5)
D15<- RunUMAP(D15, dims = 1:15)
Idents(object = D15) <- "DF.classifications_0.25_0.005_434"
pdf("UMAP_D15.pdf", width=8, height=6)
DimPlot(object = D15, reduction = 'umap', label = T)
dev.off()

saveRDS(D15,"after_doublet_detection_D15.rds")
#processing singlets ====
#remove doublets
singlets <- subset(D15, idents=c("Singlet"))
rm(D15)
saveRDS(singlets,"singlets_D15.rds")

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
pdf("UMAP_singlets_after_processing_D15.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"singlets_PCA_D15.rds")

################
C96 <- readRDS(file = 'C96_Step1.rds')
annotations <- C96@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.030*8402) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
C96<- doubletFinder_v3(C96, PCs=1:15, pN=0.25, pK=0.02, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
C96<- doubletFinder_v3(C96, PCs = 1:15, pN = 0.25, pK = 0.02, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.02_252", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Elbow_C96.pdf", width=8, height=6)
ElbowPlot(C96,ndims=50)
dev.off()
C96<- FindNeighbors(object = C96, dims = 1:15)
C96<- FindClusters(object = C96, resolution = 0.5)
C96<- RunUMAP(C96, dims = 1:15)
Idents(object = C96) <- "DF.classifications_0.25_0.02_252"
pdf("UMAP_C96.pdf", width=8, height=6)
DimPlot(object = C96, reduction = 'umap', label = T)
dev.off()

saveRDS(C96,"after_doublet_detection_C96.rds")
#processing singlets ====
#remove doublets
singlets <- subset(C96, idents=c("Singlet"))
rm(C96)
saveRDS(singlets,"singlets_C96.rds")

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
pdf("UMAP_singlets_after_processing_C96.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"singlets_PCA_C96.rds")
########################
D13 <- readRDS(file = 'D13_Step1.rds')
annotations <- D13@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.031*6280) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
D13<- doubletFinder_v3(D13, PCs=1:15, pN=0.25, pK=0.01, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
D13<- doubletFinder_v3(D13, PCs = 1:15, pN = 0.25, pK = 0.01, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.01_195", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Elbow_D13.pdf", width=8, height=6)
ElbowPlot(D13,ndims=50)
dev.off()
D13<- FindNeighbors(object = D13, dims = 1:15)
D13<- FindClusters(object = D13, resolution = 0.5)
D13<- RunUMAP(D13, dims = 1:15)
Idents(object = D13) <- "DF.classifications_0.25_0.01_195"
pdf("UMAP_D13.pdf", width=8, height=6)
DimPlot(object = D13, reduction = 'umap', label = T)
dev.off()

saveRDS(D13,"after_doublet_detection_D13.rds")
#processing singlets ====
#remove doublets
singlets <- subset(D13, idents=c("Singlet"))
rm(D13)
saveRDS(singlets,"singlets_D13.rds")

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
pdf("UMAP_singlets_after_processing_D13.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"singlets_PCA_D13.rds")
###################
D22 <- readRDS(file = 'D22_Step1.rds')
annotations <- D22@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.031*6280) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
D22<- doubletFinder_v3(D22, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
D22<- doubletFinder_v3(D22, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_195", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Elbow_D22.pdf", width=8, height=6)
ElbowPlot(D22,ndims=50)
dev.off()
D22<- FindNeighbors(object = D22, dims = 1:15)
D22<- FindClusters(object = D22, resolution = 0.5)
D22<- RunUMAP(D22, dims = 1:15)
Idents(object = D22) <- "DF.classifications_0.25_0.005_195"
pdf("UMAP_D22.pdf", width=8, height=6)
DimPlot(object = D22, reduction = 'umap', label = T)
dev.off()

saveRDS(D22,"after_doublet_detection_D22.rds")
#processing singlets ====
#remove doublets
singlets <- subset(D22, idents=c("Singlet"))
rm(D22)
saveRDS(singlets,"singlets_D22.rds")

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
pdf("UMAP_singlets_after_processing_D22.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"singlets_PCA_D22.rds")
##################
D50 <- readRDS(file = 'D50_Step1.rds')
annotations <- D50@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.046*9705) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
D50<- doubletFinder_v3(D50, PCs=1:15, pN=0.25, pK=0.01, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
D50<- doubletFinder_v3(D50, PCs = 1:15, pN = 0.25, pK = 0.01, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.01_446", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Elbow_D50.pdf", width=8, height=6)
ElbowPlot(D50,ndims=50)
dev.off()
D50<- FindNeighbors(object = D50, dims = 1:15)
D50<- FindClusters(object = D50, resolution = 0.5)
D50<- RunUMAP(D50, dims = 1:15)
Idents(object = D50) <- "DF.classifications_0.25_0.01_446"
pdf("UMAP_D50.pdf", width=8, height=6)
DimPlot(object = D50, reduction = 'umap', label = T)
dev.off()

saveRDS(D50,"after_doublet_detection_D50.rds")
#processing singlets ====
#remove doublets
singlets <- subset(D50, idents=c("Singlet"))
rm(D50)
saveRDS(singlets,"singlets_D50.rds")

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
pdf("UMAP_singlets_after_processing_D50.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"singlets_PCA_D50.rds")
#############
C92 <- readRDS(file = 'C92_Step1.rds')
annotations <- C92@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*12132) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
C92<- doubletFinder_v3(C92, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
C92<- doubletFinder_v3(C92, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_740", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Elbow_C92.pdf", width=8, height=6)
ElbowPlot(C92,ndims=50)
dev.off()
C92<- FindNeighbors(object = C92, dims = 1:15)
C92<- FindClusters(object = C92, resolution = 0.5)
C92<- RunUMAP(C92, dims = 1:15)
Idents(object = C92) <- "DF.classifications_0.25_0.005_740"
pdf("UMAP_C92.pdf", width=8, height=6)
DimPlot(object = C92, reduction = 'umap', label = T)
dev.off()

saveRDS(C92,"after_doublet_detection_C92.rds")
#processing singlets ====
#remove doublets
singlets <- subset(C92, idents=c("Singlet"))
rm(C92)
saveRDS(singlets,"singlets_C92.rds")

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
pdf("UMAP_singlets_after_processing_C92.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"singlets_PCA_C92.rds")
##############
C94 <- readRDS(file = 'C94_Step1.rds')
annotations <- C94@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.039*8393) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
C94<- doubletFinder_v3(C94, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
C94<- doubletFinder_v3(C94, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_327", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Elbow_C94.pdf", width=8, height=6)
ElbowPlot(C94,ndims=50)
dev.off()
C94<- FindNeighbors(object = C94, dims = 1:15)
C94<- FindClusters(object = C94, resolution = 0.5)
C94<- RunUMAP(C94, dims = 1:15)
Idents(object = C94) <- "DF.classifications_0.25_0.005_327"
pdf("UMAP_C94.pdf", width=8, height=6)
DimPlot(object = C94, reduction = 'umap', label = T)
dev.off()

saveRDS(C94,"after_doublet_detection_C94.rds")
#processing singlets ====
#remove doublets
singlets <- subset(C94, idents=c("Singlet"))
rm(C94)
saveRDS(singlets,"singlets_C94.rds")

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
pdf("UMAP_singlets_after_processing_C94.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"singlets_PCA_C94.rds")
##################
D49 <- readRDS(file = 'D49_Step1.rds')
annotations <- D49@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.046*8931) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
D49<- doubletFinder_v3(D49, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
D49<- doubletFinder_v3(D49, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_411", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("Elbow_D49.pdf", width=8, height=6)
ElbowPlot(D49,ndims=50)
dev.off()
D49<- FindNeighbors(object = D49, dims = 1:15)
D49<- FindClusters(object = D49, resolution = 0.5)
D49<- RunUMAP(D49, dims = 1:15)
Idents(object = D49) <- "DF.classifications_0.25_0.005_411"
pdf("UMAP_D49.pdf", width=8, height=6)
DimPlot(object = D49, reduction = 'umap', label = T)
dev.off()

saveRDS(D49,"after_doublet_detection_D49.rds")
#processing singlets ====
#remove doublets
singlets <- subset(D49, idents=c("Singlet"))
rm(D49)
saveRDS(singlets,"singlets_D49.rds")

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
pdf("UMAP_singlets_after_processing_D49.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"singlets_PCA_D49.rds")
###############
D28 <- readRDS(file = 'D28_Step1.rds')
annotations <- D28@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.046*10216) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
D28<- doubletFinder_v3(D28, PCs=1:15, pN=0.25, pK=0.011, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
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