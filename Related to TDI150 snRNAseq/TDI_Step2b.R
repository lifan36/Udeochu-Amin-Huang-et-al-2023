setwd("/athena/ganlab/scratch/jcu4001/JCUPs19cGas/TestRunJU1/TDI_sNucSeq150")
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)

C90 <- readRDS(file = 'C90_Step1.rds')
#a <- length(C90@meta.data$seurat_clusters)

D45 <- readRDS(file = 'D45_Step1.rds')
#b <- length(D45@meta.data$seurat_clusters)

D46 <- readRDS(file = 'D46_Step1.rds')
#c <- length(D46@meta.data$seurat_clusters)

D27 <- readRDS(file = 'D27_Step1.rds')
#d <- length(D27@meta.data$seurat_clusters)

C95 <- readRDS(file = 'C95_Step1.rds')
#e <- length(C95@meta.data$seurat_clusters)

D18 <- readRDS(file = 'D18_Step1.rds')
#f <- length(D18@meta.data$seurat_clusters)

D17 <- readRDS(file = 'D17_Step1.rds')
#g <- length(D17@meta.data$seurat_clusters)

D15 <- readRDS(file = 'D15_Step1.rds')
#h <- length(D15@meta.data$seurat_clusters)

C96 <- readRDS(file = 'C96_Step1.rds')
#i <- length(C96@meta.data$seurat_clusters)

D13 <- readRDS(file = 'D13_Step1.rds')
#j <- length(D13@meta.data$seurat_clusters)

D22 <- readRDS(file = 'D22_Step1.rds')
#k <- length(D22@meta.data$seurat_clusters)

D50 <- readRDS(file = 'D50_Step1.rds')
#l <- length(D50@meta.data$seurat_clusters)

C92 <- readRDS(file = 'C92_Step1.rds')
#m <- length(C92@meta.data$seurat_clusters)

C94 <- readRDS(file = 'C94_Step1.rds')
#n <- length(C94@meta.data$seurat_clusters)

D49 <- readRDS(file = 'D49_Step1.rds')
#o <- length(D49@meta.data$seurat_clusters)

D28 <- readRDS(file = 'D28_Step1.rds')
#p <- length(D28@meta.data$seurat_clusters)

#length_TDI150Cohort <- cbind(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p)
#write.csv(length_TDI150Cohort, 'length_TDI150Cohort.csv')

####################
#homotypic doublet proportion estimate
annotations <- C90@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.031*6601) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
C90<- doubletFinder_v3(C90, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
saveRDS(C90, 'C90_Step2b.rds')

annotations <- D45@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.039*8178) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
D45<- doubletFinder_v3(D45, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
saveRDS(D45, 'D45_Step2b.rds')

annotations <- D46@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.031*7138) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
D46<- doubletFinder_v3(D46, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
saveRDS(D46, 'D46_Step2b.rds')

annotations <- D27@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.046*9400) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
D27<- doubletFinder_v3(D27, PCs=1:15, pN=0.25, pK=0.01, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
saveRDS(D27, 'D27_Step2b.rds')

annotations <- C95@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.039*8587) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
C95<- doubletFinder_v3(C95, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
saveRDS(C95, 'C95_Step2b.rds')

annotations <- D18@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.031*6337) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
D18<- doubletFinder_v3(D18, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
saveRDS(D18, 'D18_Step2b.rds')

annotations <- D17@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.039*7631) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
D17<- doubletFinder_v3(D17, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
saveRDS(D17, 'D17_Step2b.rds')

annotations <- D15@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.046*9428) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
D15<- doubletFinder_v3(D15, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
saveRDS(D15, 'D15_Step2b.rds')

annotations <- C96@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.030*8402) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
C96<- doubletFinder_v3(C96, PCs=1:15, pN=0.25, pK=0.02, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
saveRDS(C96, 'C96_Step2b.rds')

annotations <- D13@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.031*6280) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
D13<- doubletFinder_v3(D13, PCs=1:15, pN=0.25, pK=0.01, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
saveRDS(D13, 'D13_Step2b.rds')

annotations <- D22@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.031*6280) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
D22<- doubletFinder_v3(D22, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
saveRDS(D22, 'D22_Step2b.rds')

annotations <- D50@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.046*9705) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
D50<- doubletFinder_v3(D50, PCs=1:15, pN=0.25, pK=0.01, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
saveRDS(D50, 'D50_Step2b.rds')

annotations <- C92@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*12132) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
C92<- doubletFinder_v3(C92, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
saveRDS(C92, 'C92_Step2b.rds')

annotations <- C94@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.039*8393) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
C94<- doubletFinder_v3(C94, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
saveRDS(C94, 'C94_Step2b.rds')

annotations <- D49@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.046*8931) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
D49<- doubletFinder_v3(D49, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
saveRDS(D49, 'D49_Step2b.rds')

annotations <- D28@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.046*10216) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
D28<- doubletFinder_v3(D28, PCs=1:15, pN=0.25, pK=0.011, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
saveRDS(D28, 'D28_Step2b.rds')