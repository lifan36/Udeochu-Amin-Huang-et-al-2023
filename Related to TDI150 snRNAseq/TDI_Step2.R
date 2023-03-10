setwd("/athena/ganlab/scratch/jcu4001/JCUPs19cGas/TestRunJU1/TDI_sNucSeq150")
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)

C90 <- readRDS(file = 'C90_Step1.rds')
a <- length(C90@meta.data$seurat_clusters)

D45 <- readRDS(file = 'D45_Step1.rds')
b <- length(D45@meta.data$seurat_clusters)

D46 <- readRDS(file = 'D46_Step1.rds')
c <- length(D46@meta.data$seurat_clusters)

D27 <- readRDS(file = 'D27_Step1.rds')
d <- length(D27@meta.data$seurat_clusters)

C95 <- readRDS(file = 'C95_Step1.rds')
e <- length(C95@meta.data$seurat_clusters)

D18 <- readRDS(file = 'D18_Step1.rds')
f <- length(D18@meta.data$seurat_clusters)

D17 <- readRDS(file = 'D17_Step1.rds')
g <- length(D17@meta.data$seurat_clusters)

D15 <- readRDS(file = 'D15_Step1.rds')
h <- length(D15@meta.data$seurat_clusters)

C96 <- readRDS(file = 'C96_Step1.rds')
i <- length(C96@meta.data$seurat_clusters)

D13 <- readRDS(file = 'D13_Step1.rds')
j <- length(D13@meta.data$seurat_clusters)

D22 <- readRDS(file = 'D22_Step1.rds')
k <- length(D22@meta.data$seurat_clusters)

D50 <- readRDS(file = 'D50_Step1.rds')
l <- length(D50@meta.data$seurat_clusters)

C92 <- readRDS(file = 'C92_Step1.rds')
m <- length(C92@meta.data$seurat_clusters)

C94 <- readRDS(file = 'C94_Step1.rds')
n <- length(C94@meta.data$seurat_clusters)

D49 <- readRDS(file = 'D49_Step1.rds')
o <- length(D49@meta.data$seurat_clusters)

D28 <- readRDS(file = 'D28_Step1.rds')
p <- length(D28@meta.data$seurat_clusters)

length_TDI150Cohort <- cbind(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p)
write.csv(length_TDI150Cohort, 'length_TDI150Cohort.csv')

