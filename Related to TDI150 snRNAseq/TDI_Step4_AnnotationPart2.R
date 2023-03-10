#setwd("/athena/ganlab/scratch/jcu4001/JCUPs19cGas/TestRunJU1/TDI_sNucSeq150")
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(RColorBrewer)

TDI150_integrated_Annotated <- readRDS(file = "TDI150_integrated_Annotated.rds")

pdf("UMAP_TDI150_integratedAnnotated.pdf", width=10, height=8)
DimPlot(TDI150_integrated_Annotated, reduction = 'umap', label = F)
dev.off()

pdf("UMAP_bySampleName_TDI150_integratedAnnotated.pdf", width=16, height=10)
DimPlot(TDI150_integrated_Annotated, reduction = 'umap', split.by = "Sample_Name", label = F, ncol = 4)
dev.off()

pdf("UMAP_byCondition1_TDI150_integratedAnnotated.pdf", width=10, height=10)
DimPlot(TDI150_integrated_Annotated, reduction = 'umap', split.by = "Condition_1", label = F, ncol = 2)
dev.off()

#Subset celltypes
Cluster_EN <- subset(TDI150_integrated_Annotated, idents = "Excitatory Neuron")
Cluster_IN <- subset(TDI150_integrated_Annotated, idents = "Inhibitory Neuron")
Cluster_MG <- subset(TDI150_integrated_Annotated, idents = "Microglia")
Cluster_AST <- subset(TDI150_integrated_Annotated, idents = "Astrocyte")
Cluster_OL <- subset(TDI150_integrated_Annotated, idents = "Oligodendrocyte")
Cluster_OPC <- subset(TDI150_integrated_Annotated, idents = "OPC")

saveRDS(Cluster_EN, file = "EN_TDI150_integrated_Annotated.rds")
saveRDS(Cluster_IN, file = "IN_TDI150_integrated_Annotated.rds")
saveRDS(Cluster_MG, file = "MG_TDI150_integrated_Annotated.rds")
saveRDS(Cluster_AST, file = "Astro_TDI150_integrated_Annotated.rds")
saveRDS(Cluster_OL, file = "Oligo_TDI150_integrated_Annotated.rds")
saveRDS(Cluster_OPC, file = "OPC_TDI150_integrated_Annotated.rds")

