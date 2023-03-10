#setwd("/athena/ganlab/scratch/jcu4001/JCUPs19cGas/TestRunJU1/TDI_sNucSeq150")
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(RColorBrewer)

TDI150_integrated <- readRDS(file = "PostClustering_TDI150integrated.rds")
Idents(TDI150_integrated) <- 'seurat_clusters'

TDI150_integrated = subset(TDI150_integrated, idents = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                                                         "11", "12", "13", "14", "15", "17", "18"))

saveRDS(TDI150_integrated, file = "TDI150_integrated_Trimmed.rds")

pdf("UMAP_by_SampleName_TDI150_integratedTrimmed.pdf", width=12, height=12)
DimPlot(TDI150_integrated, reduction = "umap", split.by = "Sample_Name", label = F, ncol = 4)
dev.off()
pdf("UMAP_byCondition1_TDI150_integratedTrimmed.pdf", width=15, height=10)
DimPlot(TDI150_integrated, reduction = "umap", split.by = "Condition_1", label = T, ncol = 3)
dev.off()

clusterSummary = as.data.frame(table(Idents(TDI150_integrated), TDI150_integrated$Sample_Name))
write.csv(clusterSummary, 'ClusterSummary_bySampleName_TDI150_integratedTrimmed.csv')

CellTypeMarkers <-c("Plp1", "Mbp", "Mobp","Slc17a7", "Nrgn", "Gad1", "Gad2", "Clu", "Aldoc", 
                    "Pla2g7", "Cx3cr1", "P2ry12", "Csf1r","Scrg1", "Pdgfra", "Vtn", "Igfbp7",
                    "Bnc2", "Slc47a1", "Ttr")
markers.to.plot <- as.matrix(CellTypeMarkers)
p = DotPlot(object = TDI150_integrated, features = rev(x = markers.to.plot)) + RotatedAxis()

pdf("CellTypeMarkers_TDI150_integratedTrimmed.pdf", width=10, height=10)
p+scale_colour_gradient2(low = "darkblue", mid = "white", high = "darkred")
dev.off()

TDI150_integrated_Annotated <- RenameIdents(TDI150_integrated,
                                            `0` = "Excitatory Neuron", `1`="Oligodendrocyte", `2`="Astrocyte", `3`="Inhibitory Neuron",
                                            `4`="Excitatory Neuron", `5`="Excitatory Neuron", `6`="Excitatory Neuron", `7`="Microglia",
                                            `8`="Excitatory Neuron", `9`="OPC", `10`="Inhibitory Neuron", `11`="Excitatory Neuron",
                                            `12`="Inhibitory Neuron", `13`="Ependymal Cell", `14`="Excitatory Neuron", `15`="Vascular Cell",
                                            `17`="Endothelia", `18`="Endothelia")

saveRDS(TDI150_integrated_Annotated, file = "TDI150_integrated_Annotated.rds")

Idents(TDI150_integrated_Annotated) <- 'seurat_clusters'

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
Cluster_EN <- subset(TDI150_integrated_Annotated, idents = "excitatory neurons")
Cluster_IN <- subset(TDI150_integrated_Annotated, idents = "inhibitory neurons")
Cluster_MG <- subset(TDI150_integrated_Annotated, idents = "microglia")
Cluster_AST <- subset(TDI150_integrated_Annotated, idents = "astrocytes")
Cluster_OL <- subset(TDI150_integrated_Annotated, idents = "oligodendrocytes")
Cluster_OPC <- subset(TDI150_integrated_Annotated, idents = "OPCs")

saveRDS(Cluster_EN, file = "EN_TDI150_integrated_Annotated.rds")
saveRDS(Cluster_IN, file = "IN_TDI150_integrated_Annotated.rds")
saveRDS(Cluster_MG, file = "MG_TDI150_integrated_Annotated.rds")
saveRDS(Cluster_AST, file = "Astro_TDI150_integrated_Annotated.rds")
saveRDS(Cluster_OL, file = "Oligo_TDI150_integrated_Annotated.rds")
saveRDS(Cluster_OPC, file = "OPC_TDI150_integrated_Annotated.rds")

