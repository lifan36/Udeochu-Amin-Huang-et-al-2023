#setwd("/athena/ganlab/scratch/jcu4001/JCUPs19cGas/TestRunJU1/TDI_sNucSeq150")
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(RColorBrewer)


TDI150_integrated <- readRDS(file = 'PCA_0.5_TDI150integrated.rds')

DefaultAssay(TDI150_integrated) <- 'RNA'
TDI150_integrated <- NormalizeData(TDI150_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
TDI150_integrated <- ScaleData(TDI150_integrated, features = rownames(TDI150_integrated))

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

saveRDS(TDI_integrated, file = "PostClustering_TDI150integrated.rds")