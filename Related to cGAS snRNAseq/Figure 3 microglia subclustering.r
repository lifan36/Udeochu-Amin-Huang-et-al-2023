#install.packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(reshape2)
library(MAST)
library(RColorBrewer)


setwd("/Users/lifan/Desktop/data_analysis/Joe_cGAS/Final_script_with_Bang_20230223/no_68")
Cluster_MG <- readRDS(file = 'MG_cGAS_integrated_Annotated.rds')
# "Cluster_MGreclusetered.rds" contains Seurat object for MG cluster as is
Idents(Cluster_MG) <- 'Genotype'
Cluster_MG <- subset(Cluster_MG, idents = c('P301S tau: +; cGAS: -/-', 'P301S tau: +; cGAS: +/-',
                                            'P301S tau: +; cGAS: +/+', 'P301S tau: -; cGAS: +/+'))

Cluster_MG <- FindVariableFeatures(object = Cluster_MG)
Cluster_MG <- RunPCA(Cluster_MG, features = VariableFeatures(object = Cluster_MG))
DimPlot(Cluster_MG, reduction = "pca")
ElbowPlot(Cluster_MG)

Cluster_MG <- FindNeighbors(Cluster_MG, dims = 1:7)
Cluster_MG <- FindClusters(Cluster_MG, resolution = 0.5)
Cluster_MG <- RunUMAP(Cluster_MG, dims = 1:7)

Idents(Cluster_MG) <- 'seurat_clusters'
#pdf("UMAP_Microglia.pdf", width=5, height=4)
DimPlot(Cluster_MG, reduction = "umap", label = T, split.by = "Genotype", ncol = 2)
#dev.off()

#find markers then regroup and rename clusters
#cluster 5, 6, 7 do not show expression of microglial genes

Cluster_MG <- subset(Cluster_MG, idents = c('0', '1', '2', '3', '4'))
Cluster_MG <- RenameIdents(Cluster_MG, '0' = "1", '3' = "1", '1' = "4", '4' = "3", '2' = "2")
levels(Cluster_MG) <- c("1", "2", "3", "4")

pdf("UMAP_Microglia1.pdf", width=5, height=4)
DimPlot(Cluster_MG, reduction = "umap", label = T, split.by = "Genotype", ncol = 2)
dev.off()

saveRDS(Cluster_MG, file = 'Cluster_MGreclusetered_03072023.rds')

pdf("VlnPlot_Microglial_supp_figure.pdf", width=14, height=2.5)
VlnPlot(Cluster_MG, features = c("P2ry12", "Siglech", "Apoe", "Itgax","Stat1","Parp14","Trim30a","Rnf213"), pt.size = 0, ncol = 8)
dev.off()


#Figure 4b
pdf("FeaturePlot_Microglial.pdf", width=6, height=5)
FeaturePlot(Cluster_MG, features = c("P2ry12", "Siglech", "Sall1", "Csf1r"), dims = c(1, 2), cols = colorRampPalette((brewer.pal(9, "Blues")) )(255))
dev.off()
#Figure 4c
DimPlot(Cluster_MG, reduction = "umap", label = F, split.by = 'Genotype', ncol = 2)

#Figure 4e
pdf("Interferon_genes_in_MG_only_DotPlot_no68.pdf", width=8, height=3)
p = DotPlot(Cluster_MG, features = c("Trim30d", "Trim30a", "Stat1","Sp100", "Slfn8", "Rnf213", "Parp14",
                                     "Nlrc5", "Herc6", "Eif2ak2", "Ddx60","Gm4951", "Apobec3","Cd300lf","Akt3"), assay = 'RNA') + RotatedAxis()
p+scale_colour_gradient2(low = "darkblue", mid = "white", high = "darkred")
dev.off()

Cluster_MG_Summary = as.data.frame(table(Idents(Cluster_MG), Cluster_MG$Sample_Name))
write.csv(Cluster_MG_Summary, 'ClusterSummary_bySample_Name_Microglia.csv')

Cluster_MG_Summary = as.data.frame(table(Idents(Cluster_MG), Cluster_MG$Genotype))
write.csv(Cluster_MG_Summary, 'ClusterSummary_byGenotype_Microglia.csv')

MGclusters.markers <- FindAllMarkers(object = Cluster_MG, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(MGclusters.markers,"ClusterMarkers_Microglia.csv")


Idents(Cluster_MG) <- 'Genotype'
PS19KOvsPS19WT_DE <- FindMarkers(Cluster_MG, ident.1 = "P301S tau: +; cGAS: -/-", ident.2 = "P301S tau: +; cGAS: +/+", logfc.threshold = 0,
                                 test.use = "MAST", assay ='RNA')
write.csv(PS19KOvsPS19WT_DE, file = 'PS19KOvsPS19WT_DE_MAST_RNAassay_Microglia.csv')

PS19HetvsPS19WT_DE <- FindMarkers(Cluster_MG, ident.1 = "P301S tau: +; cGAS: +/-", ident.2 = "P301S tau: +; cGAS: +/+", logfc.threshold = 0,
                                  test.use = "MAST", assay ='RNA')
write.csv(PS19HetvsPS19WT_DE, file = 'PS19HetvsPS19WT_DE_MAST_RNAassay_Microglia.csv')

PS19KOvsPS19Het_DE <- FindMarkers(Cluster_MG, ident.1 = "P301S tau: +; cGAS: -/-", ident.2 = "P301S tau: +; cGAS: +/-", logfc.threshold = 0,
                                  test.use = "MAST", assay ='RNA')
write.csv(PS19KOvsPS19Het_DE, file = 'PS19KOvsPS19Het_DE_MAST_RNAassay_Microglia.csv')

KOvsWT_DE <- FindMarkers(Cluster_MG, ident.1 = "P301S tau: -; cGAS: -/-", ident.2 = "P301S tau: -; cGAS: +/+", logfc.threshold = 0,
                         test.use = "MAST", assay ='RNA')
write.csv(KOvsWT_DE, file = 'KOvsWT_DE_MAST_RNAassay_Microglia.csv')

HetvsWT_DE <- FindMarkers(Cluster_MG, ident.1 = "P301S tau: -; cGAS: +/-", ident.2 = "P301S tau: -; cGAS: +/+", logfc.threshold = 0,
                          test.use = "MAST", assay ='RNA')
write.csv(HetvsWT_DE, file = 'HetvsWT_DE_MAST_RNAassay_Microglia.csv')

PS19KOvsKO_DE <- FindMarkers(Cluster_MG, ident.1 = "P301S tau: +; cGAS: -/-", ident.2 = "P301S tau: -; cGAS: -/-", logfc.threshold = 0,
                             test.use = "MAST", assay ='RNA')
write.csv(PS19KOvsKO_DE, file = 'PS19KOvsKO_DE_MAST_RNAassay_Microglia.csv')

PS19WTvsWT_DE <- FindMarkers(Cluster_MG, ident.1 = "P301S tau: +; cGAS: +/+", ident.2 = "P301S tau: -; cGAS: +/+", logfc.threshold = 0,
                             test.use = "MAST", assay ='RNA')
write.csv(PS19WTvsWT_DE, file = 'PS19WTvsWT_DE_MAST_RNAassay_Microglia.csv')

PS19HetvsHet_DE <- FindMarkers(Cluster_MG, ident.1 = "P301S tau: +; cGAS: +/-", ident.2 = "P301S tau: -; cGAS: +/-", logfc.threshold = 0,
                               test.use = "MAST", assay ='RNA')
write.csv(PS19HetvsHet_DE, file = 'PS19HetvsHet_DE_MAST_RNAassay_Microglia.csv')

top10 <- MGclusters.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
DoHeatmap(Cluster_MG, features = top10$gene) + scale_fill_gradientn(colors = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)))
p = DotPlot(Cluster_MG, features = top10$gene, assay = 'RNA') + RotatedAxis()
p+scale_colour_gradient2(low = "darkblue", mid = "white", high = "darkred")

