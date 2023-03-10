#install.packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(reshape2)

#load in data from Cell Ranger or other counts data ====
setwd("/athena/ganlab/scratch/lif4001/Joe_cGas/final")
cGAS_integrated = readRDS(file = 'cGAS_integrated_PCA_0.1_byBatch.rds')
Idents(cGAS_integrated) <- 'seurat_clusters'

cGAS_integrated = subset(cGAS_integrated, idents = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                                                     "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"))

pdf("UMAP_by_SampleName_cGAS_integrated.pdf", width=12, height=15)
DimPlot(cGAS_integrated, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 4)
dev.off()
pdf("UMAP_byGenotype_cGAS_integrated.pdf", width=15, height=10)
DimPlot(cGAS_integrated, reduction = "umap", split.by = "Genotype", label = T, ncol = 3)
dev.off()

clusterSummary = as.data.frame(table(Idents(cGAS_integrated), cGAS_integrated$Genotype))
write.csv(clusterSummary, 'Idents(cGAS_integrated_ClusterSummary_byGenotype.csv')

clusterSummary = as.data.frame(table(Idents(cGAS_integrated), cGAS_integrated$orig.ident))
write.csv(clusterSummary, 'Idents(cGAS_integrated_ClusterSummary_byOrigIdent.csv')

saveRDS(cGAS_integrated, file = "cGAS_integrated_Trimmed.rds")

Idents(cGAS_integrated) <- "seurat_clusters"
cGAS_integrated <- RenameIdents(cGAS_integrated,
                                `0` = "oligodendrocytes", `1`="astrocytes", `2`="excitatory neurons", `3`="excitatory neurons",
                                `4`="ambiguous",`5`="excitatory neurons", `6`="microglia", `7`="excitatory neurons",
                                `8`="excitatory neurons", `9`="inhibitory neurons", `10`="OPCs", `11`="excitatory neurons",
                                `12`="inhibitory neurons", `13`="excitatory neurons", `14`="excitatory neurons", `15`="excitatory neurons",
                                `16`="excitatory neurons", `17`="excitatory neurons", `18`="ependymal cells", `19`="vascular cells",
                                `20`="endothelial cells")
cGAS_integrated$celltype <- Idents(cGAS_integrated)
saveRDS(cGAS_integrated, file = "cGAS_integrated_Annotated.rds")

data <- cGAS_integrated
# calculate ratio of each genotype in each cell type cluster
a<-as.data.frame(table(data$Genotype,data$celltype))
colnames(a)<-c("clusters","cell.type","cell.no")
agg<-aggregate(cell.no~clusters,a,sum)
a$cluster.total <- agg$cell.no[match(a$clusters,agg$clusters)]
a$ratio<-a$cell.no/a$cluster.total

ggplot(a,aes(x=clusters, y=ratio, fill=cell.type))+
  geom_bar(stat="identity")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("Genotype")+
  ylab("Cell type ratio per genotype") + RotatedAxis()

ggsave("genotype_celltype_distribution.pdf",plot=last_plot(),
       width=4,height=4,units="in")


pdf("UMAP_cGAS_integrated_Annotated.pdf", width=12, height=8)
DimPlot(cGAS_integrated, reduction = 'umap', label = TRUE)
dev.off()
pdf("UMAP_byGenotype_cGAS_integrated.pdf", width=16, height=10)
DimPlot(cGAS_integrated, reduction = 'umap', split.by = "Genotype", label = TRUE, ncol = 3)
dev.off()

markers.to.plot <- c("Plp1", "Mbp", "Mobp","Clu", "Aldoc","Pla2g7","Slc17a7","Nrgn",
                     "Cx3cr1", "P2ry12", "Csf1r", "Gad1", "Gad2",  
                     "Scrg1", "Pdgfra", "Ttr", "Bnc2", "Slc47a1", "Vtn", "Igfbp7"
                     )
pdf("DotPlot_cGAS_integrated.pdf", width=7, height=4.5)
DotPlot(object = cGAS_integrated, features = rev(x = markers.to.plot)) + RotatedAxis() +scale_colour_gradient2(low = "darkblue", mid = "white", high = "darkred")
dev.off()

#Subset celltypes
Cluster_EN <- subset(cGAS_integrated, idents = "excitatory neurons")
Cluster_IN <- subset(cGAS_integrated, idents = "inhibitory neurons")
Cluster_MG <- subset(cGAS_integrated, idents = "microglia")
Cluster_AST <- subset(cGAS_integrated, idents = "astrocytes")
Cluster_OL <- subset(cGAS_integrated, idents = "oligodendrocytes")
Cluster_OPC <- subset(cGAS_integrated, idents = "OPCs")

saveRDS(Cluster_EN, file = "EN_cGAS_integrated_Annotated.rds")
saveRDS(Cluster_IN, file = "IN_cGAS_integrated_Annotated.rds")
saveRDS(Cluster_MG, file = "MG_cGAS_integrated_Annotated.rds")
saveRDS(Cluster_AST, file = "Astro_cGAS_integrated_Annotated.rds")
saveRDS(Cluster_OL, file = "Oligo_cGAS_integrated_Annotated.rds")
saveRDS(Cluster_OPC, file = "OPC_cGAS_integrated_Annotated.rds")


# We further subclustered cluster_4 and found many cells expressing Gad1 and Gad2.
Idents(cGAS_integrated) <- 'seurat_clusters'
pdf("FeaturePlot_cGAS_integrated_Gad1_Slc17a7.pdf", width=16, height=8)
FeaturePlot(cGAS_integrated, features = c("Gad1","Slc17a7"), label = T, ncol = 2)
dev.off()

cluster_4 <- subset(cGAS_integrated, idents = "4")

DefaultAssay(cluster_4) <- 'integrated'
cluster_4 <- ScaleData(cluster_4, verbose = FALSE)
cluster_4 <- RunPCA(cluster_4, features = VariableFeatures(object = cluster_4), verbose = FALSE)
ElbowPlot(cluster_4)
cluster_4 <- FindNeighbors(cluster_4, dims = 1:15)
cluster_4 <- FindClusters(cluster_4, resolution = 0.1)
cluster_4 <- RunUMAP(cluster_4, dims = 1: 15)

DefaultAssay(cluster_4) <- "RNA"
markers.to.plot <- c("Plp1", "Mbp", "Mobp","Clu", "Aldoc","Pla2g7","Slc17a7","Nrgn",
                     "Cx3cr1", "P2ry12", "Csf1r", "Gad1", "Gad2",  
                     "Scrg1", "Pdgfra", "Ttr", "Vtn", "Igfbp7",
                     "Bnc2", "Slc47a1")
pdf("DotPlot_cluster_4_reclustered.pdf", width=12, height=8)
DotPlot(object = cluster_4, features = rev(x = markers.to.plot)) + RotatedAxis()
dev.off()

IN_from_4 <- subset(cluster_4, idents=c("2","4","5","9"))

IN <- readRDS(file = "IN_cGAS_integrated_Annotated.rds")

IN <- merge(IN, y = IN_from_4)

DefaultAssay(IN) <- 'integrated'
IN <- NormalizeData(IN, normalization.method = "LogNormalize", scale.factor = 10000)
IN <- FindVariableFeatures(IN, selection.method = "vst")
#scaling the data
IN <- ScaleData(object = IN)
#perform and visualize PCA
IN <- RunPCA(object = IN, features = rownames(x = IN), verbose = FALSE)

DefaultAssay(IN) <- 'RNA'
IN <- NormalizeData(IN, normalization.method = "LogNormalize", scale.factor = 10000)
IN <- ScaleData(IN, features = rownames(IN))

saveRDS(IN, file = "IN_cGAS_integrated_Annotated_new.rds")



