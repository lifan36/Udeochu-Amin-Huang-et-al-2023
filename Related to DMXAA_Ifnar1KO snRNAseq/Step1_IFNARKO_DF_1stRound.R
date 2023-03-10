#set working directory ====
setwd("/athena/ganlab/scratch/lif4001/IFNARKO/DF_1stRound")
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(SoupX)

#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/IFNARKO/cellranger/IFNARKO_1/outs')
sc = autoEstCont(sc)
IFNARKO_1.counts = adjustCounts(sc)
IFNARKO_1 <- CreateSeuratObject(counts = IFNARKO_1.counts, project = "IFNARKO_1", min.cells = 3, min.features = 200)
IFNARKO_1[["Condition"]] = c('IFNARKO_DMXAA')
IFNARKO_1[["Sample_Name"]] = c('IFNARKO_DMXAA_1')
rm(IFNARKO_1.counts)
#vizualize QC metrics and filtering====
IFNARKO_1[["percent.mt"]] <- PercentageFeatureSet(object = IFNARKO_1, pattern = "^mt-") #recognize mitochondrial transcripts
all <- IFNARKO_1
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("IFNARKO_1_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("IFNARKO_1_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"IFNARKO_1_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("IFNARKO_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/IFNARKO/cellranger/IFNARKO_2/outs')
sc = autoEstCont(sc)
IFNARKO_2.counts = adjustCounts(sc)
IFNARKO_2 <- CreateSeuratObject(counts = IFNARKO_2.counts, project = "IFNARKO_2", min.cells = 3, min.features = 200)
IFNARKO_2[["Condition"]] = c('IFNARKO_DMXAA')
IFNARKO_2[["Sample_Name"]] = c('IFNARKO_DMXAA_2')
rm(IFNARKO_2.counts)
#vizualize QC metrics and filtering====
IFNARKO_2[["percent.mt"]] <- PercentageFeatureSet(object = IFNARKO_2, pattern = "^mt-") #recognize mitochondrial transcripts
all <- IFNARKO_2
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("IFNARKO_2_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("IFNARKO_2_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"IFNARKO_2_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("IFNARKO_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/IFNARKO/cellranger/IFNARKO_3/outs')
sc = autoEstCont(sc)
IFNARKO_3.counts = adjustCounts(sc)
IFNARKO_3 <- CreateSeuratObject(counts = IFNARKO_3.counts, project = "IFNARKO_3", min.cells = 3, min.features = 200)
IFNARKO_3[["Condition"]] = c('IFNARKO_DMXAA')
IFNARKO_3[["Sample_Name"]] = c('IFNARKO_DMXAA_3')
rm(IFNARKO_3.counts)
#vizualize QC metrics and filtering====
IFNARKO_3[["percent.mt"]] <- PercentageFeatureSet(object = IFNARKO_3, pattern = "^mt-") #recognize mitochondrial transcripts
all <- IFNARKO_3
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("IFNARKO_3_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("IFNARKO_3_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"IFNARKO_3_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("IFNARKO_3_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/IFNARKO/cellranger/IFNARKO_5/outs')
sc = autoEstCont(sc)
IFNARKO_4.counts = adjustCounts(sc)
IFNARKO_4 <- CreateSeuratObject(counts = IFNARKO_4.counts, project = "IFNARKO_5", min.cells = 3, min.features = 200)
IFNARKO_4[["Condition"]] = c('IFNARKO_DMXAA')
IFNARKO_4[["Sample_Name"]] = c('IFNARKO_DMXAA_4')
rm(IFNARKO_4.counts)
#vizualize QC metrics and filtering====
IFNARKO_4[["percent.mt"]] <- PercentageFeatureSet(object = IFNARKO_4, pattern = "^mt-") #recognize mitochondrial transcripts
all <- IFNARKO_4
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("IFNARKO_4_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("IFNARKO_4_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"IFNARKO_4_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("IFNARKO_4_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/IFNARKO/cellranger/IFNARKO_7/outs')
sc = autoEstCont(sc)
Vehicle_1.counts = adjustCounts(sc)
Vehicle_1 <- CreateSeuratObject(counts = Vehicle_1.counts, project = "IFNARKO_7", min.cells = 3, min.features = 200)
Vehicle_1[["Condition"]] = c('IFNARKO_Vehicle')
Vehicle_1[["Sample_Name"]] = c('IFNARKO_Vehicle_1')
rm(Vehicle_1.counts)
#vizualize QC metrics and filtering====
Vehicle_1[["percent.mt"]] <- PercentageFeatureSet(object = Vehicle_1, pattern = "^mt-") #recognize mitochondrial transcripts
all <- Vehicle_1
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Vehicle_1_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Vehicle_1_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"Vehicle_1_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Vehicle_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/IFNARKO/cellranger/IFNARKO_8/outs')
sc = autoEstCont(sc)
Vehicle_2.counts = adjustCounts(sc)
Vehicle_2 <- CreateSeuratObject(counts = Vehicle_2.counts, project = "IFNARKO_8", min.cells = 3, min.features = 200)
Vehicle_2[["Condition"]] = c('IFNARKO_Vehicle')
Vehicle_2[["Sample_Name"]] = c('IFNARKO_Vehicle_2')
rm(Vehicle_2.counts)
#vizualize QC metrics and filtering====
Vehicle_2[["percent.mt"]] <- PercentageFeatureSet(object = Vehicle_2, pattern = "^mt-") #recognize mitochondrial transcripts
all <- Vehicle_2
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Vehicle_2_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Vehicle_2_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"Vehicle_2_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Vehicle_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/IFNARKO/cellranger/IFNARKO_9/outs')
sc = autoEstCont(sc)
Vehicle_3.counts = adjustCounts(sc)
Vehicle_3 <- CreateSeuratObject(counts = Vehicle_3.counts, project = "IFNARKO_9", min.cells = 3, min.features = 200)
Vehicle_3[["Condition"]] = c('IFNARKO_Vehicle')
Vehicle_3[["Sample_Name"]] = c('IFNARKO_Vehicle_3')
rm(Vehicle_3.counts)
#vizualize QC metrics and filtering====
Vehicle_3[["percent.mt"]] <- PercentageFeatureSet(object = Vehicle_3, pattern = "^mt-") #recognize mitochondrial transcripts
all <- Vehicle_3
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Vehicle_3_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Vehicle_3_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"Vehicle_3_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Vehicle_3_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
sc = load10X('/athena/ganlab/scratch/lif4001/IFNARKO/cellranger/IFNARKO_10/outs')
sc = autoEstCont(sc)
Vehicle_4.counts = adjustCounts(sc)
Vehicle_4 <- CreateSeuratObject(counts = Vehicle_4.counts, project = "IFNARKO_10", min.cells = 3, min.features = 200)
Vehicle_4[["Condition"]] = c('IFNARKO_Vehicle')
Vehicle_4[["Sample_Name"]] = c('IFNARKO_Vehicle_4')
rm(Vehicle_4.counts)
#vizualize QC metrics and filtering====
Vehicle_4[["percent.mt"]] <- PercentageFeatureSet(object = Vehicle_4, pattern = "^mt-") #recognize mitochondrial transcripts
all <- Vehicle_4
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("Vehicle_4_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 1%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 15000 & percent.mt < 1)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("Vehicle_4_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
saveRDS(all,"Vehicle_4_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("Vehicle_4_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################

