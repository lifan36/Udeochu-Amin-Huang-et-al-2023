#set working directory ====
setwd("/athena/ganlab/scratch/jcu4001/JCUPs19cGas/TestRunJU1/TDI_sNucSeq150")
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(SoupX)

#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Joe_cGas_4th/Joe_C90/outs')
sc = autoEstCont(sc)
TDIC90.counts = adjustCounts(sc)
TDIC90 <- CreateSeuratObject(counts = TDIC90.counts, project = "C90", min.cells = 3, min.features = 200)
TDIC90[["Genotype"]] = c('P301S tau: -')
TDIC90[["Sample_Name"]] = c('Ntg_Ctrl_1')
TDIC90[["Condition_1"]] = c('Ntg Ctrl')
TDIC90[["Condition_2"]] = c('Ctrl diet')
TDIC90[["Batch"]] = c('Batch_2')
rm(TDIC90.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
TDIC90[["percent.mt"]] <- PercentageFeatureSet(object = TDIC90, pattern = "^mt-") #recognize mitochondrial transcripts
all <- TDIC90
pdf("TDIC90_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("TDIC90_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("TDIC90_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("TDIC90_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "TDIC90_cell_counts_Condition.csv")
saveRDS(all, 'C90_Step1.rds')
###############################################################################################

sc = load10X('/athena/ganlab/scratch/lif4001/Joe_cGas_4th/Joe_D45/outs')
sc = autoEstCont(sc)
TDID45.counts = adjustCounts(sc)
TDID45 <- CreateSeuratObject(counts = TDID45.counts, project = "D45", min.cells = 3, min.features = 200)
TDID45[["Genotype"]] = c('P301S tau: -')
TDID45[["Sample_Name"]] = c('Ntg_Ctrl_2')
TDID45[["Condition_1"]] = c('Ntg Ctrl')
TDID45[["Condition_2"]] = c('Ctrl diet')
TDID45[["Batch"]] = c('Batch_2')
rm(TDID45.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
TDID45[["percent.mt"]] <- PercentageFeatureSet(object = TDID45, pattern = "^mt-") #recognize mitochondrial transcripts
all <- TDID45
pdf("TDID45_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("TDID45_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("TDID45_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("TDID45_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "TDID45_cell_counts_Condition.csv")
saveRDS(all, 'D45_Step1.rds')
######################
sc = load10X('/athena/ganlab/scratch/lif4001/Joe_cGas_4th/Joe_D46/outs')
sc = autoEstCont(sc)
TDID46.counts = adjustCounts(sc)
TDID46 <- CreateSeuratObject(counts = TDID46.counts, project = "D46", min.cells = 3, min.features = 200)
TDID46[["Genotype"]] = c('P301S tau: -')
TDID46[["Sample_Name"]] = c('Ntg_Ctrl_3')
TDID46[["Condition_1"]] = c('Ntg Ctrl')
TDID46[["Condition_2"]] = c('Ctrl diet')
TDID46[["Batch"]] = c('Batch_2')
rm(TDID46.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
TDID46[["percent.mt"]] <- PercentageFeatureSet(object = TDID46, pattern = "^mt-") #recognize mitochondrial transcripts
all <- TDID46
pdf("TDID46_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("TDID46_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("TDID46_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("TDID46_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "TDID46_cell_counts_Condition.csv")
saveRDS(all, 'D46_Step1.rds')
#############################

sc = load10X('/athena/ganlab/scratch/lif4001/Joe_cGas_4th/Joe_D27/outs')
sc = autoEstCont(sc)
TDID27.counts = adjustCounts(sc)
TDID27 <- CreateSeuratObject(counts = TDID27.counts, project = "D27", min.cells = 3, min.features = 200)
TDID27[["Genotype"]] = c('P301S tau: -')
TDID27[["Sample_Name"]] = c('Ntg_Ctrl_4')
TDID27[["Condition_1"]] = c('Ntg Ctrl')
TDID27[["Condition_2"]] = c('Ctrl diet')
TDID27[["Batch"]] = c('Batch_2')
rm(TDID27.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
TDID27[["percent.mt"]] <- PercentageFeatureSet(object = TDID27, pattern = "^mt-") #recognize mitochondrial transcripts
all <- TDID27
pdf("TDID27_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("TDID27_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("TDID27_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("TDID27_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "TDID27_cell_counts_Condition.csv")
saveRDS(all, 'D27_Step1.rds')
################

sc = load10X('/athena/ganlab/scratch/lif4001/Joe_cGas_4th/Joe_C95/outs')
sc = autoEstCont(sc)
TDIC95.counts = adjustCounts(sc)
TDIC95 <- CreateSeuratObject(counts = TDIC95.counts, project = "C95", min.cells = 3, min.features = 200)
TDIC95[["Genotype"]] = c('P301S tau: -')
TDIC95[["Sample_Name"]] = c('Ntg_TDI_1')
TDIC95[["Condition_1"]] = c('Ntg TDI')
TDIC95[["Condition_2"]] = c('TDI diet')
TDIC95[["Batch"]] = c('Batch_2')
rm(TDIC95.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
TDIC95[["percent.mt"]] <- PercentageFeatureSet(object = TDIC95, pattern = "^mt-") #recognize mitochondrial transcripts
all <- TDIC95
pdf("TDIC95_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("TDIC95_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("TDIC95_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("TDIC95_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "TDIC95_cell_counts_Condition.csv")
saveRDS(all, 'C95_Step1.rds')
##################

sc = load10X('/athena/ganlab/scratch/lif4001/Joe_cGas_4th/Joe_D18/outs')
sc = autoEstCont(sc)
TDID18.counts = adjustCounts(sc)
TDID18 <- CreateSeuratObject(counts = TDID18.counts, project = "D18", min.cells = 3, min.features = 200)
TDID18[["Genotype"]] = c('P301S tau: -')
TDID18[["Sample_Name"]] = c('Ntg_TDI_2')
TDID18[["Condition_1"]] = c('Ntg TDI')
TDID18[["Condition_2"]] = c('TDI diet')
TDID18[["Batch"]] = c('Batch_2')
rm(TDID18.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
TDID18[["percent.mt"]] <- PercentageFeatureSet(object = TDID18, pattern = "^mt-") #recognize mitochondrial transcripts
all <- TDID18
pdf("TDID18_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("TDID18_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("TDID18_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("TDID18_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "TDID18_cell_counts_Condition.csv")
saveRDS(all, 'D18_Step1.rds')
##########################
sc = load10X('/athena/ganlab/scratch/lif4001/Joe_cGas_4th/Joe_D17/outs')
sc = autoEstCont(sc)
TDID17.counts = adjustCounts(sc)
TDID17 <- CreateSeuratObject(counts = TDID17.counts, project = "D17", min.cells = 3, min.features = 200)
TDID17[["Genotype"]] = c('P301S tau: -')
TDID17[["Sample_Name"]] = c('Ntg_TDI_3')
TDID17[["Condition_1"]] = c('Ntg TDI')
TDID17[["Condition_2"]] = c('TDI diet')
TDID17[["Batch"]] = c('Batch_2')
rm(TDID17.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
TDID17[["percent.mt"]] <- PercentageFeatureSet(object = TDID17, pattern = "^mt-") #recognize mitochondrial transcripts
all <- TDID17
pdf("TDID17_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("TDID17_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("TDID17_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("TDID17_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "TDID17_cell_counts_Condition.csv")
saveRDS(all, 'D17_Step1.rds')
#################
sc = load10X('/athena/ganlab/scratch/lif4001/Joe_cGas_4th/Joe_D15/outs')
sc = autoEstCont(sc)
TDID15.counts = adjustCounts(sc)
TDID15 <- CreateSeuratObject(counts = TDID15.counts, project = "D15", min.cells = 3, min.features = 200)
TDID15[["Genotype"]] = c('P301S tau: -')
TDID15[["Sample_Name"]] = c('Ntg_TDI_4')
TDID15[["Condition_1"]] = c('Ntg TDI')
TDID15[["Condition_2"]] = c('TDI diet')
TDID15[["Batch"]] = c('Batch_2')
rm(TDID15.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
TDID15[["percent.mt"]] <- PercentageFeatureSet(object = TDID15, pattern = "^mt-") #recognize mitochondrial transcripts
all <- TDID15
pdf("TDID15_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("TDID15_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("TDID15_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("TDID15_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "TDID15_cell_counts_Condition.csv")
saveRDS(all, 'D15_Step1.rds')
###########################
sc = load10X('/athena/ganlab/scratch/lif4001/Joe_cGas_4th/Joe_C96/outs')
sc = autoEstCont(sc)
TDIC96.counts = adjustCounts(sc)
TDIC96 <- CreateSeuratObject(counts = TDIC96.counts, project = "C96", min.cells = 3, min.features = 200)
TDIC96[["Genotype"]] = c('P301S tau: +')
TDIC96[["Sample_Name"]] = c('P301S_TDI_1')
TDIC96[["Condition_1"]] = c('P301S TDI')
TDIC96[["Condition_2"]] = c('TDI diet')
TDIC96[["Batch"]] = c('Batch_1')
rm(TDIC96.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
TDIC96[["percent.mt"]] <- PercentageFeatureSet(object = TDIC96, pattern = "^mt-") #recognize mitochondrial transcripts
all <- TDIC96
pdf("TDIC96_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("TDIC96_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("TDIC96_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("TDIC96_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "TDIC96_cell_counts_Condition.csv")
saveRDS(all, 'C96_Step1.rds')
################
sc = load10X('/athena/ganlab/scratch/lif4001/Joe_cGas_4th/Joe_D13/outs')
sc = autoEstCont(sc)
TDID13.counts = adjustCounts(sc)
TDID13 <- CreateSeuratObject(counts = TDID13.counts, project = "D13", min.cells = 3, min.features = 200)
TDID13[["Genotype"]] = c('P301S tau: +')
TDID13[["Sample_Name"]] = c('P301S_TDI_2')
TDID13[["Condition_1"]] = c('P301S TDI')
TDID13[["Condition_2"]] = c('TDI diet')
TDID13[["Batch"]] = c('Batch_1')
rm(TDID13.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
TDID13[["percent.mt"]] <- PercentageFeatureSet(object = TDID13, pattern = "^mt-") #recognize mitochondrial transcripts
all <- TDID13
pdf("TDID13_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("TDID13_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("TDID13_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("TDID13_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "TDID13_cell_counts_Condition.csv")
saveRDS(all, 'D13_Step1.rds')
##############
sc = load10X('/athena/ganlab/scratch/lif4001/Joe_cGas_4th/Joe_D22/outs')
sc = autoEstCont(sc)
TDID22.counts = adjustCounts(sc)
TDID22 <- CreateSeuratObject(counts = TDID22.counts, project = "D22", min.cells = 3, min.features = 200)
TDID22[["Genotype"]] = c('P301S tau: +')
TDID22[["Sample_Name"]] = c('P301S_TDI_3')
TDID22[["Condition_1"]] = c('P301S TDI')
TDID22[["Condition_2"]] = c('TDI diet')
TDID22[["Batch"]] = c('Batch_1')
rm(TDID22.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
TDID22[["percent.mt"]] <- PercentageFeatureSet(object = TDID22, pattern = "^mt-") #recognize mitochondrial transcripts
all <- TDID22
pdf("TDID22_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("TDID22_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("TDID22_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("TDID22_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "TDID22_cell_counts_Condition.csv")
saveRDS(all, 'D22_Step1.rds')
####################
sc = load10X('/athena/ganlab/scratch/lif4001/Joe_cGas_4th/Joe_D50/outs')
sc = autoEstCont(sc)
TDID50.counts = adjustCounts(sc)
TDID50 <- CreateSeuratObject(counts = TDID50.counts, project = "D50", min.cells = 3, min.features = 200)
TDID50[["Genotype"]] = c('P301S tau: +')
TDID50[["Sample_Name"]] = c('P301S_TDI_4')
TDID50[["Condition_1"]] = c('P301S TDI')
TDID50[["Condition_2"]] = c('TDI diet')
TDID50[["Batch"]] = c('Batch_1')
rm(TDID50.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
TDID50[["percent.mt"]] <- PercentageFeatureSet(object = TDID50, pattern = "^mt-") #recognize mitochondrial transcripts
all <- TDID50
pdf("TDID50_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("TDID50_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("TDID50_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("TDID50_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "TDID50_cell_counts_Condition.csv")
saveRDS(all, 'D50_Step1.rds')
##########
sc = load10X('/athena/ganlab/scratch/lif4001/Joe_cGas_4th/Joe_C92/outs')
sc = autoEstCont(sc)
TDIC92.counts = adjustCounts(sc)
TDIC92 <- CreateSeuratObject(counts = TDIC92.counts, project = "C92", min.cells = 3, min.features = 200)
TDIC92[["Genotype"]] = c('P301S tau: +')
TDIC92[["Sample_Name"]] = c('P301S_Ctrl_1')
TDIC92[["Condition_1"]] = c('P301S Ctrl')
TDIC92[["Condition_2"]] = c('Ctrl diet')
TDIC92[["Batch"]] = c('Batch_1')
rm(TDIC92.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
TDIC92[["percent.mt"]] <- PercentageFeatureSet(object = TDIC92, pattern = "^mt-") #recognize mitochondrial transcripts
all <- TDIC92
pdf("TDIC92_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("TDIC92_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("TDIC92_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("TDIC92_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "TDIC92_cell_counts_Condition.csv")
saveRDS(all, 'C92_Step1.rds')
###################
sc = load10X('/athena/ganlab/scratch/lif4001/Joe_cGas_4th/Joe_C94/outs')
sc = autoEstCont(sc)
TDIC94.counts = adjustCounts(sc)
TDIC94 <- CreateSeuratObject(counts = TDIC94.counts, project = "C94", min.cells = 3, min.features = 200)
TDIC94[["Genotype"]] = c('P301S tau: +')
TDIC94[["Sample_Name"]] = c('P301S_Ctrl_2')
TDIC94[["Condition_1"]] = c('P301S Ctrl')
TDIC94[["Condition_2"]] = c('Ctrl diet')
TDIC94[["Batch"]] = c('Batch_1')
rm(TDIC94.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
TDIC94[["percent.mt"]] <- PercentageFeatureSet(object = TDIC94, pattern = "^mt-") #recognize mitochondrial transcripts
all <- TDIC94
pdf("TDIC94_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("TDIC94_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("TDIC94_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("TDIC94_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "TDIC94_cell_counts_Condition.csv")
saveRDS(all, 'C94_Step1.rds')
#####################

sc = load10X('/athena/ganlab/scratch/lif4001/Joe_cGas_4th/Joe_D49/outs')
sc = autoEstCont(sc)
TDID49.counts = adjustCounts(sc)
TDID49 <- CreateSeuratObject(counts = TDID49.counts, project = "D49", min.cells = 3, min.features = 200)
TDID49[["Genotype"]] = c('P301S tau: +')
TDID49[["Sample_Name"]] = c('P301S_Ctrl_3')
TDID49[["Condition_1"]] = c('P301S Ctrl')
TDID49[["Condition_2"]] = c('Ctrl diet')
TDID49[["Batch"]] = c('Batch_1')
rm(TDID49.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
TDID49[["percent.mt"]] <- PercentageFeatureSet(object = TDID49, pattern = "^mt-") #recognize mitochondrial transcripts
all <- TDID49
pdf("TDID49_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("TDID49_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("TDID49_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("TDID49_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "TDID49_cell_counts_Condition.csv")
saveRDS(all, 'D49_Step1.rds')
##############

sc = load10X('/athena/ganlab/scratch/lif4001/Joe_cGas_4th/Joe_D28/outs')
sc = autoEstCont(sc)
TDID28.counts = adjustCounts(sc)
TDID28 <- CreateSeuratObject(counts = TDID28.counts, project = "D28", min.cells = 3, min.features = 200)
TDID28[["Genotype"]] = c('P301S tau: +')
TDID28[["Sample_Name"]] = c('P301S_Ctrl_4')
TDID28[["Condition_1"]] = c('P301S Ctrl')
TDID28[["Condition_2"]] = c('Ctrl diet')
TDID28[["Batch"]] = c('Batch_1')
rm(TDID28.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
TDID28[["percent.mt"]] <- PercentageFeatureSet(object = TDID28, pattern = "^mt-") #recognize mitochondrial transcripts
all <- TDID28
pdf("TDID28_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("TDID28_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("TDID28_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("TDID28_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "TDID28_cell_counts_Condition.csv")
saveRDS(all, 'D28_Step1.rds')