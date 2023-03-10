#set working directory ====
setwd("/athena/ganlab/scratch/lif4001/Joe_cGas/data_analysis/DF_2ndRound")
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(SoupX)


#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/NovaSeqS2/cellranger/Joe-7_25/outs')
sc = autoEstCont(sc)
cGAS_25.counts = adjustCounts(sc)
cGAS_25 <- CreateSeuratObject(counts = cGAS_25.counts, project = "25", min.cells = 3, min.features = 200)
cGAS_25[["Genotype"]] = c('P301S tau: -; cGAS: -/-')
cGAS_25[["Sample_Name"]] = c('NT_cGAS_Homo_1')
cGAS_25[["Condition_1"]] = c('NT_cGAS_Homo_M')
cGAS_25[["Condition_2"]] = c('NT_cGAS_Homo')
cGAS_25[["Condition_3"]] = c('NT_M')
cGAS_25[["Condition_4"]] = c('cGAS_Homo_M')
cGAS_25[["PS19"]] = c('NT')
cGAS_25[["cGAS"]] = c('cGAS_Homo')
cGAS_25[["Sex"]] = c('M')
cGAS_25[["Batch"]] = c('Batch_2')
rm(cGAS_25.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
cGAS_25[["percent.mt"]] <- PercentageFeatureSet(object = cGAS_25, pattern = "^mt-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/NovaSeqS2/cellranger/Joe-1_23/outs')
sc = autoEstCont(sc)
cGAS_23.counts = adjustCounts(sc)
cGAS_23 <- CreateSeuratObject(counts = cGAS_23.counts, project = "23", min.cells = 3, min.features = 200)
cGAS_23[["Genotype"]] = c('P301S tau: -; cGAS: -/-')
cGAS_23[["Sample_Name"]] = c('NT_cGAS_Homo_2')
cGAS_23[["Condition_1"]] = c('NT_cGAS_Homo_M')
cGAS_23[["Condition_2"]] = c('NT_cGAS_Homo')
cGAS_23[["Condition_3"]] = c('NT_M')
cGAS_23[["Condition_4"]] = c('cGAS_Homo_M')
cGAS_23[["PS19"]] = c('NT')
cGAS_23[["cGAS"]] = c('cGAS_Homo')
cGAS_23[["Sex"]] = c('M')
cGAS_23[["Batch"]] = c('Batch_1')
rm(cGAS_23.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
cGAS_23[["percent.mt"]] <- PercentageFeatureSet(object = cGAS_23, pattern = "^mt-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/NovaSeqS2/cellranger/Joe-2_18/outs')
sc = autoEstCont(sc)
cGAS_18.counts = adjustCounts(sc)
cGAS_18 <- CreateSeuratObject(counts = cGAS_18.counts, project = "18", min.cells = 3, min.features = 200)
cGAS_18[["Genotype"]] = c('P301S tau: -; cGAS: +/-')
cGAS_18[["Sample_Name"]] = c('NT_cGAS_Hetero_1')
cGAS_18[["Condition_1"]] = c('NT_cGAS_Hetero_M')
cGAS_18[["Condition_2"]] = c('NT_cGAS_Hetero')
cGAS_18[["Condition_3"]] = c('NT_M')
cGAS_18[["Condition_4"]] = c('cGAS_Hetero_M')
cGAS_18[["PS19"]] = c('NT')
cGAS_18[["cGAS"]] = c('cGAS_Hetero')
cGAS_18[["Sex"]] = c('M')
cGAS_18[["Batch"]] = c('Batch_1')
rm(cGAS_18.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
cGAS_18[["percent.mt"]] <- PercentageFeatureSet(object = cGAS_18, pattern = "^mt-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/NovaSeqS2/cellranger/Joe-8_20/outs')
sc = autoEstCont(sc)
cGAS_20.counts = adjustCounts(sc)
cGAS_20 <- CreateSeuratObject(counts = cGAS_20.counts, project = "20", min.cells = 3, min.features = 200)
cGAS_20[["Genotype"]] = c('P301S tau: -; cGAS: +/-')
cGAS_20[["Sample_Name"]] = c('NT_cGAS_Hetero_2')
cGAS_20[["Condition_1"]] = c('NT_cGAS_Hetero_M')
cGAS_20[["Condition_2"]] = c('NT_cGAS_Hetero')
cGAS_20[["Condition_3"]] = c('NT_M')
cGAS_20[["Condition_4"]] = c('cGAS_Hetero_M')
cGAS_20[["PS19"]] = c('NT')
cGAS_20[["cGAS"]] = c('cGAS_Hetero')
cGAS_20[["Sex"]] = c('M')
cGAS_20[["Batch"]] = c('Batch_2')
rm(cGAS_20.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
cGAS_20[["percent.mt"]] <- PercentageFeatureSet(object = cGAS_20, pattern = "^mt-") #recognize mitochondrial transcripts


sc = load10X('/athena/ganlab/scratch/lif4001/NovaSeqS2/cellranger/Joe-3_55/outs')
sc = autoEstCont(sc)
cGAS_55.counts = adjustCounts(sc)
cGAS_55 <- CreateSeuratObject(counts = cGAS_55.counts, project = "55", min.cells = 3, min.features = 200)
cGAS_55[["Genotype"]] = c('P301S tau: -; cGAS: +/+')
cGAS_55[["Sample_Name"]] = c('NT_cGAS_WT_1')
cGAS_55[["Condition_1"]] = c('NT_cGAS_WT_M')
cGAS_55[["Condition_2"]] = c('NT_cGAS_WT')
cGAS_55[["Condition_3"]] = c('NT_M')
cGAS_55[["Condition_4"]] = c('cGAS_WT_M')
cGAS_55[["PS19"]] = c('NT')
cGAS_55[["cGAS"]] = c('cGAS_WT')
cGAS_55[["Sex"]] = c('M')
cGAS_55[["Batch"]] = c('Batch_1')
rm(cGAS_55.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
cGAS_55[["percent.mt"]] <- PercentageFeatureSet(object = cGAS_55, pattern = "^mt-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/NovaSeqS2/cellranger/Joe-9_61/outs')
sc = autoEstCont(sc)
cGAS_61.counts = adjustCounts(sc)
cGAS_61 <- CreateSeuratObject(counts = cGAS_61.counts, project = "61", min.cells = 3, min.features = 200)
cGAS_61[["Genotype"]] = c('P301S tau: -; cGAS: +/+')
cGAS_61[["Sample_Name"]] = c('NT_cGAS_WT_2')
cGAS_61[["Condition_1"]] = c('NT_cGAS_WT_M')
cGAS_61[["Condition_2"]] = c('NT_cGAS_WT')
cGAS_61[["Condition_3"]] = c('NT_M')
cGAS_61[["Condition_4"]] = c('cGAS_WT_M')
cGAS_61[["PS19"]] = c('NT')
cGAS_61[["cGAS"]] = c('cGAS_WT')
cGAS_61[["Sex"]] = c('M')
cGAS_61[["Batch"]] = c('Batch_2')
rm(cGAS_61.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
cGAS_61[["percent.mt"]] <- PercentageFeatureSet(object = cGAS_61, pattern = "^mt-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/NovaSeqS2/cellranger/Joe-4_21/outs')
sc = autoEstCont(sc)
cGAS_21.counts = adjustCounts(sc)
cGAS_21 <- CreateSeuratObject(counts = cGAS_21.counts, project = "21", min.cells = 3, min.features = 200)
cGAS_21[["Genotype"]] = c('P301S tau: +; cGAS: -/-')
cGAS_21[["Sample_Name"]] = c('PS19_cGAS_Homo_1')
cGAS_21[["Condition_1"]] = c('PS19_cGAS_Homo_M')
cGAS_21[["Condition_2"]] = c('PS19_cGAS_Homo')
cGAS_21[["Condition_3"]] = c('PS19_M')
cGAS_21[["Condition_4"]] = c('cGAS_Homo_M')
cGAS_21[["PS19"]] = c('PS19')
cGAS_21[["cGAS"]] = c('cGAS_Homo')
cGAS_21[["Sex"]] = c('M')
cGAS_21[["Batch"]] = c('Batch_1')
rm(cGAS_21.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
cGAS_21[["percent.mt"]] <- PercentageFeatureSet(object = cGAS_21, pattern = "^mt-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/NovaSeqS2/cellranger/Joe-10_22/outs')
sc = autoEstCont(sc)
cGAS_22.counts = adjustCounts(sc)
cGAS_22 <- CreateSeuratObject(counts = cGAS_22.counts, project = "22", min.cells = 3, min.features = 200)
cGAS_22[["Genotype"]] = c('P301S tau: +; cGAS: -/-')
cGAS_22[["Sample_Name"]] = c('PS19_cGAS_Homo_2')
cGAS_22[["Condition_1"]] = c('PS19_cGAS_Homo_M')
cGAS_22[["Condition_2"]] = c('PS19_cGAS_Homo')
cGAS_22[["Condition_3"]] = c('PS19_M')
cGAS_22[["Condition_4"]] = c('cGAS_Homo_M')
cGAS_22[["PS19"]] = c('PS19')
cGAS_22[["cGAS"]] = c('cGAS_Homo')
cGAS_22[["Sex"]] = c('M')
cGAS_22[["Batch"]] = c('Batch_2')
rm(cGAS_22.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
cGAS_22[["percent.mt"]] <- PercentageFeatureSet(object = cGAS_22, pattern = "^mt-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/NovaSeqS2/cellranger/Joe-5_19/outs')
sc = autoEstCont(sc)
cGAS_19.counts = adjustCounts(sc)
cGAS_19 <- CreateSeuratObject(counts = cGAS_19.counts, project = "19", min.cells = 3, min.features = 200)
cGAS_19[["Genotype"]] = c('P301S tau: +; cGAS: +/-')
cGAS_19[["Sample_Name"]] = c('PS19_cGAS_Hetero_1')
cGAS_19[["Condition_1"]] = c('PS19_cGAS_Hetero_M')
cGAS_19[["Condition_2"]] = c('PS19_cGAS_Hetero')
cGAS_19[["Condition_3"]] = c('PS19_M')
cGAS_19[["Condition_4"]] = c('cGAS_Hetero_M')
cGAS_19[["PS19"]] = c('PS19')
cGAS_19[["cGAS"]] = c('cGAS_Hetero')
cGAS_19[["Sex"]] = c('M')
cGAS_19[["Batch"]] = c('Batch_1')
rm(cGAS_19.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
cGAS_19[["percent.mt"]] <- PercentageFeatureSet(object = cGAS_19, pattern = "^mt-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/NovaSeqS2/cellranger/Joe-11_30/outs')
sc = autoEstCont(sc)
cGAS_30.counts = adjustCounts(sc)
cGAS_30 <- CreateSeuratObject(counts = cGAS_30.counts, project = "30", min.cells = 3, min.features = 200)
cGAS_30[["Genotype"]] = c('P301S tau: +; cGAS: +/-')
cGAS_30[["Sample_Name"]] = c('PS19_cGAS_Hetero_2')
cGAS_30[["Condition_1"]] = c('PS19_cGAS_Hetero_M')
cGAS_30[["Condition_2"]] = c('PS19_cGAS_Hetero')
cGAS_30[["Condition_3"]] = c('PS19_M')
cGAS_30[["Condition_4"]] = c('cGAS_Hetero_M')
cGAS_30[["PS19"]] = c('PS19')
cGAS_30[["cGAS"]] = c('cGAS_Hetero')
cGAS_30[["Sex"]] = c('M')
cGAS_30[["Batch"]] = c('Batch_2')
rm(cGAS_30.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
cGAS_30[["percent.mt"]] <- PercentageFeatureSet(object = cGAS_30, pattern = "^mt-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/NovaSeqS2/cellranger/Joe-6_29/outs')
sc = autoEstCont(sc)
cGAS_29.counts = adjustCounts(sc)
cGAS_29 <- CreateSeuratObject(counts = cGAS_29.counts, project = "29", min.cells = 3, min.features = 200)
cGAS_29[["Genotype"]] = c('P301S tau: +; cGAS: +/+')
cGAS_29[["Sample_Name"]] = c('PS19_cGAS_WT_1')
cGAS_29[["Condition_1"]] = c('PS19_cGAS_WT_M')
cGAS_29[["Condition_2"]] = c('PS19_cGAS_WT')
cGAS_29[["Condition_3"]] = c('PS19_M')
cGAS_29[["Condition_4"]] = c('cGAS_WT_M')
cGAS_29[["PS19"]] = c('PS19')
cGAS_29[["cGAS"]] = c('cGAS_WT')
cGAS_29[["Sex"]] = c('M')
cGAS_29[["Batch"]] = c('Batch_1')
rm(cGAS_29.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
cGAS_29[["percent.mt"]] <- PercentageFeatureSet(object = cGAS_29, pattern = "^mt-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/NovaSeqS2/cellranger/Joe-12_63/outs')
sc = autoEstCont(sc)
cGAS_63.counts = adjustCounts(sc)
cGAS_63 <- CreateSeuratObject(counts = cGAS_63.counts, project = "63", min.cells = 3, min.features = 200)
cGAS_63[["Genotype"]] = c('P301S tau: +; cGAS: +/+')
cGAS_63[["Sample_Name"]] = c('PS19_cGAS_WT_2')
cGAS_63[["Condition_1"]] = c('PS19_cGAS_WT_M')
cGAS_63[["Condition_2"]] = c('PS19_cGAS_WT')
cGAS_63[["Condition_3"]] = c('PS19_M')
cGAS_63[["Condition_4"]] = c('cGAS_WT_M')
cGAS_63[["PS19"]] = c('PS19')
cGAS_63[["cGAS"]] = c('cGAS_WT')
cGAS_63[["Sex"]] = c('M')
cGAS_63[["Batch"]] = c('Batch_2')
rm(cGAS_63.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
cGAS_63[["percent.mt"]] <- PercentageFeatureSet(object = cGAS_63, pattern = "^mt-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Joe_cGas_2nd/Joe_13_36/outs')
sc = autoEstCont(sc)
cGAS_36.counts = adjustCounts(sc)
cGAS_36 <- CreateSeuratObject(counts = cGAS_36.counts, project = "36", min.cells = 3, min.features = 200)
cGAS_36[["Genotype"]] = c('P301S tau: +; cGAS: +/+')
cGAS_36[["Sample_Name"]] = c('PS19_cGAS_WT_3')
cGAS_36[["Condition_1"]] = c('PS19_cGAS_WT_F')
cGAS_36[["Condition_2"]] = c('PS19_cGAS_WT')
cGAS_36[["Condition_3"]] = c('PS19_F')
cGAS_36[["Condition_4"]] = c('cGAS_WT_F')
cGAS_36[["PS19"]] = c('PS19')
cGAS_36[["cGAS"]] = c('cGAS_WT')
cGAS_36[["Sex"]] = c('F')
cGAS_36[["Batch"]] = c('Batch_3')
rm(cGAS_36.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
cGAS_36[["percent.mt"]] <- PercentageFeatureSet(object = cGAS_36, pattern = "^mt-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Joe_cGas_2nd/Joe_14_43/outs')
sc = autoEstCont(sc)
cGAS_43.counts = adjustCounts(sc)
cGAS_43 <- CreateSeuratObject(counts = cGAS_43.counts, project = "43", min.cells = 3, min.features = 200)
cGAS_43[["Genotype"]] = c('P301S tau: +; cGAS: +/+')
cGAS_43[["Sample_Name"]] = c('PS19_cGAS_WT_4')
cGAS_43[["Condition_1"]] = c('PS19_cGAS_WT_F')
cGAS_43[["Condition_2"]] = c('PS19_cGAS_WT')
cGAS_43[["Condition_3"]] = c('PS19_F')
cGAS_43[["Condition_4"]] = c('cGAS_WT_F')
cGAS_43[["PS19"]] = c('PS19')
cGAS_43[["cGAS"]] = c('cGAS_WT')
cGAS_43[["Sex"]] = c('F')
cGAS_43[["Batch"]] = c('Batch_3')
rm(cGAS_43.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
cGAS_43[["percent.mt"]] <- PercentageFeatureSet(object = cGAS_43, pattern = "^mt-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Joe_cGas_2nd/Joe_15_33/outs')
sc = autoEstCont(sc)
cGAS_33.counts = adjustCounts(sc)
cGAS_33 <- CreateSeuratObject(counts = cGAS_33.counts, project = "33", min.cells = 3, min.features = 200)
cGAS_33[["Genotype"]] = c('P301S tau: +; cGAS: +/-')
cGAS_33[["Sample_Name"]] = c('PS19_cGAS_Hetero_3')
cGAS_33[["Condition_1"]] = c('PS19_cGAS_Hetero_M')
cGAS_33[["Condition_2"]] = c('PS19_cGAS_Hetero')
cGAS_33[["Condition_3"]] = c('PS19_M')
cGAS_33[["Condition_4"]] = c('cGAS_Hetero_M')
cGAS_33[["PS19"]] = c('PS19')
cGAS_33[["cGAS"]] = c('cGAS_Hetero')
cGAS_33[["Sex"]] = c('M')
cGAS_33[["Batch"]] = c('Batch_3')
rm(cGAS_33.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
cGAS_33[["percent.mt"]] <- PercentageFeatureSet(object = cGAS_33, pattern = "^mt-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Joe_cGas_2nd/Joe_16_59/outs')
sc = autoEstCont(sc)
cGAS_59.counts = adjustCounts(sc)
cGAS_59 <- CreateSeuratObject(counts = cGAS_59.counts, project = "59", min.cells = 3, min.features = 200)
cGAS_59[["Genotype"]] = c('P301S tau: +; cGAS: +/-')
cGAS_59[["Sample_Name"]] = c('PS19_cGAS_Hetero_4')
cGAS_59[["Condition_1"]] = c('PS19_cGAS_Hetero_F')
cGAS_59[["Condition_2"]] = c('PS19_cGAS_Hetero')
cGAS_59[["Condition_3"]] = c('PS19_F')
cGAS_59[["Condition_4"]] = c('cGAS_Hetero_F')
cGAS_59[["PS19"]] = c('PS19')
cGAS_59[["cGAS"]] = c('cGAS_Hetero')
cGAS_59[["Sex"]] = c('F')
cGAS_59[["Batch"]] = c('Batch_3')
rm(cGAS_59.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
cGAS_59[["percent.mt"]] <- PercentageFeatureSet(object = cGAS_59, pattern = "^mt-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Joe_cGas_2nd/Joe_17_42/outs')
sc = autoEstCont(sc)
cGAS_42.counts = adjustCounts(sc)
cGAS_42 <- CreateSeuratObject(counts = cGAS_42.counts, project = "42", min.cells = 3, min.features = 200)
cGAS_42[["Genotype"]] = c('P301S tau: +; cGAS: -/-')
cGAS_42[["Sample_Name"]] = c('PS19_cGAS_Homo_3')
cGAS_42[["Condition_1"]] = c('PS19_cGAS_Homo_F')
cGAS_42[["Condition_2"]] = c('PS19_cGAS_Homo')
cGAS_42[["Condition_3"]] = c('PS19_F')
cGAS_42[["Condition_4"]] = c('cGAS_Homo_F')
cGAS_42[["PS19"]] = c('PS19')
cGAS_42[["cGAS"]] = c('cGAS_Homo')
cGAS_42[["Sex"]] = c('F')
cGAS_42[["Batch"]] = c('Batch_3')
rm(cGAS_42.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
cGAS_42[["percent.mt"]] <- PercentageFeatureSet(object = cGAS_42, pattern = "^mt-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Joe_cGas_2nd/Joe_18_74/outs')
sc = autoEstCont(sc)
cGAS_74.counts = adjustCounts(sc)
cGAS_74 <- CreateSeuratObject(counts = cGAS_74.counts, project = "74", min.cells = 3, min.features = 200)
cGAS_74[["Genotype"]] = c('P301S tau: +; cGAS: -/-')
cGAS_74[["Sample_Name"]] = c('PS19_cGAS_Homo_4')
cGAS_74[["Condition_1"]] = c('PS19_cGAS_Homo_F')
cGAS_74[["Condition_2"]] = c('PS19_cGAS_Homo')
cGAS_74[["Condition_3"]] = c('PS19_F')
cGAS_74[["Condition_4"]] = c('cGAS_Homo_F')
cGAS_74[["PS19"]] = c('PS19')
cGAS_74[["cGAS"]] = c('cGAS_Homo')
cGAS_74[["Sex"]] = c('F')
cGAS_74[["Batch"]] = c('Batch_3')
rm(cGAS_74.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
cGAS_74[["percent.mt"]] <- PercentageFeatureSet(object = cGAS_74, pattern = "^mt-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Joe_cGas_2nd/Joe_19_45/outs')
sc = autoEstCont(sc)
cGAS_45.counts = adjustCounts(sc)
cGAS_45 <- CreateSeuratObject(counts = cGAS_45.counts, project = "45", min.cells = 3, min.features = 200)
cGAS_45[["Genotype"]] = c('P301S tau: -; cGAS: +/+')
cGAS_45[["Sample_Name"]] = c('NT_cGAS_WT_3')
cGAS_45[["Condition_1"]] = c('NT_cGAS_WT_M')
cGAS_45[["Condition_2"]] = c('NT_cGAS_WT')
cGAS_45[["Condition_3"]] = c('NT_M')
cGAS_45[["Condition_4"]] = c('cGAS_WT_M')
cGAS_45[["PS19"]] = c('NT')
cGAS_45[["cGAS"]] = c('cGAS_WT')
cGAS_45[["Sex"]] = c('M')
cGAS_45[["Batch"]] = c('Batch_3')
rm(cGAS_45.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
cGAS_45[["percent.mt"]] <- PercentageFeatureSet(object = cGAS_45, pattern = "^mt-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/lif4001/Joe_cGas_2nd/Joe_20_46/outs')
sc = autoEstCont(sc)
cGAS_46.counts = adjustCounts(sc)
cGAS_46 <- CreateSeuratObject(counts = cGAS_46.counts, project = "46", min.cells = 3, min.features = 200)
cGAS_46[["Genotype"]] = c('P301S tau: -; cGAS: +/+')
cGAS_46[["Sample_Name"]] = c('NT_cGAS_WT_4')
cGAS_46[["Condition_1"]] = c('NT_cGAS_WT_M')
cGAS_46[["Condition_2"]] = c('NT_cGAS_WT')
cGAS_46[["Condition_3"]] = c('NT_M')
cGAS_46[["Condition_4"]] = c('cGAS_WT_M')
cGAS_46[["PS19"]] = c('NT')
cGAS_46[["cGAS"]] = c('cGAS_WT')
cGAS_46[["Sex"]] = c('M')
cGAS_46[["Batch"]] = c('Batch_3')
rm(cGAS_46.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
cGAS_46[["percent.mt"]] <- PercentageFeatureSet(object = cGAS_46, pattern = "^mt-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/jcu4001/JCUPs19cGas/TestRunJU1/Joe_cGas_3rd/Joe_21_27/outs')
sc = autoEstCont(sc)
cGAS_27.counts = adjustCounts(sc)
cGAS_27 <- CreateSeuratObject(counts = cGAS_27.counts, project = "27", min.cells = 3, min.features = 620)
cGAS_27[["Genotype"]] = c('P301S tau: +; cGAS: -/-')
cGAS_27[["Sample_Name"]] = c('PS19_cGAS_Homo_5')
cGAS_27[["Condition_1"]] = c('PS19_cGAS_Homo_M')
cGAS_27[["Condition_2"]] = c('PS19_cGAS_Homo')
cGAS_27[["Condition_3"]] = c('PS19_M')
cGAS_27[["Condition_4"]] = c('cGAS_Homo_M')
cGAS_27[["PS19"]] = c('PS19')
cGAS_27[["cGAS"]] = c('cGAS_Homo')
cGAS_27[["Sex"]] = c('M')
cGAS_27[["Batch"]] = c('Batch_4')
rm(cGAS_27.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
cGAS_27[["percent.mt"]] <- PercentageFeatureSet(object = cGAS_27, pattern = "^mt-") #recognize mitochondrial transcripts

###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/jcu4001/JCUPs19cGas/TestRunJU1/Joe_cGas_3rd/Joe_22_41/outs')
sc = autoEstCont(sc)
cGAS_41.counts = adjustCounts(sc)
cGAS_41 <- CreateSeuratObject(counts = cGAS_41.counts, project = "41", min.cells = 3, min.features = 620)
cGAS_41[["Genotype"]] = c('P301S tau: +; cGAS: -/-')
cGAS_41[["Sample_Name"]] = c('PS19_cGAS_Homo_6')
cGAS_41[["Condition_1"]] = c('PS19_cGAS_Homo_F')
cGAS_41[["Condition_2"]] = c('PS19_cGAS_Homo')
cGAS_41[["Condition_3"]] = c('PS19_F')
cGAS_41[["Condition_4"]] = c('cGAS_Homo_F')
cGAS_41[["PS19"]] = c('PS19')
cGAS_41[["cGAS"]] = c('cGAS_Homo')
cGAS_41[["Sex"]] = c('F')
cGAS_41[["Batch"]] = c('Batch_4')
rm(cGAS_41.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
cGAS_41[["percent.mt"]] <- PercentageFeatureSet(object = cGAS_41, pattern = "^mt-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/jcu4001/JCUPs19cGas/TestRunJU1/Joe_cGas_3rd/Joe_23_57/outs')
sc = autoEstCont(sc)
cGAS_57.counts = adjustCounts(sc)
cGAS_57 <- CreateSeuratObject(counts = cGAS_57.counts, project = "57", min.cells = 3, min.features = 620)
cGAS_57[["Genotype"]] = c('P301S tau: +; cGAS: +/-')
cGAS_57[["Sample_Name"]] = c('PS19_cGAS_Hetero_5')
cGAS_57[["Condition_1"]] = c('PS19_cGAS_Hetero_F')
cGAS_57[["Condition_2"]] = c('PS19_cGAS_Hetero')
cGAS_57[["Condition_3"]] = c('PS19_F')
cGAS_57[["Condition_4"]] = c('cGAS_Hetero_F')
cGAS_57[["PS19"]] = c('PS19')
cGAS_57[["cGAS"]] = c('cGAS_Hetero')
cGAS_57[["Sex"]] = c('F')
cGAS_57[["Batch"]] = c('Batch_4')
rm(cGAS_57.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
cGAS_57[["percent.mt"]] <- PercentageFeatureSet(object = cGAS_57, pattern = "^mt-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/jcu4001/JCUPs19cGas/TestRunJU1/Joe_cGas_3rd/Joe_24_62/outs')
sc = autoEstCont(sc)
cGAS_62.counts = adjustCounts(sc)
cGAS_62 <- CreateSeuratObject(counts = cGAS_62.counts, project = "62", min.cells = 3, min.features = 620)
cGAS_62[["Genotype"]] = c('P301S tau: -; cGAS: +/+')
cGAS_62[["Sample_Name"]] = c('NT_cGAS_WT_5')
cGAS_62[["Condition_1"]] = c('NT_cGAS_WT_M')
cGAS_62[["Condition_2"]] = c('NT_cGAS_WT')
cGAS_62[["Condition_3"]] = c('NT_M')
cGAS_62[["Condition_4"]] = c('cGAS_WT_M')
cGAS_62[["PS19"]] = c('NT')
cGAS_62[["cGAS"]] = c('cGAS_WT')
cGAS_62[["Sex"]] = c('M')
cGAS_62[["Batch"]] = c('Batch_4')
rm(cGAS_62.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
cGAS_62[["percent.mt"]] <- PercentageFeatureSet(object = cGAS_62, pattern = "^mt-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/jcu4001/JCUPs19cGas/TestRunJU1/Joe_cGas_3rd/Joe_25_65/outs')
sc = autoEstCont(sc)
cGAS_65.counts = adjustCounts(sc)
cGAS_65 <- CreateSeuratObject(counts = cGAS_65.counts, project = "65", min.cells = 3, min.features = 620)
cGAS_65[["Genotype"]] = c('P301S tau: +; cGAS: +/+')
cGAS_65[["Sample_Name"]] = c('PS19_cGAS_WT_5')
cGAS_65[["Condition_1"]] = c('PS19_cGAS_WT_F')
cGAS_65[["Condition_2"]] = c('PS19_cGAS_WT')
cGAS_65[["Condition_3"]] = c('PS19_F')
cGAS_65[["Condition_4"]] = c('cGAS_WT_F')
cGAS_65[["PS19"]] = c('PS19')
cGAS_65[["cGAS"]] = c('cGAS_WT')
cGAS_65[["Sex"]] = c('F')
cGAS_65[["Batch"]] = c('Batch_4')
rm(cGAS_65.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
cGAS_65[["percent.mt"]] <- PercentageFeatureSet(object = cGAS_65, pattern = "^mt-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/jcu4001/JCUPs19cGas/TestRunJU1/Joe_cGas_3rd/Joe_26_66/outs')
sc = autoEstCont(sc)
cGAS_66.counts = adjustCounts(sc)
cGAS_66 <- CreateSeuratObject(counts = cGAS_66.counts, project = "66", min.cells = 3, min.features = 620)
cGAS_66[["Genotype"]] = c('P301S tau: +; cGAS: +/-')
cGAS_66[["Sample_Name"]] = c('PS19_cGAS_Hetero_6')
cGAS_66[["Condition_1"]] = c('PS19_cGAS_Hetero_F')
cGAS_66[["Condition_2"]] = c('PS19_cGAS_Hetero')
cGAS_66[["Condition_3"]] = c('PS19_F')
cGAS_66[["Condition_4"]] = c('cGAS_Hetero_F')
cGAS_66[["PS19"]] = c('PS19')
cGAS_66[["cGAS"]] = c('cGAS_Hetero')
cGAS_66[["Sex"]] = c('F')
cGAS_66[["Batch"]] = c('Batch_4')
rm(cGAS_66.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
cGAS_66[["percent.mt"]] <- PercentageFeatureSet(object = cGAS_66, pattern = "^mt-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/scratch/jcu4001/JCUPs19cGas/TestRunJU1/Joe_cGas_3rd/Joe_28_77/outs')
sc = autoEstCont(sc)
cGAS_77.counts = adjustCounts(sc)
cGAS_77 <- CreateSeuratObject(counts = cGAS_77.counts, project = "77", min.cells = 3, min.features = 620)
cGAS_77[["Genotype"]] = c('P301S tau: +; cGAS: +/+')
cGAS_77[["Sample_Name"]] = c('PS19_cGAS_WT_6')
cGAS_77[["Condition_1"]] = c('PS19_cGAS_WT_M')
cGAS_77[["Condition_2"]] = c('PS19_cGAS_WT')
cGAS_77[["Condition_3"]] = c('PS19_M')
cGAS_77[["Condition_4"]] = c('cGAS_WT_M')
cGAS_77[["PS19"]] = c('PS19')
cGAS_77[["cGAS"]] = c('cGAS_Homo')
cGAS_77[["Sex"]] = c('M')
cGAS_77[["Batch"]] = c('Batch_4')
rm(cGAS_77.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
cGAS_77[["percent.mt"]] <- PercentageFeatureSet(object = cGAS_77, pattern = "^mt-") #recognize mitochondrial transcripts



###############################################################################################
all <- cGAS_25
pdf("cGAS_25_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("cGAS_25_FeatureScatter.pdf", width=12, height=4)
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
pdf("cGAS_25_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_25_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("cGAS_25_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "cGAS_25_cell_counts_Condition.csv")
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8385) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_511", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("cGAS_25_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_25_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_511" #visualizing the singlet vs doublet cells
pdf("cGAS_25_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"cGAS_25_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"cGAS_25_singlets.rds")
singlets<-readRDS("cGAS_25_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("cGAS_25_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("cGAS_25_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"cGAS_25_singlets_PCA.rds")
###############################################################################################
all <- cGAS_23
pdf("cGAS_23_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("cGAS_23_FeatureScatter.pdf", width=12, height=4)
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
pdf("cGAS_23_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_23_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("cGAS_23_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "cGAS_23_cell_counts_Condition.csv")
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.054*7246) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_391", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("cGAS_23_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_23_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_391" #visualizing the singlet vs doublet cells
pdf("cGAS_23_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"cGAS_23_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"cGAS_23_singlets.rds")
singlets<-readRDS("cGAS_23_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("cGAS_23_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("cGAS_23_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"cGAS_23_singlets_PCA.rds")
###############################################################################################
###############################################################################################
all <- cGAS_18
pdf("cGAS_18_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("cGAS_18_FeatureScatter.pdf", width=12, height=4)
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
pdf("cGAS_18_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_18_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("cGAS_18_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "cGAS_18_cell_counts_Condition.csv")
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8418) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.26, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.26, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.26_513", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("cGAS_18_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_18_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.26_513" #visualizing the singlet vs doublet cells
pdf("cGAS_18_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"cGAS_18_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"cGAS_18_singlets.rds")
singlets<-readRDS("cGAS_18_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("cGAS_18_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("cGAS_18_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"cGAS_18_singlets_PCA.rds")
###############################################################################################
###############################################################################################
all <- cGAS_20
pdf("cGAS_20_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("cGAS_20_FeatureScatter.pdf", width=12, height=4)
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
pdf("cGAS_20_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_20_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("cGAS_20_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "cGAS_20_cell_counts_Condition.csv")
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.054*7903) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.14, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.14, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.14_427", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("cGAS_20_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_20_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.14_427" #visualizing the singlet vs doublet cells
pdf("cGAS_20_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"cGAS_20_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"cGAS_20_singlets.rds")
singlets<-readRDS("cGAS_20_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("cGAS_20_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("cGAS_20_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"cGAS_20_singlets_PCA.rds")
###############################################################################################

###############################################################################################
all <- cGAS_55
pdf("cGAS_55_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("cGAS_55_FeatureScatter.pdf", width=12, height=4)
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
pdf("cGAS_55_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_55_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("cGAS_55_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8631) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_526", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("cGAS_55_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_55_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_526" #visualizing the singlet vs doublet cells
pdf("cGAS_55_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"cGAS_55_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"cGAS_55_singlets.rds")
singlets<-readRDS("cGAS_55_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("cGAS_55_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("cGAS_55_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"cGAS_55_singlets_PCA.rds")

###############################################################################################
###############################################################################################
all <- cGAS_61
pdf("cGAS_61_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("cGAS_61_FeatureScatter.pdf", width=12, height=4)
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
pdf("cGAS_61_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_61_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("cGAS_61_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.054*7770) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_420", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("cGAS_61_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_61_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_420" #visualizing the singlet vs doublet cells
pdf("cGAS_61_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"cGAS_61_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"cGAS_61_singlets.rds")
singlets<-readRDS("cGAS_61_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("cGAS_61_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("cGAS_61_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"cGAS_61_singlets_PCA.rds")

###############################################################################################
###############################################################################################
all <- cGAS_21
pdf("cGAS_21_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("cGAS_21_FeatureScatter.pdf", width=12, height=4)
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
pdf("cGAS_21_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_21_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("cGAS_21_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.069*9281) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.17, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.17, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.17_640", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("cGAS_21_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_21_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.17_640" #visualizing the singlet vs doublet cells
pdf("cGAS_21_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"cGAS_21_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"cGAS_21_singlets.rds")
singlets<-readRDS("cGAS_21_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("cGAS_21_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("cGAS_21_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"cGAS_21_singlets_PCA.rds")

###############################################################################################
###############################################################################################
all <- cGAS_22
pdf("cGAS_22_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("cGAS_22_FeatureScatter.pdf", width=12, height=4)
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
pdf("cGAS_22_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_22_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("cGAS_22_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.069*9293) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.01, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.01, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.01_641", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("cGAS_22_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_22_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.01_641" #visualizing the singlet vs doublet cells
pdf("cGAS_22_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"cGAS_22_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"cGAS_22_singlets.rds")
singlets<-readRDS("cGAS_22_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("cGAS_22_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("cGAS_22_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"cGAS_22_singlets_PCA.rds")

###############################################################################################
###############################################################################################
all <- cGAS_19
pdf("cGAS_19_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("cGAS_19_FeatureScatter.pdf", width=12, height=4)
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
pdf("cGAS_19_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_19_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("cGAS_19_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.069*9725) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.01, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.01, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.01_671", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("cGAS_19_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_19_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.01_671" #visualizing the singlet vs doublet cells
pdf("cGAS_19_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"cGAS_19_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"cGAS_19_singlets.rds")
singlets<-readRDS("cGAS_19_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("cGAS_19_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("cGAS_19_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"cGAS_19_singlets_PCA.rds")

###############################################################################################
###############################################################################################
all <- cGAS_30
pdf("cGAS_30_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("cGAS_30_FeatureScatter.pdf", width=12, height=4)
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
pdf("cGAS_30_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_30_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("cGAS_30_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8319) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.19, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.19, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.19_507", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("cGAS_30_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_30_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.19_507" #visualizing the singlet vs doublet cells
pdf("cGAS_30_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"cGAS_30_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"cGAS_30_singlets.rds")
singlets<-readRDS("cGAS_30_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("cGAS_30_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("cGAS_30_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"cGAS_30_singlets_PCA.rds")

###############################################################################################
###############################################################################################
all <- cGAS_29
pdf("cGAS_29_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("cGAS_29_FeatureScatter.pdf", width=12, height=4)
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
pdf("cGAS_29_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_29_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("cGAS_29_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8563) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.13, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.13, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.13_522", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("cGAS_29_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_29_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.13_522" #visualizing the singlet vs doublet cells
pdf("cGAS_29_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"cGAS_29_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"cGAS_29_singlets.rds")
singlets<-readRDS("cGAS_29_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("cGAS_29_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("cGAS_29_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"cGAS_29_singlets_PCA.rds")

###############################################################################################
###############################################################################################
all <- cGAS_63
pdf("cGAS_63_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("cGAS_63_FeatureScatter.pdf", width=12, height=4)
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
pdf("cGAS_63_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_63_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("cGAS_63_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.069*9335) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.02, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.02, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.02_644", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("cGAS_63_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_63_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.02_644" #visualizing the singlet vs doublet cells
pdf("cGAS_63_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"cGAS_63_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"cGAS_63_singlets.rds")
singlets<-readRDS("cGAS_63_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("cGAS_63_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("cGAS_63_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"cGAS_63_singlets_PCA.rds")

###############################################################################################
###############################################################################################
all <- cGAS_36
pdf("cGAS_36_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("cGAS_36_FeatureScatter.pdf", width=12, height=4)
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
pdf("cGAS_36_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_36_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("cGAS_36_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.076*10297) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.3, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.3, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.3_783", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("cGAS_36_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_36_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.3_783" #visualizing the singlet vs doublet cells
pdf("cGAS_36_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"cGAS_36_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"cGAS_36_singlets.rds")
singlets<-readRDS("cGAS_36_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("cGAS_36_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("cGAS_36_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"cGAS_36_singlets_PCA.rds")

###############################################################################################
###############################################################################################
all <- cGAS_43
pdf("cGAS_43_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("cGAS_43_FeatureScatter.pdf", width=12, height=4)
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
pdf("cGAS_43_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_43_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("cGAS_43_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.054*7906) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.08, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.08, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.08_427", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("cGAS_43_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_43_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.08_427" #visualizing the singlet vs doublet cells
pdf("cGAS_43_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"cGAS_43_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"cGAS_43_singlets.rds")
singlets<-readRDS("cGAS_43_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("cGAS_43_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("cGAS_43_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"cGAS_43_singlets_PCA.rds")

###############################################################################################
###############################################################################################
all <- cGAS_33
pdf("cGAS_33_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("cGAS_33_FeatureScatter.pdf", width=12, height=4)
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
pdf("cGAS_33_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_33_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("cGAS_33_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.054*7940) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.01, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.01, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.01_429", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("cGAS_33_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_33_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.01_429" #visualizing the singlet vs doublet cells
pdf("cGAS_33_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"cGAS_33_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"cGAS_33_singlets.rds")
singlets<-readRDS("cGAS_33_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("cGAS_33_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("cGAS_33_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"cGAS_33_singlets_PCA.rds")

###############################################################################################
###############################################################################################
all <- cGAS_59
pdf("cGAS_59_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("cGAS_59_FeatureScatter.pdf", width=12, height=4)
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
pdf("cGAS_59_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_59_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("cGAS_59_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.054*7765) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.21, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.21, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.21_419", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("cGAS_59_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_59_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.21_419" #visualizing the singlet vs doublet cells
pdf("cGAS_59_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"cGAS_59_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"cGAS_59_singlets.rds")
singlets<-readRDS("cGAS_59_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("cGAS_59_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("cGAS_59_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"cGAS_59_singlets_PCA.rds")

###############################################################################################
###############################################################################################
all <- cGAS_42
pdf("cGAS_42_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("cGAS_42_FeatureScatter.pdf", width=12, height=4)
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
pdf("cGAS_42_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_42_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("cGAS_42_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.054*7709) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_416", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("cGAS_42_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_42_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_416" #visualizing the singlet vs doublet cells
pdf("cGAS_42_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"cGAS_42_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"cGAS_42_singlets.rds")
singlets<-readRDS("cGAS_42_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("cGAS_42_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("cGAS_42_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"cGAS_42_singlets_PCA.rds")

###############################################################################################
###############################################################################################
all <- cGAS_74
pdf("cGAS_74_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("cGAS_74_FeatureScatter.pdf", width=12, height=4)
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
pdf("cGAS_74_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_74_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("cGAS_74_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8317) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.2, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.2, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.2_507", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("cGAS_74_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_74_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.2_507" #visualizing the singlet vs doublet cells
pdf("cGAS_74_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"cGAS_74_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"cGAS_74_singlets.rds")
singlets<-readRDS("cGAS_74_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("cGAS_74_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("cGAS_74_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"cGAS_74_singlets_PCA.rds")

###############################################################################################
###############################################################################################
all <- cGAS_45
pdf("cGAS_45_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("cGAS_45_FeatureScatter.pdf", width=12, height=4)
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
pdf("cGAS_45_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_45_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("cGAS_45_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8612) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.01, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.01, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.01_525", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("cGAS_45_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_45_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.01_525" #visualizing the singlet vs doublet cells
pdf("cGAS_45_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"cGAS_45_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"cGAS_45_singlets.rds")
singlets<-readRDS("cGAS_45_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("cGAS_45_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("cGAS_45_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"cGAS_45_singlets_PCA.rds")

###############################################################################################
###############################################################################################
all <- cGAS_46
pdf("cGAS_46_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("cGAS_46_FeatureScatter.pdf", width=12, height=4)
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
pdf("cGAS_46_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_46_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("cGAS_46_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.054*7849) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_424", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("cGAS_46_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_46_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_424" #visualizing the singlet vs doublet cells
pdf("cGAS_46_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"cGAS_46_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"cGAS_46_singlets.rds")
singlets<-readRDS("cGAS_46_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("cGAS_46_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("cGAS_46_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"cGAS_46_singlets_PCA.rds")

###############################################################################################

all <- cGAS_27
pdf("cGAS_27_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("cGAS_27_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 6200)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("cGAS_27_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_27_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("cGAS_27_ggplot_pK.pdf", width=57, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "cGAS_27_cell_counts_Condition.csv")
saveRDS(all, file = 'cGAS_27.rds')

###############################################################################################
all <- cGAS_41
pdf("cGAS_41_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("cGAS_41_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 6200)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("cGAS_41_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_41_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("cGAS_41_ggplot_pK.pdf", width=57, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "cGAS_41_cell_counts_Condition.csv")
saveRDS(all, file = 'cGAS_41.rds')
###############################################################################################
###############################################################################################
all <- cGAS_57
pdf("cGAS_57_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("cGAS_57_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 6200)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("cGAS_57_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_57_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("cGAS_57_ggplot_pK.pdf", width=57, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "cGAS_57_cell_counts_Condition.csv")
saveRDS(all, file = 'cGAS_57.rds')
###############################################################################################
###############################################################################################
all <- cGAS_62
pdf("cGAS_62_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("cGAS_62_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 6200)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("cGAS_62_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_62_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("cGAS_62_ggplot_pK.pdf", width=57, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "cGAS_62_cell_counts_Condition.csv")
saveRDS(all, file = 'cGAS_62.rds')
###############################################################################################
###############################################################################################
all <- cGAS_65
pdf("cGAS_65_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("cGAS_65_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 6200)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("cGAS_65_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_65_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("cGAS_65_ggplot_pK.pdf", width=57, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "cGAS_65_cell_counts_Condition.csv")
saveRDS(all, file = 'cGAS_65.rds')
###############################################################################################
###############################################################################################
all <- cGAS_66
pdf("cGAS_66_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("cGAS_66_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 6200)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("cGAS_66_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_66_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("cGAS_66_ggplot_pK.pdf", width=57, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "cGAS_66_cell_counts_Condition.csv")
saveRDS(all, file = 'cGAS_66.rds')
###############################################################################################

###############################################################################################
all <- cGAS_77
pdf("cGAS_77_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("cGAS_77_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 6200)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("cGAS_77_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_77_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("cGAS_77_ggplot_pK.pdf", width=57, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
write.csv(table(all$seurat_clusters), "cGAS_77_cell_counts_Condition.csv")
saveRDS(all, file = 'cGAS_77.rds')
###############################################################################################
###############################################################################################

all <- readRDS('cGAS_27.rds')
a = length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.039*a) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
#visualizing clusters and multiplet cells====
pdf("cGAS_27_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_27_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

saveRDS(all,"cGAS_27_after_doublet_detection.rds")

###############################################################################################
###############################################################################################
all <- readRDS('cGAS_41.rds')

b = length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.046*b) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.01, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
#visualizing clusters and multiplet cells====
pdf("cGAS_41_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_41_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
saveRDS(all,"cGAS_41_after_doublet_detection.rds")

###############################################################################################
###############################################################################################
all <- readRDS('cGAS_57.rds')

c = length(all@meta.data$seurat_clusters)
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.039*c) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
#visualizing clusters and multiplet cells====
pdf("cGAS_57_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_57_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
saveRDS(all,"cGAS_57_after_doublet_detection.rds")

###############################################################################################
###############################################################################################
all <- readRDS('cGAS_62.rds')

d = length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.054*d) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
#visualizing clusters and multiplet cells====
pdf("cGAS_62_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_62_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
saveRDS(all,"cGAS_62_after_doublet_detection.rds")

###############################################################################################
###############################################################################################
all <- readRDS('cGAS_65.rds')

e = length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.069*e) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
#visualizing clusters and multiplet cells====
pdf("cGAS_65_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_65_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
saveRDS(all,"cGAS_65_after_doublet_detection.rds")

###############################################################################################
###############################################################################################
all <- readRDS('cGAS_66.rds')

f = length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.046*f) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.01, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
#visualizing clusters and multiplet cells====
pdf("cGAS_66_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_66_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_420" #visualizing the singlet vs doublet cells
saveRDS(all,"cGAS_66_after_doublet_detection.rds")

###############################################################################################

##########################################################################

all <- readRDS('cGAS_77.rds')

h = length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.046*h) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
#visualizing clusters and multiplet cells====
pdf("cGAS_77_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("cGAS_77_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
saveRDS(all,"cGAS_77_after_doublet_detection.rds")

###############################################################################################
all <- readRDS('cGAS_27_after_doublet_detection.rds')
Idents(object = all) <- "DF.classifications_0.25_0.005_311" #visualizing the singlet vs doublet cells
pdf("cGAS_27_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"cGAS_27_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"cGAS_27_singlets.rds")
singlets<-readRDS("cGAS_27_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("cGAS_27_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("cGAS_27_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"cGAS_27_singlets_PCA.rds")

################
all <- readRDS('cGAS_41_after_doublet_detection.rds')
Idents(object = all) <- "DF.classifications_0.25_0.01_455" #visualizing the singlet vs doublet cells
pdf("cGAS_41_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"cGAS_41_after_doublet_detection.rds")

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"cGAS_41_singlets.rds")
singlets<-readRDS("cGAS_41_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("cGAS_41_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("cGAS_41_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"cGAS_41_singlets_PCA.rds")

##########
all <- readRDS('cGAS_57_after_doublet_detection.rds')
Idents(object = all) <- "DF.classifications_0.25_0.005_294" #visualizing the singlet vs doublet cells
pdf("cGAS_57_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"cGAS_57_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"cGAS_57_singlets.rds")
singlets<-readRDS("cGAS_57_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("cGAS_57_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("cGAS_57_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"cGAS_57_singlets_PCA.rds")

###########
all <- readRDS('cGAS_62_after_doublet_detection.rds')
Idents(object = all) <- "DF.classifications_0.25_0.005_617" #visualizing the singlet vs doublet cells
pdf("cGAS_62_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"cGAS_62_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"cGAS_62_singlets.rds")
singlets<-readRDS("cGAS_62_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("cGAS_62_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("cGAS_62_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"cGAS_62_singlets_PCA.rds")

#####
all <- readRDS('cGAS_65_after_doublet_detection.rds')
Idents(object = all) <- "DF.classifications_0.25_0.005_993" #visualizing the singlet vs doublet cells
pdf("cGAS_65_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"cGAS_65_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"cGAS_65_singlets.rds")
singlets<-readRDS("cGAS_65_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("cGAS_65_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("cGAS_65_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"cGAS_65_singlets_PCA.rds")
############

all <- readRDS('cGAS_66_after_doublet_detection.rds')
Idents(object = all) <- "DF.classifications_0.25_0.01_497" #visualizing the singlet vs doublet cells
pdf("cGAS_66_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"cGAS_66_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"cGAS_66_singlets.rds")
singlets<-readRDS("cGAS_66_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("cGAS_66_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("cGAS_66_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"cGAS_66_singlets_PCA.rds")

#####################
all <- readRDS('cGAS_77_after_doublet_detection.rds')
Idents(object = all) <- "DF.classifications_0.25_0.005_412" #visualizing the singlet vs doublet cells
pdf("cGAS_77_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"cGAS_77_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"cGAS_77_singlets.rds")
singlets<-readRDS("cGAS_77_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("cGAS_77_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("cGAS_77_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"cGAS_77_singlets_PCA.rds")








