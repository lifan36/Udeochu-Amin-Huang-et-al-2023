library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(reshape2)
library(MAST)

setwd()

Cluster_EN <- readRDS(file = 'EN_cGAS_reclustered.rds')
Cluster_IN <- readRDS(file = 'IN_cGAS_reclustered.rds')
Cluster_EN_sub <- subset(Cluster_EN, subset= Genotype =='P301S tau: -; cGAS: +/+'| Genotype == 'P301S tau: +; cGAS: -/-' | Genotype == 'P301S tau: +; cGAS: +/+')
saveRDS(Cluster_EN_sub, 'Cluster_EN_3genotypes.rds')
Cluster_IN_sub <- subset(Cluster_IN, subset= Genotype =='P301S tau: -; cGAS: +/+'| Genotype == 'P301S tau: +; cGAS: -/-' | Genotype == 'P301S tau: +; cGAS: +/+')
saveRDS(Cluster_IN_sub,'Cluster_IN_3genotypes.rds')


DefaultAssay(Cluster_EN_sub) <- 'RNA'
Idents(Cluster_EN_sub) <- 'Genotype'
pdf("Fig5.i.pdf", width=6.5, height=4)
DotPlot(Cluster_EN_sub, features = c('Mef2c','Vwa5b2','Sfxn5','Srcin1','Ptprt','Peg3','Spock1','Tshz2','Ncald','Rasgef1b')) + scale_colour_gradient2(low = "darkblue", mid = "white", high = "darkred")+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

DefaultAssay(Cluster_IN_sub) <- 'RNA'
Idents(Cluster_IN_sub) <- 'Genotype'
pdf("Fig5.j.pdf", width=6.5, height=4)
DotPlot(Cluster_IN_sub, features = c('Mef2c','Vwa5b2','Grin1','Edil3','Peg3','Ncald','Igsf3','R3hdm2'), 
        scale.min = 0) + scale_colour_gradient2(low = "darkblue", mid = "white", high = "darkred")+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()
