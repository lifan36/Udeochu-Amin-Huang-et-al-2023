library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(reshape2)
library(MAST)
library(EnhancedVolcano)

### Figure 5k 
Mathys <- readRDS('Mathys2019_processedHuman.rds')
Mathys_MG <- subset(Mathys2019_processedHuman, subset = broad.cell.type=='Mic')
Mathys_EN <- subset(Mathys2019_processedHuman, subset = broad.cell.type=='Ex')
Mathys_IN <- subset(Mathys2019_processedHuman, subset = broad.cell.type=='In')

Sayed_EN <- readRDS('human_all_54_EN.rds')
Sayed_IN <- readRDS('human_all_54_IN.rds')
Sayed_MG <- readRDS('human_all_54_MG.rds')

#get average expression of miroglial IRF3, RNF213 and neuronal MEF2C in each human sample
Idents(Mathys_MG) <- 'projid'
cluster.averages <- AverageExpression(Mathys_MG, return.seurat = FALSE)
cluster.averages_mat <- cluster.averages$RNA
IFN <- c('IRF3','RNF213')
cluster.averages_mat <- cluster.averages_mat[rownames(cluster.averages_mat) %in% IFN,]
write.csv(cluster.averages_mat, 'Mathys_MG_IFN.csv')

Idents(Mathys_EN) <- 'projid'
cluster.averages <- AverageExpression(Mathys_EN, return.seurat = FALSE)
cluster.averages_mat <- cluster.averages$RNA
cluster.averages_mat <- cluster.averages_mat[rownames(cluster.averages_mat) %in% 'MEF2C',]
write.csv(cluster.averages_mat, 'Mathys_EN_MEF2C.csv')

Idents(Mathys_IN) <- 'projid'
cluster.averages <- AverageExpression(Mathys_IN, return.seurat = FALSE)
cluster.averages_mat <- cluster.averages$RNA
cluster.averages_mat <- cluster.averages_mat[rownames(cluster.averages_mat) %in% 'MEF2C',]
write.csv(cluster.averages_mat, 'Mathys_IN_MEF2C.csv')

Idents(Sayed_MG) <- 'orig.ident'
cluster.averages <- AverageExpression(Sayed_MG, return.seurat = FALSE)
cluster.averages_mat <- cluster.averages$RNA
cluster.averages_mat <- cluster.averages_mat[rownames(cluster.averages_mat) %in% IFN,]
write.csv(cluster.averages_mat, 'Sayed_MG_IFN.csv')

Idents(Sayed_EN) <- 'orig.ident'
cluster.averages <- AverageExpression(Sayed_EN, return.seurat = FALSE)
cluster.averages_mat <- cluster.averages$RNA
cluster.averages_mat <- cluster.averages_mat[rownames(cluster.averages_mat) %in% 'MEF2C',]
write.csv(cluster.averages_mat, 'Sayed_EN_MEF2C.csv')

Idents(Sayed_IN) <- 'orig.ident'
cluster.averages <- AverageExpression(Sayed_IN, return.seurat = FALSE)
cluster.averages_mat <- cluster.averages$RNA
cluster.averages_mat <- cluster.averages_mat[rownames(cluster.averages_mat) %in% 'MEF2C',]
write.csv(cluster.averages_mat, 'Sayed_IN_MEF2C.csv')

# linear regression plot
data <- read.csv('combined_human_MEF2C_IFN.csv')
p1 <- ggplot(data, aes(x=MEF2C_EN, y=RNF213)) + geom_point(size=0.1) +
  geom_smooth(method="lm", se=T)+ylab('RNF213')+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.line = element_line(colour = "black",size = 0.5, linetype = "solid"))

P2 <- ggplot(data, aes(x=MEF2C_EN, y=IRF3)) + geom_point(size=0.1) +
  geom_smooth(method="lm", se=T)+ylab('IRF3')+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.line = element_line(colour = "black",size = 0.5, linetype = "solid"))

p3 <- ggplot(data, aes(x=MEF2C_IN, y=RNF213)) + geom_point(size=0.1) +
  geom_smooth(method="lm", se=T)+ylab('RNF213')+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.line = element_line(colour = "black",size = 0.5, linetype = "solid"))

P4 <- ggplot(data, aes(x=MEF2C_IN, y=IRF3)) + geom_point(size=0.1) +
  geom_smooth(method="lm", se=T)+ylab('IRF3')+
  theme(panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.line = element_line(colour = "black",size = 0.5, linetype = "solid"))

