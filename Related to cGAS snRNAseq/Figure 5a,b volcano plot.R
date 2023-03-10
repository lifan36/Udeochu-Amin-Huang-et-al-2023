library(ggplot2)
library(Seurat)
library(ggrepel)
library(RColorBrewer)
library(dplyr)
library(cowplot)
library(reshape2)
library(MAST)
library(patchwork)

setwd()

#Figure 5a
data<-read.csv("PS19KOvsPS19WT_DE_MAST_RNAassay_ExcitatoryNeuron.csv", header = T)
#refer to supplementary table S7

# Volcano plot of marker genes =====
data$colours <- c("NC")
data$colours[data$avg_log2FC >= 0.1 & data$p_val_adj <= 0.05] <- c("UP")
data$colours[data$avg_log2FC <= -0.1 & data$p_val_adj <= 0.05] <- c("DN")

# Selected genes to highlight
genes_select_up<- c("Pdzrn3", "Mef2c", "Satb1", "Satb2","Nrg1",'Pcdh7','Pcdh15','Dpp10')
genes_to_plot_up <- data[data$X %in% genes_select_up, ]
genes_to_plot_up$Cluster <- "up"

genes_select_down <- c("Gria1", "Grin2a","Cacnb2","Cacna1c","Grm1")
genes_to_plot_down <- data[data$X %in% genes_select_down, ]
genes_to_plot_down$Cluster <- c("down")

genes_to_plot <- rbind(genes_to_plot_up, genes_to_plot_down)

# Set color palette
my_color_1 <- c("#D7301F","Grey", "#2B8CBE")

# Plot volcano plot

pdf("Fig5.a.pdf", width=9.5, height=4.5)
ggplot() + 
  geom_point(data=data, aes(x=avg_log2FC, y=-log10(p_val_adj), colour=colours),
             shape=19, alpha=1, size=0.5) +
  scale_color_manual(values = my_color_1,
                     name="DEGs",
                     breaks=rev(names(table(data$colours))),
                     labels=rev(names(table(data$colours)))) +
  geom_point(data=genes_to_plot,
             aes(x=avg_log2FC, y=-log10(p_val_adj)),
             shape=19, alpha=1, size=1.5) + geom_text_repel(data=genes_to_plot,
                  aes(x=avg_log2FC, y=-log10(p_val_adj), label = genes_to_plot$X), 
                  max.overlaps = 20, color="black",size = 5, box.padding = 0.5,
                  point.padding = 0.5, segment.size=0.25, segment.colour="black") +
  ylab("-Log10[FDR]") + xlab("Log2FC") +
  ggtitle("Cgas-/- P301S vs P301S EN")+
  theme_bw()+
  theme(panel.grid.major.x  = element_blank(),
        panel.grid.major.y  = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.text.x = element_text(colour = "black", size=15),
        axis.text.y = element_text(colour = "black", size=15),
        axis.title.x = element_text(colour = "black", size=15),
        axis.title.y = element_text(colour = "black", size=15),
        plot.title = element_text(size = 15)) +
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks=seq(-1, 1, 0.5), limits=c(-1, 1))
dev.off()

#Figure 5b
data<-read.csv("PS19KOvsPS19WT_DE_MAST_RNAassay_InhibitoryNeuron.csv")
#refer to supplementary table S7

# Volcano plot of marker genes =====
data$colours <- c("NC")
data$colours[data$avg_log2FC >= 0.1 & data$p_val_adj <= 0.05] <- c("UP")
data$colours[data$avg_log2FC <= -0.1 & data$p_val_adj <= 0.05] <- c("DN")

# Selected genes to highlight
genes_select_up<- c("Gabbr2", "Mef2c", "Rnf144b", "Kcnc1","Kcnc2",'Scn1a','Slc6a1','Erbb4')
genes_to_plot_up <- data[data$X %in% genes_select_up, ]
genes_to_plot_up$Cluster <- "up"

genes_select_down <- c("Ryr3", "Cacnb2","Cacna1e","Atp2b1")
genes_to_plot_down <- data[data$X %in% genes_select_down, ]
genes_to_plot_down$Cluster <- c("down")

genes_to_plot <- rbind(genes_to_plot_up, genes_to_plot_down)

# Set color palette
my_color_1 <- c("#D7301F","Grey", "#2B8CBE")

# Plot volcano plot

pdf("Fig5.b.pdf", width=10, height=4.5)
ggplot() + 
  geom_point(data=data, aes(x=avg_log2FC, y=-log10(p_val_adj), colour=colours),
             shape=19, alpha=1, size=0.5) +
  scale_color_manual(values = my_color_1,
                     name="DEGs",
                     breaks=rev(names(table(data$colours))),
                     labels=rev(names(table(data$colours)))) +
  geom_point(data=genes_to_plot,
             aes(x=avg_log2FC, y=-log10(p_val_adj)),
             shape=19, alpha=1, size=1) + 
  geom_text_repel(data=genes_to_plot,
                  aes(x=avg_log2FC, y=-log10(p_val_adj), label = genes_to_plot$X), 
                  max.overlaps = 20, color="black",size = 5, box.padding = 0.5,
                  point.padding = 0.5, segment.size=0.25, segment.colour="black") +
  ylab("-Log10[FDR]") + xlab("Log2FC") +
  ggtitle("Cgas-/- P301S vs P301S IN")+
  theme_bw()+
  theme(panel.grid.major.x  = element_blank(),
        panel.grid.major.y  = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        axis.text.x = element_text(colour = "black", size=15),
        axis.text.y = element_text(colour = "black", size=15),
        axis.title.x = element_text(colour = "black", size=15),
        axis.title.y = element_text(colour = "black", size=15),
        plot.title = element_text(size = 15)) +
  ylim(0,40)+
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks=seq(-1, 1, 0.5), limits=c(-1, 1))
dev.off()
