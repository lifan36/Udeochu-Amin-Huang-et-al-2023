
DefaultAssay(Cluster_AST) <- 'integrated'
all.gASTes <- rownames(Cluster_AST)
Cluster_AST <- ScaleData(Cluster_AST, features = all.gASTes)
Cluster_AST <- FASTdVariableFeatures(object = Cluster_AST)
Cluster_AST <- RunPCA(Cluster_AST, features = VariableFeatures(object = Cluster_AST))

DimPlot(Cluster_AST, reduction = "pca")
ElbowPlot(Cluster_AST)

Cluster_AST <- FASTdNeighbors(Cluster_AST, dims = 1:10)
Cluster_AST <- FASTdClusters(Cluster_AST, resolution = 0.1)
Cluster_AST <- RunUMAP(Cluster_AST, dims = 1:10)

Idents(Cluster_AST) <- 'Genotype'
PS19KOvsPS19WT_DE <- FindMarkers(Cluster_AST, ident.1 = "P301S tau: +; cGAS: -/-", ident.2 = "P301S tau: +; cGAS: +/+", logfc.threshold = 0,
                                 test.use = "MAST", assay ='RNA')
write.csv(PS19KOvsPS19WT_DE, file = 'PS19KOvsPS19WT_DE_MAST_RNAassay_Astrocyte.csv')

PS19WTvsWT_DE <- FindMarkers(Cluster_AST, ident.1 = "P301S tau: +; cGAS: +/+", ident.2 = "P301S tau: -; cGAS: +/+", logfc.threshold = 0,
                                 test.use = "MAST", assay ='RNA')
write.csv(PS19WTvsWT_DE, file = 'PS19WTvsWT_DE_MAST_RNAassay_Astrocyte.csv')

# Volcano plot of marker genes =====
data <- PS19KOvsPS19WT_DE
data$colours <- c("NC")
data$colours[data$avg_log2FC >= 0.1 & data$p_val_adj <= 0.05] <- c("UP")
data$colours[data$avg_log2FC <= -0.1 & data$p_val_adj <= 0.05] <- c("DN")

# Selected genes to highlight
genes_select_up<- c("Akap6",'Kcnq3','Cntn1','Unc13c','Fgr13','Gabbr2','Clu','Apoe','Ldlr','App','mt-Co1','mt-Co3','mt-Atp6')
genes_to_plot_up <- data[rownames(data) %in% genes_select_up, ]
genes_to_plot_up$Cluster <- "up"

genes_select_down <- c("Csmd1", "Cacna1a","Robo1","Robo2","Slit2",'Kcnd2','Gfap','Ddx60','Trim30a','Stat1','Rnf213')
genes_to_plot_down <- data[rownames(data) %in% genes_select_down, ]
genes_to_plot_down$Cluster <- c("down")

genes_to_plot <- rbind(genes_to_plot_up, genes_to_plot_down)

# Set color palette
my_color_1 <- c("#D7301F","Grey", "#2B8CBE")

# Plot volcano plot

pdf("Extended Fig4.g.pdf", width=7, height=4.5)
ggplot() + 
  geom_point(data=data, aes(x=avg_log2FC, y=-log10(p_val_adj), colour=colours),
             shape=19, alpha=1, size=0.5) +
  scale_color_manual(values = my_color_1,
                     name="DEGs",
                     breaks=rev(names(table(data$colours))),
                     labels=rev(names(table(data$colours)))) +
  geom_point(data=genes_to_plot,
             aes(x=avg_log2FC, y=-log10(p_val_adj)),
             shape=19, alpha=1, size=1) + geom_text_repel(data=genes_to_plot,
                                                            aes(x=avg_log2FC, y=-log10(p_val_adj), label = rownames(genes_to_plot)), 
                                                            max.overlaps = 20, color="black",size = 3, box.padding = 0.5,
                                                            point.padding = 0.5, segment.size=0.25, segment.colour="black") +
  ylab("-Log10[FDR]") + xlab("Log2FC") +
  ggtitle("Cgas-/- P301S vs P301S AST")+
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

##### Extended Data Fig4.j
Idents(Cluster_AST) <- 'Genotype'
Cluster_AST_sub <- subset(Cluster_AST, subset= Genotype == 'P301S tau: +; cGAS: +/+'| Genotype =='P301S tau: +; cGAS: +/-'| Genotype == 'P301S tau: +; cGAS: -/-' | Genotype =='P301S tau: -; cGAS: +/+' )

pdf("Extended Data Fig4.j.pdf", width=8, height=4)
DotPlot(Cluster_AST_sub, features = c('Gfap','Ddx60','Rnf213','Stat1','Herc6','Trim30a'), assay = 'RNA') + RotatedAxis()+
  scale_colour_gradient2(low = "darkblue", mid = "white", high = "darkred")
dev.off()

### Extended Fig 4i upstream regulator ======

ipa<-read.csv("AST_upstream_final.csv",header=T)
ipa <- ipa[order(ipa$Activation.z.score), ]

ggplot(data=ipa, aes(x=reorder(Upstream.Regulator,-log(p.value.of.overlap)), y= -log(p.value.of.overlap), fill=Activation.z.score)) +
  geom_dotplot(position ="identity",binaxis='y', stackdir='center', dotsize = 1.6)+ coord_flip()+ scale_fill_continuous(low = 'blue', high='red')+
  theme(legend.position="top")+
  theme(axis.text.y = element_text(size = 15))+
  theme(axis.line = element_line(colour = "black",size = 0.5, linetype = "solid"),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  labs(x=' ', y ="-log(p.value.of.overlap)")

#### bubble plot
pdf("Extended Data Fig4.i.pdf", width=8, height=4)
ggplot(data=ipa, aes(x=reorder(Upstream.Regulator,-log(p.value.of.overlap)), y= -log(p.value.of.overlap), label=round(Activation.z.score, digits=1))) + 
  geom_point(stat='identity', aes(col=Predicted.Activation.State, size=abs(Activation.z.score)))  +
  scale_size_continuous(range = c(5, 8))+
  scale_color_manual(name="Predicted Activation State", 
                     labels = c("Activated", "Inhibited"), 
                     values = c("Activated"="red", "Inhibited"="blue")) + 
  geom_text(color="white", size=2.5) +
  facet_wrap(~Predicted.Activation.State, scales='free_y',ncol=2 )+
  labs(title="Astrocyte Upstream Regulator") +
  theme_bw()+
  theme(axis.text.y = element_text(size = 10))+
  coord_flip()
dev.off()

##### oligodendrocyte

DefaultAssay(Cluster_Oligo) <- 'integrated'
all.genes <- rownames(Cluster_Oligo)
Cluster_Oligo <- ScaleData(Cluster_Oligo, features = all.genes)
Cluster_Oligo <- FindVariableFeatures(object = Cluster_Oligo)
Cluster_Oligo <- RunPCA(Cluster_Oligo, features = VariableFeatures(object = Cluster_Oligo))

DimPlot(Cluster_Oligo, reduction = "pca")
ElbowPlot(Cluster_Oligo)

Cluster_Oligo <- FindNeighbors(Cluster_Oligo, dims = 1:10)
Cluster_Oligo <- FindClusters(Cluster_Oligo, resolution = 0.2)
Cluster_Oligo <- RunUMAP(Cluster_Oligo, dims = 1:10)

Idents(Cluster_Oligo) <- 'Genotype'
PS19KOvsPS19WT_DE <- FindMarkers(Cluster_Oligo, ident.1 = "P301S tau: +; cGAS: -/-", ident.2 = "P301S tau: +; cGAS: +/+", logfc.threshold = 0,
                                 test.use = "MAST", assay ='RNA')
write.csv(PS19KOvsPS19WT_DE, file = 'PS19KOvsPS19WT_DE_MAST_RNAassay_Oligodendrocyte.csv')

# Volcano plot of marker genes =====
data <- PS19KOvsPS19WT_DE
data$colours <- c("NC")
data$colours[data$avg_log2FC >= 0.1 & data$p_val_adj <= 0.05] <- c("UP")
data$colours[data$avg_log2FC <= -0.1 & data$p_val_adj <= 0.05] <- c("DN")

# Selected genes to highlight
genes_select_up<- c('Mylk','Zbtb16','Spock1','mt-Atp6','mt-Co3','mt-Co2','Mef2c','Gria2')
genes_to_plot_up <- data[rownames(data) %in% genes_select_up, ]
genes_to_plot_up$Cluster <- "up"

genes_select_down <- c('Gria1','Cacnb2','Grik1','Camk2a','Grin2a','Gabrb3','Camk2d','Grin2b')
genes_to_plot_down <- data[rownames(data) %in% genes_select_down, ]
genes_to_plot_down$Cluster <- c("down")

genes_to_plot <- rbind(genes_to_plot_up, genes_to_plot_down)

# Set color palette
my_color_1 <- c("#D7301F","Grey", "#2B8CBE")

# Plot volcano plot

pdf("Extended Fig4.h.pdf", width=7, height=4.5)
ggplot() + 
  geom_point(data=data, aes(x=avg_log2FC, y=-log10(p_val_adj), colour=colours),
             shape=19, alpha=1, size=0.5) +
  scale_color_manual(values = my_color_1,
                     name="DEGs",
                     breaks=rev(names(table(data$colours))),
                     labels=rev(names(table(data$colours)))) +
  geom_point(data=genes_to_plot,
             aes(x=avg_log2FC, y=-log10(p_val_adj)),
             shape=19, alpha=1, size=1) + geom_text_repel(data=genes_to_plot,
                                                          aes(x=avg_log2FC, y=-log10(p_val_adj), label = rownames(genes_to_plot)), 
                                                          max.overlaps = 20, color="black",size = 3, box.padding = 0.5,
                                                          point.padding = 0.5, segment.size=0.25, segment.colour="black") +
  ylab("-Log10[FDR]") + xlab("Log2FC") +
  ggtitle("Cgas-/- P301S vs P301S Oligo")+
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
  ylim(0,150)+
  scale_x_continuous(breaks=seq(-1, 1, 0.5), limits=c(-1, 1))
dev.off()
