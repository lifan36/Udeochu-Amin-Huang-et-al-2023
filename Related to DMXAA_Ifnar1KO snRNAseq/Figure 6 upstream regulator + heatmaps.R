library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(reshape2)
library(MAST)

#### Figure 6a,c - refer to Figure 5a,b codes
# DEG list refer to supplementary table S9

#### Figure 6b

ipa<-read.csv("MG_IPA_Upstream_Reg.csv",header=T)
ipa <- ipa[order(ipa$Activation.z.score), ]

#bubble plot
ggplot(data=ipa, aes(x=reorder(Upstream.Regulator,-log(p.value.of.overlap)), y= -log(p.value.of.overlap), label=round(Activation.z.score, digits=1))) + 
  geom_point(stat='identity', aes(col=Predicted.Activation.State, size=abs(Activation.z.score)))  +
  scale_size_continuous(range = c(5, 8))+
  scale_color_manual(name="Predicted Activation State", 
                     labels = c("Activated", "Inhibited"), 
                     values = c("Activated"="red", "Inhibited"="blue")) + 
  geom_text(color="white", size=2.5) +
  facet_wrap(~Predicted.Activation.State, scales='free_y',ncol=2 )+
  labs(title="Microglia Upstream Regulator") +
  theme_bw()+
  theme(axis.text.y = element_text(size = 10))+
  coord_flip()

#### Figure 6d - refer to Figure 5g codes

#### Figure 6e

IFN <- read.csv('CGAS_DMXAA_IFN_DEG_logFC_curated.csv')
x <- PS19KOvsPS19WT_DE[rownames(PS19KOvsPS19WT_DE) %in% IFN$logFC.of.DEG, ]
mat <- as.matrix(IFN)
mat_num <- matrix(as.numeric(mat), ncol = ncol(mat))
mat_num <- mat_num[,2:3]
rownames(mat_num) <- mat[,1]
colnames(mat_num) <- colnames(mat[,2:3])

mypal = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256))

p=pheatmap::pheatmap(t(mat_num),
                   scale="none", clustering_method="ward.D", color = mypal,border_color = NA,
                   angle_col = 45, fontsize = 11, main='MG IFN genes', legend_labels = 'Average expression')
save_pheatmap_pdf <- function(x, filename, width=4, height=2.5) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(p, "Fig.6e MG.pdf")



Mef2c <- read.csv('CGAS_DMXAA_IFN_Mef2c_DEG_logFC_curated.csv')
mat <- as.matrix(Mef2c)
mat_num <- matrix(as.numeric(mat), ncol = ncol(mat))
mat_num <- mat_num[,2:3]
rownames(mat_num) <- mat[,1]
colnames(mat_num) <- colnames(mat[,2:3])

mypal = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256))


p=pheatmap::pheatmap(t(mat_num),
                   scale="none", clustering_method="ward.D", color = mypal,border_color = NA,
                   angle_col = 45, fontsize = 11, main='EN Mef2c genes', legend_labels = 'Average expression')
save_pheatmap_pdf <- function(x, filename, width=6, height=2.5) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(p, "Fig.6e.pdf")


#### Figure 6g

DMXAA_EN <- readRDS('EN_reclusted_res0.1.rds')
cluster.averages <- AverageExpression(DMXAA_EN, return.seurat = FALSE)
cluster.averages_mat <- cluster.averages$RNA
Mef2c_targets <- read.csv('IFNARKO_DMXAA_Mef2c_Targets_Expression_Curated.csv')

mypal = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256))

pdf("4geno_DMXAA_EN_Mef2c_heat.pdf", width=10, height=8)
pheatmap::pheatmap(t(cluster.averages_mat),
                   scale="none", clustering_method="ward.D", color = mypal,border_color = NA,
                   angle_col = 45, fontsize = 11, main='EN Mef2c genes', legend_labels = 'Average expression')
dev.off()

