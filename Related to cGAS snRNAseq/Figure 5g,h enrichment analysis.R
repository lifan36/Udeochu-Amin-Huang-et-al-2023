
# Figure 5h ###################################################################
Mef2c  = read.csv('Mef2c Target genes.csv', header = F)
Mef2a  = read.csv('Mef2a Target genes.csv', header = F)
JunB  = read.csv('ms JUNB Targets.csv',header = F)
ARG = read.csv('ARG_Tyssowski2018.csv',header = F)
Activity_genes = read.csv('ActivityGenes_Lacar2016.csv',header = F)
FOSL2 = read.csv('FOSL2 Targets.csv',header = F)

PS19KOVSWT_EN = read.csv('PS19KOvsPS19WT_DE_MAST_RNAassay_ExcitatoryNeuron.csv', header = T)
PS19KOVSWT_IN = read.csv('PS19KOvsPS19WT_DE_MAST_RNAassay_InhibitoryNeuron.csv', header = T)

PS19KOVSWT_EN_deg = filter(PS19KOVSWT_EN, p_val_adj < 0.05&(avg_log2FC>=0.1|avg_log2FC<=-0.1))
PS19KOVSWT_IN_deg = filter(PS19KOVSWT_IN, p_val_adj < 0.05&(avg_log2FC>=0.1|avg_log2FC<=-0.1))

library(GeneOverlap)

go.obj_EN <- newGeneOverlap(PS19KOVSWT_EN_deg$X,Mef2c$V3, genome.size=21988)
go.obj_IN <- newGeneOverlap(PS19KOVSWT_IN_deg$X,Mef2c$V3, genome.size=21988)
go.obj_EN <- testGeneOverlap(go.obj_EN)
go.obj_IN <- testGeneOverlap(go.obj_IN)
print(go.obj_EN)
print(go.obj_IN)


#save DEGs that are Mef2c targets ---> supplemental table S7
write.csv(go.obj_EN@intersection, 'Mef2c target_ExN DEG.cvs')
write.csv(go.obj_IN@intersection, 'Mef2c target_IN DEG.cvs')

# Figure 5g ###################################################################
library(VennDiagram)

set3 <-Mef2c$V3
set1 <- PS19KOVSWT_EN_deg$X
set2 <- PS19KOVSWT_IN_deg$X

library(ggVennDiagram)
x <- list('set1'=Mef2c$V3,'set2'=PS19KOVSWT_EN_deg$X,'set3'=PS19KOVSWT_IN_deg$X)
Venn <- ggVennDiagram(x,label = "count",category.names = c("Mef2c targets",
                                                           "PS19 Cgas KO VS PS19 EN",
                                                           "PS19 Cgas KO VS PS19 IN"),
                      label_alpha = 0)
Venn + scale_x_continuous(expand = expansion(mult = .2))+ scale_fill_distiller(palette = "RdBu")

# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

# Chart
library(VennDiagram)
grid.newpage()
venn<-venn.diagram(
  x = list(set2,set3,set1),
  category.names = c("P301S_Cgas-/-vsP301S_EN","P301S_Cgas-/-vsP301S_IN","Mef2c targets"),
  filename = NULL,
  output=FALSE,
  show.plot=TRUE,
  
  height = 3600, 
  width = 1800 , 
  resolution = 300,
  #compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill =c('blue','red','yellow'),
  
  # Numbers
  cex = 1.5,
  # Set names
  cat.cex = 1,
  cat.pos = c(-27, 27,135),
  cat.dist = c(0.055, 0.055,0.085),
  cat.default.pos = "outer",
  cat.fontfamily = "sans",
)
grid.newpage()
grid.draw(venn)
library(grDevices)
pdf(file="venn.pdf")
grid.draw(venn)
dev.off()

