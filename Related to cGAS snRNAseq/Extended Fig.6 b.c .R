
EN_Marker <- c('Prox1','Pdzd2','Stxbp6',
               'Mpped1','Wfs1','Dcn','Fibcd1','Pex5l','Pou3f1','Satb2',
               'Cacng5','Gpr12',
               'Cpne4','Ly6e','Npy2r','Grik4','Elavl2')

DefaultAssay(Cluster_EN) <- 'RNA'
Idents(Cluster_EN) <- 'seurat_clusters'
pdf("Extended Data Fig.6.b.pdf", width=6.5, height=6.5)
DotPlot(Cluster_EN, features = EN_Marker) + 
  scale_colour_gradient2(low = "darkblue", mid = "white", high = "darkred")+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))
dev.off()

EN_DEG <- c('Pdzrn3','Dpp10','Efna5','Mef2c','Lrrc4c','Satb2','Satb1','Nrg1','Pcdh15','Pcdh7')
pdf("Extended Data Fig.6.c.pdf", width=6.5, height=6.5)
DotPlot(Cluster_EN, features = EN_DEG) + 
  scale_colour_gradient2(low = "darkblue", mid = "white", high = "darkred")+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))
dev.off()
