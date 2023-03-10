Mef2c_EN <- read.csv('Mef2c target_ExN DEG.cvs', header = T)
Idents(Cluster_EN_sub) <- 'Genotype'
cluster.averages <- AverageExpression(Cluster_EN_sub, return.seurat = FALSE)
cluster.averages_mat <- cluster.averages$RNA
cluster.averages_mat <- cluster.averages_mat[rownames(cluster.averages_mat) %in% Mef2c_EN$x,]
cluster.averages_mat <- as.matrix(cluster.averages_mat)
colnames(cluster.averages_mat ) <- c('P301S tau: -; cGAS: +/+','P301S tau: +; cGAS: -/-','P301S tau: +; cGAS: +/+')

mypal = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256))


p=pheatmap::pheatmap(t(cluster.averages_mat),
                   scale="column", clustering_method="ward.D", color = mypal,border_color = NA,
                   angle_col = 90, fontsize = 11, main='Mef2cTarget genes_EN', legend_labels = 'Average expression',cluster_rows=FALSE)
save_pheatmap_pdf <- function(x, filename, width=32, height=3) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(p, "Extended Data Fig7.a test.pdf")

Mef2c_IN <- read.csv('Mef2c target_IN DEG.cvs', header = T)
Idents(Cluster_IN_sub) <- 'Genotype'
cluster.averages <- AverageExpression(Cluster_IN_sub, return.seurat = FALSE)
cluster.averages_mat <- cluster.averages$RNA
cluster.averages_mat <- cluster.averages_mat[rownames(cluster.averages_mat) %in% Mef2c_IN$x,]
cluster.averages_mat <- as.matrix(cluster.averages_mat)
colnames(cluster.averages_mat ) <- c('P301S tau: -; cGAS: +/+','P301S tau: +; cGAS: -/-','P301S tau: +; cGAS: +/+')

mypal = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256))

p=pheatmap::pheatmap(t(cluster.averages_mat),
                     scale="column", clustering_method="ward.D", color = mypal,border_color = NA,
                     angle_col = 90, fontsize = 11, main='Mef2cTarget genes_IN', legend_labels = 'Average expression',cluster_rows=FALSE)
save_pheatmap_pdf <- function(x, filename, width=18, height=3) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(p, "Extended Data Fig7.b test.pdf")
