library(monocle3)
Cluster_MG = readRDS('Cluster_MGreclusetered_03072023.rds')
gene_annotation <- as.data.frame(rownames(Cluster_MG@reductions[["pca"]]@feature.loadings), row.names = rownames(Cluster_MG@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"

# part two, cell information

cell_metadata <- as.data.frame(Cluster_MG@assays[["RNA"]]@counts@Dimnames[[2]], row.names = Cluster_MG@assays[["RNA"]]@counts@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"

# part three, counts sparse matrix

New_matrix <- Cluster_MG@assays[["RNA"]]@counts
New_matrix <- New_matrix[rownames(Cluster_MG@reductions[["pca"]]@feature.loadings), ]
expression_matrix <- New_matrix
colnames(expression_matrix) = expression_matrix@Dimnames[[2]]
# or use this to name colnames(expression_matrix) = row.names(Genotype_metadata)

Genotype_metadata <- as.data.frame(Cluster_MG@meta.data$Genotype, row.names = Cluster_MG@assays$RNA@counts@Dimnames[[2]])
colnames(Genotype_metadata) <- "Genotype"
Genotype_metadata$barcode = row.names(Genotype_metadata)

### Construct the basic cds object

cds_from_seurat <- new_cell_data_set(expression_matrix,
                                     cell_metadata = Genotype_metadata,
                                     gene_metadata = gene_annotation)

#reduce dimension required before learning object
cds_from_seurat <- preprocess_cds(cds_from_seurat, num_dim = 50)
cds_from_seurat <- reduce_dimension(cds_from_seurat, preprocess_method = "PCA")
plot_cells(cds_from_seurat,label_cell_groups=FALSE)
#note: cds_from_seurat structure not adequeate at this point. RowData only contains gene_short_name but not gene id
#also colData only has size factor, no genotype

cds_from_seurat = cluster_cells(cds_from_seurat, resolution=1e-4)

plot_cells(cds_from_seurat,label_cell_groups=FALSE, show_trajectory_graph = T,
           label_leaves=FALSE,
           label_branch_points=FALSE, color_cells_by = 'cluster') + facet_wrap(~Genotype)

cds_from_seurat <- learn_graph(cds_from_seurat)
cds_from_seurat <- order_cells(cds_from_seurat)
#chose branch point within homeostatic cluster
p = plot_cells(cds_from_seurat,
               color_cells_by = "pseudotime",
               label_cell_groups=FALSE,
               label_leaves=FALSE,
               label_branch_points=FALSE,
               graph_label_size=1.5)+facet_wrap(~Genotype, ncol = 4)
mypal = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256))

#Suppl Fig7a
p + scale_color_gradientn(colors = mypal)

cds_from_seurat@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(cds_from_seurat[["RNA"]])
rowData(cds_from_seurat)$gene_name <- rownames(cds_from_seurat)
rowData(cds_from_seurat)$gene_short_name <- rowData(cds_from_seurat)$gene_name
Pseudotime_cds_pr_test_Seurat <- graph_test(cds_from_seurat, neighbor_graph="principal_graph", cores=4)
pr_deg_ids_Seurat <- row.names(subset(Pseudotime_cds_pr_test_Seurat, q_value < 0.05))
gene_module_Seurat <- find_gene_modules(cds_from_seurat[pr_deg_ids_Seurat,], resolution=9e-3)
plot_cells(cds_from_seurat, genes=gene_module_Seurat, 
           show_trajectory_graph=F, 
           label_cell_groups=FALSE) 
# remove modules that do not show discrete spatial clustering
#Suppl Fig7b
plot_cells(cds_from_seurat, genes=gene_module_Seurat, 
           show_trajectory_graph=F, 
           label_cell_groups=FALSE,
           scale_to_range = F,
           norm_method = 'size_only') 

plot_cells(cds_from_seurat, genes=gene_module_Seurat %>% filter(module %in% c(2,3,4)), 
           show_trajectory_graph=F, 
           label_cell_groups=FALSE,
           scale_to_range = F,
           norm_method = 'size_only') 

gene_module_Seurat <- filter(gene_module_Seurat, gene_module_Seurat$module == 4 | gene_module_Seurat$module == 2 |gene_module_Seurat$module == 3)

#Fig. 3i
cell_group_df_Seurat <- tibble::tibble(cell=row.names(colData(cds_from_seurat)), 
                                       cell_group=colData(cds_from_seurat)$Genotype)
agg_mat_Seurat <- aggregate_gene_expression(cds_from_seurat, gene_module_Seurat, cell_group_df_Seurat)
row.names(agg_mat_Seurat) <- stringr::str_c("Module ", row.names(agg_mat_Seurat))

pheatmap::pheatmap(agg_mat_Seurat,
                   scale="row", clustering_method="ward.D2", color = mypal,
                   )
Pseudotime_cds_pr_test_Seurat$id = row.names(Pseudotime_cds_pr_test_Seurat)
x = merge(gene_module_Seurat,Pseudotime_cds_pr_test_Seurat, by = 'id')
write.csv(x, 'Cluster_MG_PseudotimeModule.csv')

#Fig. 3j
d1 = subset(gene_module_Seurat, subset = gene_module_Seurat$module == 3)
d1 = d1$id
d2 = subset(gene_module_Seurat, subset = gene_module_Seurat$module == 4)
d2 = d2$id
MGoverlapList <- read.csv('MG_diseaseStats_OverlapComp.csv', header = T)
erm = MGoverlapList$Early.Response.microglia
lrm = MGoverlapList$Late.Response.Microglia
dam = MGoverlapList$DAM
maxl = max(c(length(d1), length(d2), length(dam), length(erm), length(lrm)))
data = data.frame(DAM = c(dam, rep(NA, maxl - length(dam))),
                  Module1 = c(d1, rep(NA, maxl - length(d1))),
                  Module3 = c(d2, rep(NA, maxl - length(d2))),
                  Early.Response.microglia = c(erm, rep(NA, maxl - length(erm))),
                  Late.Response.Microglia = c(lrm, rep(NA, maxl - length(lrm))))
MGoverlapList = data
BiocManager::install("GeneOverlap")
library(GeneOverlap)

## D1 vs. DAM
test_DAMvsDiseaseModule1 <- newGeneOverlap(MGoverlapList$DAM, MGoverlapList$Module1 )
test_DAMvsDiseaseModule1 <- testGeneOverlap(test_DAMvsDiseaseModule1)
print(test_DAMvsDiseaseModule1)
## D2 vs. DAM 
test_DAMvsDiseaseModule1 <- newGeneOverlap(MGoverlapList$DAM, MGoverlapList$Module3 )
test_DAMvsDiseaseModule1 <- testGeneOverlap(test_DAMvsDiseaseModule1)
print(test_DAMvsDiseaseModule1)
## D1 vs. ERM
test_DAMvsDiseaseModule1 <- newGeneOverlap(MGoverlapList$Early.Response.microglia, MGoverlapList$Module1 )
test_DAMvsDiseaseModule1 <- testGeneOverlap(test_DAMvsDiseaseModule1)
print(test_DAMvsDiseaseModule1)
## D2 vs. ERM
test_DAMvsDiseaseModule1 <- newGeneOverlap(MGoverlapList$Early.Response.microglia, MGoverlapList$Module3 )
test_DAMvsDiseaseModule1 <- testGeneOverlap(test_DAMvsDiseaseModule1)
print(test_DAMvsDiseaseModule1)
## D1 vs. LRM
test_DAMvsDiseaseModule1 <- newGeneOverlap(MGoverlapList$Late.Response.Microglia, MGoverlapList$Module1 )
test_DAMvsDiseaseModule1 <- testGeneOverlap(test_DAMvsDiseaseModule1)
print(test_DAMvsDiseaseModule1)
## D2 vs. LRM
test_DAMvsDiseaseModule1 <- newGeneOverlap(MGoverlapList$Late.Response.Microglia, MGoverlapList$Module3 )
test_DAMvsDiseaseModule1 <- testGeneOverlap(test_DAMvsDiseaseModule1)
print(test_DAMvsDiseaseModule1)
# Enter data to Prism to plot





