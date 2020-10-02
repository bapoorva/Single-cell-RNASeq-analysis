library(monocle)
library(garnett)
library(org.Mm.eg.db)

#Load data
scrna=readRDS("~/Desktop/Server/Morrisey/Jarod/scRNA/JZ_lung_timeseries/P3_SCT/P3.RDS")

#Convert to monocle
#cds= importCDS(scrna, import_all=FALSE)
gene_annotation <- as.data.frame(rownames(scrna@reductions[["pca"]]@feature.loadings), row.names = rownames(scrna@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"
featureData <- new("AnnotatedDataFrame", data=gene_annotation)

cell_metadata <- as.data.frame(scrna@assays[["SCT"]]@counts@Dimnames[[2]], row.names = scrna@assays[["SCT"]]@counts@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"
phenoData <- new("AnnotatedDataFrame", data=cell_metadata)

New_matrix <- scrna@assays[["SCT"]]@counts
New_matrix <- New_matrix[rownames(scrna@reductions[["pca"]]@feature.loadings), ]
expression_matrix <- New_matrix


### Construct the basic cds object

cds= newCellDataSet(expression_matrix,phenoData  = phenoData,featureData  = featureData)

cds <- estimateSizeFactors(cds)

#test marker files
marker_file_path <- "/Users/bapoorva/Desktop/Server/Morrisey/JohnLeach/scRNA/markers_garnett.txt"
marker_check <- check_markers(cds, marker_file_path,
                              db=org.Mm.eg.db,
                              cds_gene_id_type = "SYMBOL",
                              marker_file_gene_id_type = "SYMBOL")

pdf("../../Jarod/scRNA/JZ_lung_timeseries/P3_SCT/Marker_check.pdf", height=15,width=8)
plot_markers(marker_check)
dev.off()


#train classifier
set.seed(260)

classifier <- train_cell_classifier(cds = cds,
                                         marker_file = marker_file_path,
                                         db=org.Mm.eg.db,
                                         cds_gene_id_type = "SYMBOL",
                                         num_unknown = 50,
                                         marker_file_gene_id_type = "SYMBOL")


feature_genes <- get_feature_genes(classifier,
                                   node = "root",
                                   db = org.Mm.eg.db)
head(feature_genes)

#Load test data
new_sc=readRDS("~/Desktop/Server/Morrisey/JohnLeach/scRNA/CD45_depleted_adult/Seurat/CD45_depleted_adult.RDS")

#Make cds object from test data
gene_annotation2 <- as.data.frame(rownames(new_sc@reductions[["pca"]]@feature.loadings), row.names = rownames(new_sc@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation2) <- "gene_short_name"
featureData2 <- new("AnnotatedDataFrame", data=gene_annotation2)

cell_metadata2 <- as.data.frame(new_sc@assays[["SCT"]]@counts@Dimnames[[2]], row.names = new_sc@assays[["SCT"]]@counts@Dimnames[[2]])
colnames(cell_metadata2) <- "barcode"
phenoData2 <- new("AnnotatedDataFrame", data=cell_metadata2)

New_matrix2 <- new_sc@assays[["SCT"]]@counts
New_matrix2 <- New_matrix2[rownames(new_sc@reductions[["pca"]]@feature.loadings), ]
expression_matrix2 <- New_matrix2

new_cds= newCellDataSet(expression_matrix2,phenoData  = phenoData2,featureData  = featureData2)
new_cds <- estimateSizeFactors(new_cds)
#Apply classifier

new_cds <- classify_cells(new_cds, classifier,
                           db = org.Mm.eg.db,
                           cluster_extend = TRUE,
                           cds_gene_id_type = "SYMBOL")







