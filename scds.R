library(Seurat)
library(scds)

scrna= readRDS("~/Desktop/AGIQ424_RUL_D_CD45neg_20190922.RDS")
sce= as.SingleCellExperiment(scrna)
sce = cxds(sce)
sce = bcds(sce,retRes = TRUE,verb=TRUE)
sce = cxds_bcds_hybrid(sce)
seurat=as.Seurat(sce)
scrna=seurat
saveRDS(scrna, file="~/Desktop/Maria_scds_test.RDS")

process<- function(scrna){
  sce= as.SingleCellExperiment(scrna)
  sce = cxds(sce)
  sce = bcds(sce,retRes = TRUE,verb=TRUE)
  sce = cxds_bcds_hybrid(sce)
  seurat=as.Seurat(sce)
  meta= seurat@meta.data
  return(meta)
}

setwd("~/Desktop/labproject/projects/Morrisey/Maria/lungMAP/")
scrna=readRDS(file = 'R710000511_COPD_D_CD45neg_20181109/Seurat/R710000511_COPD_D_CD45neg_20181109.RDS')
sce= as.SingleCellExperiment(scrna)
sce = cxds(sce)
sce = bcds(sce,retRes = TRUE,verb=TRUE)
sce = cxds_bcds_hybrid(sce)
seurat=as.Seurat(sce)
R710000511_meta= seurat@meta.data

scrna=readRDS(file = 'R710000546_COPD_D_CD45neg_20190517/Seurat/R710000546_COPD_D_CD45neg_20190517.RDS')
R710000546_meta= process(scrna)


