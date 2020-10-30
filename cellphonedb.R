library(data.table)
library(igraph)
library(scExtras)
library(RColorBrewer)
library(NMF)
######################################################################
autocurve.edges2 <-function (graph, start = 0.5)
{
  cm <- count.multiple(graph)
  mut <-is.mutual(graph)  #are connections mutual?
  el <- apply(get.edgelist(graph, names = FALSE), 1, paste,
              collapse = ":")
  ord <- order(el)
  res <- numeric(length(ord))
  p <- 1
  while (p <= length(res)) {
    m <- cm[ord[p]]
    mut.obs <-mut[ord[p]] #are the connections mutual for this point?
    idx <- p:(p + m - 1)
    if (m == 1 & mut.obs==FALSE) { #no mutual conn = no curve
      r <- 0
    }
    else {
      r <- seq(-start, start, length = m)
    }
    res[ord[idx]] <- r
    p <- p + m
  }
  res
}
######################################################################
setwd("~/Desktop/linux/home/bapoorva/Desktop/server/Morrisey/Maria/lungMAP")
scrna= readRDS("Periph_5sample_Intergrate/epi_subset/SeuratAnn.RDS")
scrna2 = readRDS("COPD_4sample_Intergrate/epi_subset/SeuratAnn.RDS")

peri_counts = as.data.frame(as.matrix(scrna@assays$integrated@data)) %>% rownames_to_column() %>% rename('Gene'='rowname')
copd_counts = as.data.frame(as.matrix(scrna2@assays$integrated@data)) %>% rownames_to_column() %>% rename('Gene'='rowname')
write.table(peri_counts, file="../cellphoneDB/peri_counts2.txt",row.names = F,col.names=T,quote = F,sep="\t")
write.table(copd_counts, file="../cellphoneDB/copd_counts2.txt",row.names = F,col.names=T,quote = F,sep="\t")

peri_meta= scrna@meta.data %>% rownames_to_column() %>%
          select(rowname,var_celltype) %>%
          rename('Cell'='rowname','cell_type'= 'var_celltype')

copd_meta= scrna2@meta.data %>% rownames_to_column() %>%
  select(rowname,var_celltype) %>%
  rename('Cell'='rowname','cell_type'= 'var_celltype')

write.table(peri_meta, file="../cellphoneDB/peri_meta.txt", row.names = F,col.names=T,quote = F,sep="\t")
write.table(copd_meta, file="../cellphoneDB/copd_meta.txt", row.names = F,col.names=T,quote = F,sep="\t")

df_pv= fread("cellphoneDB/copd_out_stats/pvalues.txt")
df_pv2 = df_pv %>% gather("interacting_clust","pval",`AT1_AT2|AT1_AT2`:`Sec/Gob|Sec/Gob`)
df_pv2 = df_pv2[df_pv2$pval < 0.05,]

peri_pv= fread("cellphoneDB/peri_out_stats/pvalues.txt")
peri_pv2 = peri_pv %>% gather("interacting_clust","pval",`AT1_AT2|AT1_AT2`:`Sec/Gob|Sec/Gob`)
peri_pv2 = peri_pv2[peri_pv2$pval < 0.05,]


copd_res= df_pv2[df_pv2$receptor_a != df_pv2$receptor_b,]
copd_res = copd_res %>% dplyr::select(interacting_pair,gene_a,gene_b,receptor_a,receptor_b,interacting_clust) %>% 
            separate(interacting_clust,c("A","B"),sep="[|]")
copd_res$Ligand_cluster=ifelse(copd_res$receptor_a == "FALSE",copd_res$A,copd_res$B)
copd_res$Receptor_cluster=ifelse(copd_res$receptor_a == "TRUE",copd_res$A,copd_res$B)
copd_res$Ligand=ifelse(copd_res$receptor_a == "FALSE",copd_res$gene_a,copd_res$gene_b)
copd_res$Receptor=ifelse(copd_res$receptor_a == "TRUE",copd_res$gene_a,copd_res$gene_b)

copd_res = copd_res %>% dplyr::select(interacting_pair,Ligand_cluster:Receptor)

write.csv(result,file="cellphoneDB/copd_LRpairs.csv",row.names = F)

png("~/Desktop/COPD_LR_network.png", width = 8, height = 8,res=400, units = "in")
plotLRnetwork(copd_res)
dev.off()
png("~/Desktop/COPD_LR_Heatmap.png", width = 8, height = 8,res=400, units = "in")
plotLRHeatmap(copd_res, clusterby = "both")
dev.off()

result= peri_pv2[peri_pv2$receptor_a != peri_pv2$receptor_b,]
result = result %>% dplyr::select(interacting_pair,gene_a,gene_b,receptor_a,receptor_b,interacting_clust)
result = result %>% separate(interacting_clust,c("A","B"),sep="[|]")
result$Ligand_cluster=ifelse(result$receptor_a == "FALSE",result$A,result$B)
result$Receptor_cluster=ifelse(result$receptor_a == "TRUE",result$A,result$B)
result$Ligand=ifelse(result$receptor_a == "FALSE",result$gene_a,result$gene_b)
result$Receptor=ifelse(result$receptor_a == "TRUE",result$gene_a,result$gene_b)
result = result %>% dplyr::select(interacting_pair,Ligand_cluster:Receptor)

result = result %>% dplyr::select(-A,-B) %>% arrange(receptor_a)
write.csv(result,file="cellphoneDB/periph_LRpairs.csv",row.names = F)

png("~/Desktop/Peri_LR_network.png", width = 8, height = 8,res=400, units = "in")
plotLRnetwork(result, filtermin = 8,filtermax = 20)
dev.off()

png("~/Desktop/Peri_LR_Heatmap.png", width = 8, height = 8,res=400, units = "in")
plotLRHeatmap(result, clusterby = "both")
dev.off()

