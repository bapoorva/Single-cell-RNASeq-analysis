library(DoubletFinder)
library(dplyr)
library(tidyr)

ctrl= readRDS("~/Desktop/vera/VK_Healthy.RDS")
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list <- paramSweep_v3(ctrl, PCs = 1:15, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

pK=as.numeric(as.character(bcmvn$pK))
BCmetric=bcmvn$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose,lwd=2,col='red',lty=2)
title("The BCmvn distributions")
text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations=ctrl@meta.data$var_cluster
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*length(rownames(ctrl)))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
ctrl <- doubletFinder_v3(ctrl, PCs = 1:15, pN = 0.25, pK = pK_choose, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
ctrl2 <- doubletFinder_v3(ctrl, PCs = 1:15, pN = 0.25, pK = pK_choose, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)

m1=ctrl@meta.data
m2= ctrl2@meta.data
m1 = m1 %>% dplyr::select(var_cluster,DF.classifications_0.25_0.19_1281)
m2 = m2 %>% dplyr::select(var_cluster,DF.classifications_0.25_0.19_1124)

m1$cells = rownames(m1)
m2$cells = rownames(m2)
# df= inner_join(m1,m2, by= "cells")
# df2= df[df$DF.classifications_0.25_0.19_1281 == "Doublet",]
# 
# df2 = df2 %>% select(-var_cluster.y)
# a1= df2 %>% select(var_cluster.x:cells)
# a2= df2 %>% select(var_cluster.x,cells:DF.classifications_0.25_0.19_1124)
 a1 = m1[m1$DF.classifications_0.25_0.19_1281 =="Doublet",]
 a2= m2[m2$DF.classifications_0.25_0.19_1124 == "Doublet",]
 a1= a1 %>% group_by(var_cluster) %>% summarise(n())
 colnames(a1)=c("var_cluster","DF.classifications_0.25_0.19_1281")
 a2= a2 %>% group_by(var_cluster) %>% summarise(n())
 colnames(a2)=c("var_cluster","DF.classifications_0.25_0.19_1124")
 
cnt= inner_join(a1,a2,by="var_cluster")
c.cnt=as.data.frame(table(Idents(object=ctrl)))

cnt= left_join(c.cnt,cnt, by=c("Var1"="var_cluster"))
cnt$DF.classifications_0.25_0.19_1281[is.na(cnt$DF.classifications_0.25_0.19_1281)]<- 0
cnt$DF.classifications_0.25_0.19_1124[is.na(cnt$DF.classifications_0.25_0.19_1124)]<- 0
cnt$perc= (cnt$DF.classifications_0.25_0.19_1281/cnt$Freq)*100
cnt = cnt %>% select(-DF.classifications_0.25_0.19_1124)
colnames(cnt)= c("Cluster","Number of Cells in cluster","Number of doublets","Percentage of cells in cluster identified as doublets")

write.csv(cnt, file="~/Desktop/Healthy_Doublet.csv", row.names = F)
