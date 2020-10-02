library(Seurat)
library(dplyr,quietly = T)
library(cowplot,quietly = T)
library(scExtras)
library(readxl)
library(ggplot2)

cpallette=c("#64B2CE", "#DA5724", "#74D944", "#CE50CA", "#C0717C", "#CBD588", "#5F7FC7", 
            "#673770", "#D3D93E", "#8569D5", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
            "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
            "#8A7C64", "#599861")


#################################################################
#################################################################
######Subset mesenchyme and perform trajectory analysis#########
#################################################################
#################################################################
scrna= readRDS(file = 'MBJK_WT_I73T_Intergrated.RDS')

#Subset data
mes= subset(scrna,idents = c(0,1,3,6,8,10,13))

DefaultAssay(mes) <- "integrated"

scrna <- RunPCA(mes, npcs = 50)
ElbowPlot(scrna, ndims = 50)

dims <- 1:20

#Recluster
mes= ClusterDR(mes,dims,n.neighbors=20,findallmarkers=F,resolution=c(0, 0.2,0.4,0.6,0.8,1.0),min.dist=0.1)
mes@misc[["findallmarkers"]] <- FindAllMarkers(object = mes, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

mes$var_cluster <- mes@active.ident

DefaultAssay(mes) <- "SCT"

#Trajectory
mes <- runSlingshot(mes,reduction='umap',approx_points = 200)
png(file="~/Desktop/slingshot_I72T.png", height = 10, width = 10, units = "in", res=500)
CurvePlot(mes,sds = mes@misc$sds$data,group.by = "var_cluster", cols = cpallette)
dev.off()

saveRDS(mes, file="~/Desktop/MBJK_WT_I73T_Integrated_mes.RDS")

meta= mes@meta.data
anno = read.csv("~/Desktop/Jeremy/subset/I73T_anno.csv")
meta$var_cluster= as.numeric(as.character(meta$var_cluster))
meta$id= rownames(meta)
meta = left_join(meta, anno, by= c("var_cluster"="cluster"))
rownames(meta)= meta$id
meta = meta %>% select(-id)
mes@meta.data = meta
Idents(mes) <- 'celltype'

mes <- runSlingshot(mes,reduction='umap',group.by = "celltype",approx_points = 200,start.clus=NULL,end.clus='Myofibroblast',extend= "n",stretch=0)
png(file="~/Desktop/slingshot_I73T_endmyo.png", height = 10, width = 10, units = "in", res=500)
CurvePlot(mes,sds = mes@misc$sds$data,group.by = "celltype", cols = cpallette)
dev.off()

mes <- runSlingshot(mes,reduction='umap',group.by = "celltype",approx_points = 200,start.clus='MANC',end.clus='Myofibroblast',extend= "n",stretch=0)
png(file="~/Desktop/slingshot_I73T_manc_myo.png", height = 10, width = 10, units = "in", res=500)
CurvePlot(mes,sds = mes@misc$sds$data,group.by = "celltype", cols = cpallette)
dev.off()

mes <- runSlingshot(mes,reduction='umap',group.by = "celltype",approx_points = 200,start.clus='AMP',end.clus='Myofibroblast',extend= "n",stretch=0)
png(file="~/Desktop/slingshot_I73T_amp_myo.png", height = 10, width = 10, units = "in", res=500)
CurvePlot(mes,sds = mes@misc$sds$data,group.by = "celltype", cols = cpallette)
dev.off()

#Subset clusters of interes, run slingshot and plot it back in full umap

mes_sub = subset(mes,idents = c("AMP","MANC","MANC like","Wnt2","Myofibroblast","Proliferating","VSM","Cluster10_from_I73T_dataset"))

mes_sub <- runSlingshot(mes_sub,reduction='umap',group.by = "celltype",approx_points = 200,start.clus='AMP',end.clus='Myofibroblast',extend= "n",stretch=0)
png(file="~/Desktop/slingshot_I73T_amp_myo_sub.png", height = 10, width = 10, units = "in", res=500)
CurvePlot(mes,sds = mes_sub@misc$sds$data,group.by = "celltype", cols = cpallette)
dev.off()

mes_sub <- runSlingshot(mes_sub,reduction='umap',group.by = "celltype",approx_points = 200,start.clus='MANC',end.clus='Myofibroblast',extend= "n",stretch=0)
png(file="~/Desktop/slingshot_I73T_manc_myo_sub.png", height = 10, width = 10, units = "in", res=500)
CurvePlot(mes,sds = mes_sub@misc$sds$data,group.by = "celltype", cols = cpallette)
dev.off()

mes_sub <- runSlingshot(mes_sub,reduction='umap',group.by = "celltype",approx_points = 200,start.clus='Wnt2',end.clus='Myofibroblast',extend= "n",stretch=0)
png(file="~/Desktop/slingshot_I73T_wnt2_myo_sub.png", height = 10, width = 10, units = "in", res=500)
CurvePlot(mes,sds = mes_sub@misc$sds$data,group.by = "celltype", cols = cpallette)
dev.off()


mes_sub <- runSlingshot(mes_sub,reduction='umap',group.by = "celltype",approx_points = 200,end.clus='Proliferating',start.clus='Myofibroblast',extend= "n",stretch=0)
png(file="~/Desktop/slingshot_I73T_myo_prof_sub.png", height = 10, width = 10, units = "in", res=500)
CurvePlot(mes,sds = mes_sub@misc$sds$data,group.by = "celltype", cols = cpallette)
dev.off()


