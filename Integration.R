library(Seurat)
library(dplyr,quietly = T)
library(cowplot,quietly = T)
library(scExtras)
library(readxl)
library(ggplot2)
options(future.globals.maxSize = 4000 * 1024^2)

cpallette=c("#64B2CE", "#DA5724", "#74D944", "#CE50CA", "#C0717C", "#CBD588", "#5F7FC7", 
            "#673770", "#D3D93E", "#8569D5", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
            "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
            "#8A7C64", "#599861")

outdir<-'plots' # All results will be saved here plots, RData file 
name1<-'WT' # specify project name,this will also be the Rdata file name
name2<-'C2C' # specify project name,this will also be the Rdata file name

org= "mouse"
input1 <- c('../Adult_WT_CD45neg/filtered')
input2 <- ("../Adult_C2C_CD45neg/filtered")

#Read input datasets
WT= RunQC(dir=outdir,org=org,name=name1,files=input1 ,filter=T, doubletdetection = T,UpperMitoCutoff=10)
C2C= RunQC(dir=outdir,org=org,name=name2,files=input2 ,filter=T, doubletdetection = T,UpperMitoCutoff=10)

#create list of data
obj=list()
obj[["WT"]] =WT
obj[["C2C"]] = C2C

#Sctransform data
for (i in 1:length(obj)) {
  obj[[i]] <- SCTransform(obj[[i]], verbose = FALSE)
}

#select features for downstream integration
features <- SelectIntegrationFeatures(object.list = obj, nfeatures = 3000)
obj <- PrepSCTIntegration(object.list = obj, anchor.features = features,verbose = FALSE)

#identify anchors and integrate the datasets
anchors <- FindIntegrationAnchors(object.list = obj, normalization.method = "SCT", anchor.features = features, verbose = FALSE)
scrna <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = FALSE)

#Cluster and visualize
npcs=50
scrna= PCATools(scrna, npcs=npcs, jackstraw=T, plotdir = outdir)
png(file=paste0(outdir,"/elbowplot_C2C.png",sep=""), height = 15, width = 15, units = "in", res=500)
ElbowPlot(scrna, ndims = npcs)
dev.off()

png(file=paste0(outdir,"/jackstraw_C2C.png",sep=""), height = 15, width = 15, units = "in", res=500)
JackStrawPlot(scrna, dims = 1:npcs)
dev.off()


dims <- 1:25

scrna <- ClusterDR(scrna,dims=dims,n.neighbors =30,findallmarkers = F)
scrna@meta.data$orig.ident=factor(scrna@meta.data$orig.ident,levels=c("C2C","WT"))
markergenes =c("Sox17",	"Vwf","Cxcl12","Fabp4",	"Car4","Prox1","Hopx",	"Ager","Sftpc","Sftpb",
               "Scgb3a2","Foxj1",	"Wnt2","Mfap5","Hhip","Notch3","Pdgfrb","Stc1","Alas2",	"Wt1"	,"Msln",
               "Ptprc",	"Lgals3","Mki67","Pecam1","Pdgfra","Acta2","Ascl2","Scgb1a1")


DefaultAssay(scrna) <- "SCT"

for(g in markergenes){
  try({
    cowplot::plot_grid(
      DimPlot(object = scrna, reduction = "umap",cols = cpallette,label=T,pt.size = 1) + NoLegend(),
      DimPlot(scrna,group.by ='orig.ident')+ theme(legend.position = "bottom"),
      FeaturePlot(scrna,g,order = T,),
      VlnPlot(scrna,g,pt.size = 0,cols = cpallette,split.by = 'orig.ident') + theme(legend.position = "bottom"),
      ncol=2
    ) %>%  cowplot::save_plot(filename = paste0(outdir,'/',g,'_UMAP_Vln.png'),base_height = 12,base_width = 16)
  })
}

#scrna <- RunLigRec(scrna,org=org, group.by = "var_cluster")
saveRDS(scrna,file = 'IntergratedSeurat_C2C.RDS')

DefaultAssay(scrna) <- "integrated"

scrna_old=scrna
prof=c("Top2a","Ccnb2","Hmgb2","Mki67")
other = c("Epcam", "Sftpc", "Hopx", "Pecam", "Acta2", "Spp1", "Dcn", "Postn", "Cnn")
types=c("prof","other")
for(dim in seq(7,27,2)){
  scrna <- ClusterDR(scrna,dims=1:dim,n.neighbors =30,findallmarkers = F)
  for (i in types){
    j=eval(parse(text=i))
    cowplot::plot_grid(
      DimPlot(object = scrna,reduction = "umap",cols = cpallette,label=T,pt.size = 1) + NoLegend(),
      FeaturePlot(scrna,j,order = T,),
      ncol=2)%>%  cowplot::save_plot(filename = paste0('plots/dim_',dim,'_',i,'_c2c_cd45neg.png'),base_height = 10,base_width = 20)
  }
}

scrna_old = scrna

scrna= DietSeurat(scrna_old,assays=c("integrated","SCT"),dimreducs = names(scrna_old@reductions),scale.data = T)
saveRDS(scrna,file = paste0(outdir,'/MBJK_WT_C2C_Intergrated_reduced.RDS'))

