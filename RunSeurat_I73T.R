library(Seurat)
library(dplyr,quietly = T)
library(cowplot,quietly = T)
library(scExtras)
library(readxl)
library(ggplot2)


cpallette =cpallette=c("#64B2CE", "#DA5724", "#74D944", "#CE50CA", "#C0717C", "#CBD588", "#5F7FC7", 
                       "#673770", "#D3D93E", "#8569D5", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
                       "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
                       "#8A7C64", "#599861")

####################### Global Vars ############################
outdir<-'Seurat' # All results will be saved here plots, RData file 
projectname<-'Adult_I73T_CD45neg' # specify project name,this will also be the Rdata file name
input10x <- c('filtered') # dir(s) of the 10x output files, genes.tsv,barcodes.tsv

mouseorthologfile <- '~/Desktop/Apoorva/NGSshare/homologs/mouse_human.csv'
org<-'mouse' 
npcs<-40 #How many inital PC dimensions to compute. 
k=30 #This for nearest neighbors, 30 is default

##########################################################################
#
# Create output dirs
#
##########################################################################

dir.create(outdir,recursive = T,showWarnings = F)
plotdir <- paste0(outdir,'/plots')
dir.create(plotdir,showWarnings = F)
qcdir <-paste0(outdir,'/qc')
dir.create(qcdir,showWarnings = F)


##########################################################################
#
# Preprocess the data and create a seurat object
#
##########################################################################
scrna= RunQC(dir=outdir,org=org,name=projectname,files=input10x ,filter=T, doubletdetection = T,UpperMitoCutoff=10)
scrna = processExper(scrna ,ccscale = T, sc.transform = T)

# scrna[["percent.immune"]] <- PercentageFeatureSet(scrna,assay='RNA',features=immunegenes[immunegenes %in% rownames(scrna@assays$RNA@meta.features)])
# write.csv(scrna@meta.data, file=paste(qcdir,"/FilteredMetaDataImmune.csv",sep=""))

scrna= PCATools(scrna, npcs=npcs, jackstraw=T, plotdir = qcdir)

png(file=paste0(qcdir,"/elbowplot.png",sep=""), height = 15, width = 15, units = "in", res=500)
ElbowPlot(scrna, ndims = npcs)
dev.off()

png(file=paste0(qcdir,"/jackstraw.png",sep=""), height = 15, width = 15, units = "in", res=500)
JackStrawPlot(scrna, dims = 1:npcs)
dev.off()

npcs= 16
scrna <- ClusterDR(scrna,dims=1:npcs,n.neighbors =30,min.dist=0.3, findallmarkers = T)
cowplot::plot_grid(
  DimPlot(scrna),
  FeaturePlot(scrna,c("Hopx","Epcam","Sftpc","Foxj1","Pdgfrb","Pecam1","Car4","Prox1","Wt1")),ncol=2)

png(file=paste0(plotdir,"/UMAP.png",sep=""), height = 15, width = 15, units = "in", res=500)
DimPlot(scrna, cols = cpallette, pt.size = 1, label = T)
dev.off()

scrna <- RunLigRec(scrna,org=org)

saveRDS(scrna,paste0(outdir,'/',projectname,'.RDS'))



