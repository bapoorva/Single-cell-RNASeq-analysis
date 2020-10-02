library(tradeSeq)
library(dplyr)
library(tidyr)
library(slingshot)
library(scExtras)
library(tibble)
library(Seurat)
library(ComplexHeatmap)

cpallette=c("#64B2CE", "#DA5724" ,"#74D944" ,"#CE50CA", "#C0717C" ,"#CBD588", "#5F7FC7", "#673770" ,"#D3D93E",
"#8569D5", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD" ,"#D14285", "#6DDE88", "#652926",
"#7FDCC0","#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861")

scrna= readRDS("~/Desktop/Jeremy/subset/MBJK_WT_C2C_Integrated_mes.RDS")
counts <- scrna@assays$RNA@counts %>% as.matrix(.)
counts_var <- scrna@assays$RNA@counts[scrna@assays$SCT@var.features, ] %>% as.matrix(.)


meta= scrna@meta.data
anno = read.csv("~/Desktop/C2C_anno.csv")
meta$var_cluster= as.numeric(as.character(meta$var_cluster))
meta$id= rownames(meta)
meta = left_join(meta, anno, by= c("var_cluster"="cluster"))
rownames(meta)= meta$id
meta = meta %>% select(-id)
scrna@meta.data = meta
Idents(scrna) <- 'celltype'

mes_sub = subset(scrna,idents = c("AMP","MANC","Wnt2","Myofibroblast","Proliferative"))
mes_sub <- runSlingshot(mes_sub,reduction='umap',group.by = "celltype",approx_points = 200,start.clus='MANC',end.clus='Myofibroblast',extend= "n",stretch=0)
png(file="~/Desktop/slingshot_C2CT_manc_myo_sub.png", height = 10, width = 10, units = "in", res=500)
CurvePlot(scrna,sds = mes_sub@misc$sds$data,group.by = "celltype", cols = cpallette)
dev.off()


crv= mes_sub@misc$sds$data

#Evaluate the optimal number of knots required for fitGAM
set.seed(5)
icMat <- evaluateK(counts = counts_var, sds = crv, k = 3:20,
                   nGenes = 200, verbose = T, plot=T)

set.seed(7)
pseudotime <- slingPseudotime(crv, na=FALSE)
cellWeights <- slingCurveWeights(crv)

sc= list(counts=counts,pseudotime=pseudotime,cellWeights=cellWeights)
saveRDS(sc, "~/Desktop/linux/home/bapoorva/Desktop/server/C2C_sc.RDS")

sc= list(counts=counts_var,pseudotime=pseudotime,cellWeights=cellWeights)
saveRDS(sc, "~/Desktop/linux/home/bapoorva/Desktop/server/C2C_sc_var.RDS")

##############RUN THIS SECTION ALONE ON CLUSTER##################

BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- 20
sc= readRDS("C2C_sc_var.RDS")
counts= sc$counts
pseudotime= sc$pseudotime
cellWeights= sc$cellWeights
sce <- fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,nknots = 8, verbose = FALSE,parallel=TRUE, BPPARAM = BPPARAM)
saveRDS(sce, file="~/sce_c2c_var.RDS")

################################################################################################

list= readRDS("~/Desktop/linux/home/bapoorva/Desktop/server/C2C_sc_var.RDS")
counts=list$counts
pseudotime=list$pseudotime
cellWeights=list$cellWeights

#Read in results of fitgam after running it on cluster
sce= readRDS("~/Desktop/linux/home/bapoorva/Desktop/server/sce_c2c_var.RDS")
assoRes <- associationTest(sce)
head(assoRes)


plotGeneCount(curve = crv, counts = counts,
              clusters = apply(slingClusterLabels(crv), 1, which.max),
              models = sce)


startRes <- startVsEndTest(sce, lineages = T) %>% rownames_to_column('gene')
oStart <- order(startRes$waldStat_lineage2, decreasing = TRUE)
sigGeneStart <- names(sce)[oStart[3]]

plotSmoothers(sce, counts, gene = sigGeneStart)
plotGeneCount(crv, counts, gene = sigGeneStart)


# earlyDERes <- earlyDETest(sce, knots = c(1, 5))
# oEarly <- order(earlyDERes$waldStat, decreasing = TRUE)
# head(rownames(earlyDERes)[oEarly])
# 
# plotSmoothers(sce[[sigGeneStart]], counts, gene = rownames(earlyDERes)[oEarly][2])

#1st run heatmap to select pseudotime values
cellcolorpal<- cpallette[1:length(levels(mes_sub$celltype))]
names(cellcolorpal) <- levels(mes_sub$celltype)

  l =names(slingLineages(crv))[2]
    l <- tolower(l)
    pvalue <- sym(paste0('pvalue_',l))
    #logFC <- sym(paste0('logFC',l))
    res= startRes
    topgenes <- res %>% 
      top_n(n = -10, wt =!!pvalue) %>%
      #top_n(n = 50, wt =abs(!!logFC)) %>% 
      #arrange(!!logFC) %>%
      pull(gene)
    #map(topgenes,genePlot,plotdir = currplotdir, cols=cellcolorpal)
    png(filename = '~/Desktop/TopGenes_Heatmap.png',width = 12,height=8,units = 'in',res=300)
    ComplexHeatmap::draw(plotLineageHeatMap2(mes_sub,sdsname = 'sds', features=topgenes,lineage=l,col=cellcolorpal))
    dev.off()
    

#Proliferative to Myo - no fc sort
testRes <- startVsEndTest(sce,lineages=T,pseudotimeValues=c(5,15)) %>% rownames_to_column('gene')
res= testRes
topgenes2 <- res %>% 
  top_n(n = -10, wt =!!pvalue) %>%
  #top_n(n = 50, wt =abs(!!logFC)) %>% 
  #arrange(!!logFC) %>%
  pull(gene)
topgenes2 = topgenes2[1:50]
#map(topgenes,genePlot,plotdir = currplotdir, cols=cellcolorpal)
png(filename = '~/Desktop/TopGenes_Heatmap_new.png',width = 12,height=8,units = 'in',res=300)
ComplexHeatmap::draw(plotLineageHeatMap2(mes_sub,sdsname = 'sds', features=topgenes2,lineage=l,col=cellcolorpal))
dev.off()
#####3
#Proliferative to Myo - fc sorted
pvalue <- sym(paste0('pvalue_',l))
logFC <- sym(paste0('logFC',l))

testRes <- startVsEndTest(sce,lineages=T,pseudotimeValues=c(5,15)) %>% rownames_to_column('gene')
res= testRes
topgenes2 <- res %>% 
  top_n(n = -10, wt =!!pvalue) %>%
  top_n(n = 50, wt =abs(!!logFC)) %>% 
  arrange(!!logFC) %>%
  pull(gene)
topgenes2 = topgenes2[1:50]
#map(topgenes,genePlot,plotdir = currplotdir, cols=cellcolorpal)
png(filename = '~/Desktop/TopGenes_Heatmap_new3.png',width = 12,height=8,units = 'in',res=300)
ComplexHeatmap::draw(plotLineageHeatMap2(mes_sub,sdsname = 'sds', features=topgenes2,lineage=l,col=cellcolorpal))
dev.off()

#MANC to Proliferative- fc sorted
testRes <- startVsEndTest(sce,lineages=T,pseudotimeValues=c(0,10)) %>% rownames_to_column('gene')
res= testRes
topgenes2 <- res %>% 
  top_n(n = -10, wt =!!pvalue) %>%
  top_n(n = 50, wt =abs(!!logFC)) %>% 
  arrange(!!logFC) %>%
  pull(gene)
topgenes2 = topgenes2[1:50]
#map(topgenes,genePlot,plotdir = currplotdir, cols=cellcolorpal)
png(filename = '~/Desktop/TopGenes_Heatmap_new2.png',width = 12,height=8,units = 'in',res=300)
ComplexHeatmap::draw(plotLineageHeatMap2(mes_sub,sdsname = 'sds', features=topgenes2,lineage=l,col=cellcolorpal))
dev.off()

write.csv(testRes, file= "testRes_manc_to_prol.csv", row.names = F)

#Proliferative to Myo - fc sorted - top 100
testRes2 <- startVsEndTest(sce,lineages=T,pseudotimeValues=c(5,15)) %>% rownames_to_column('gene')
res= testRes2
topgenes <- res %>% 
  top_n(n = -10, wt =!!pvalue) %>%
  top_n(n = 100, wt =abs(!!logFC)) %>% 
  arrange(!!logFC) %>%
  pull(gene)
#topgenes2 = topgenes2[1:50]
#map(topgenes,genePlot,plotdir = currplotdir, cols=cellcolorpal)
png(filename = '~/Desktop/TopGenes_Heatmap_new4.png',width = 12,height=8,units = 'in',res=300)
ComplexHeatmap::draw(plotLineageHeatMap2(mes_sub,sdsname = 'sds', features=topgenes,lineage=l,col=cellcolorpal))
dev.off()

################################################################################
################################################################################
################################################################################
################################################################################

plotLineageHeatMap2 <- function(object,sdsname,features,lineage='lineage1',col, group.by="celltype"){
  if (is.null(x = group.by)) {
    stop("Please Enter the variable that was used to define groups in Slingshot, ie var_celltype, var_cluster etc.")
  }
  
  ## Maybe add curve is user puts in integer
  
  sds <- object@misc[[sdsname]]$data
  clusterinlineage <- slingLineages(sds)[[stringr::str_to_title(lineage)]]
  ### Should add a check for lineages in model
  group.by <- sym(group.by)
  qlineage <- quo(lineage)
  
  cells <- inner_join(
    slingCurveWeights(sds,as.probs=T) %>% as.data.frame() %>%
      setNames(gsub('curve','lineage',names(.))) %>%
      rownames_to_column('cellid') %>%
      gather(curve,w,-cellid) %>%
      group_by(cellid) %>%
      top_n(n=1,wt=w) %>%
      filter(curve==!!lineage),
    slingPseudotime(sds) %>% as.data.frame %>%
      setNames(gsub('curve','lineage',names(.))) %>%
      rownames_to_column('cellid') %>%
      select(cellid,!!qlineage) %>%
      dplyr::rename('time'=lineage)
  ) %>% arrange(time) %>%
    inner_join(., object@meta.data %>% rownames_to_column('cellid')) %>%
    filter(!!group.by %in% clusterinlineage) %>%
    mutate(group_by=factor(group.by,levels=clusterinlineage))
  
  data <- FetchData(object=object, vars=features,cells=cells$cellid) %>% t(.)
  mat_scaled = t(scale(t(data)))
  
  f1=circlize::colorRamp2(c(-2,0,2), c('skyblue1', "grey10","yellow"))
  
  
  col_fun =circlize::colorRamp2(c(0, 20), c("blue", "red"))
  
  ha = HeatmapAnnotation(
    #pseudotime=anno_lines(cells$time),
    pseudotime=anno_barplot(cells$time, gp = gpar(col = "#296EFA"),border = F,bar_width = 1),
    celltype = cells %>% pull(group.by),
    col = list(celltype = col
    ))
  
  ht <- Heatmap(mat_scaled,
                col=f1,
                show_row_dend = F,
                row_names_side='left',
                show_column_names = F,
                cluster_columns = F,
                cluster_rows = F,
                top_annotation = ha
  )
  return(ht)
  
  
  
}
