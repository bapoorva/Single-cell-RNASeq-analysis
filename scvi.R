Sys.setenv(PATH= paste("/Users/bapoorva/opt/miniconda3/envs/scvi-env/bin/",Sys.getenv()["PATH"],sep=";"))
Sys.setenv(RETICULATE_PYTHON = "/Users/bapoorva/opt/miniconda3/envs/scvi-env/bin/python3.7")
library(reticulate)
library(Seurat)
library(SeuratData)
library(cowplot)

use_condaenv("scvi-env")
py_config()


sc <- import('scanpy', convert = FALSE)
scvi <- import('scvi', convert = FALSE)
scvi$settings$progress_bar_style = 'tqdm'

outdir<-'plots' # All results will be saved here plots, RData file 
name1<-'WT' # specify project name,this will also be the Rdata file name
name2<-'I73T' # specify project name,this will also be the Rdata file name

org= "mouse"
input1 <- c('Adult_WT_CD45neg/filtered')
input2 <- ("Adult_I73T_CD45neg/filtered/")

#Read input datasets
WT= RunQC(dir=outdir,org=org,name=name1,files=input1 ,filter=T, doubletdetection = F,UpperMitoCutoff=10)
WT$var_sample <- 'WT'
I73T= RunQC(dir=outdir,org=org,name=name2,files=input2 ,filter=T, doubletdetection = F,UpperMitoCutoff=10)
I73T$var_sample <- 'I73T'

scrna <- merge(x = WT, y = I73T)

# use seurat for variable gene selection
scrna <- NormalizeData(scrna, normalization.method = "LogNormalize", scale.factor = 10000)
#scrna[["percent.mt"]] <- PercentageFeatureSet(scrna, pattern = "^MT-")
#scrna <- subset(scrna, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
scrna <- FindVariableFeatures(scrna, selection.method = "vst", nfeatures = 3000)
top3000 <- head(VariableFeatures(scrna), 3000)
scrna <- scrna[top3000] #use only top 2000 var genes

adata <- sc$AnnData(
  X   = t(as.matrix(GetAssayData(scrna, slot='counts'))),
  obs = scrna[[]]
)

print(adata)

# run seteup_anndata, use column stim for batch
scvi$data$setup_anndata(adata, batch_key = 'var_sample')

# create the model
model = scvi$model$SCVI(adata)

# train the model
model$train()

latent = model$get_latent_representation()

# put it back in our original Seurat object
latent <- as.matrix(latent)
rownames(latent) = colnames(scrna)
scrna[['scvi']] <- CreateDimReducObject(embeddings = latent, key = "scvi_", assay = DefaultAssay(scrna))

#Plot data
scrna <- RunUMAP(scrna, dims = 1:10, reduction = 'scvi', n.components = 2)
saveRDS(scrna,"~/Desktop/scrna_scvi.RDS")
scrna2 = readRDS("~/Desktop/IntergratedSeurat_I73T.RDS")

p1 <- DimPlot(scrna, reduction = "umap", group.by = "var_sample", pt.size=1)
p2 <- DimPlot(scrna, reduction = "umap", split.by = "var_sample", pt.size=1)
p3 <- DimPlot(scrna2, reduction = "umap", group.by = 'orig.ident', pt.size=1)
p4 <- DimPlot(scrna2, reduction = "umap", group.by = 'orig.ident',split.by = "orig.ident", pt.size=1)
plot_grid(p1,p2,p3,p4, nrow = 2, labels = c("SCVI","","Seurat","")) %>% ggsave(file="~/Desktop/scvi_seurat.png",height=10, width=15, units= "in", dpi = 300)



