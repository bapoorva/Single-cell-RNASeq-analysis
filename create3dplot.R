library(Seurat)
library(rgl)
require(slingshot)
devtools::source_url('https://raw.githubusercontent.com/mpmorley/scRNA/master/SeuratV3_Extras.R')
library(dplyr)
library(tidyr)
library(plot3D)

load("~/Desktop/labproject/projects/Morrisey/Jarod/scRNA/JZ_lung_timeseries/SeuratV3_integration_E12.5_to_Adult/E15.5_to_P7_epi.RData")
object=scrna.sub
#rm(scrna.sub)
reduction='dm'
groupby='var_cluster'
dims=1:3
dims <- paste0(Key(object = object[[reduction]]), dims)
data <- FetchData(object = object, vars = c(dims,groupby))


if(is.numeric(data[,4])){
  data$color = map2color(data[,4],c(-3, 3))
}else{
  data$color=cpallette[as.numeric(data[,4])]
}

dims <- 1:3
object <- runSlingshot(object)

maxdim<-max(dims)

groupby='curve1'

curved <- bind_rows(lapply(names(object@misc$sds$data@curves), function(x){c <- slingCurves(object@misc$sds$data)[[x]]
d <- as.data.frame(c$s[c$ord,seq_len(maxdim)])
d$curve<-x
return(d)
})
) 

a3=aggregate(data$DM_3, by=list(data$var_cluster), FUN=mean)
a2=aggregate(data$DM_2, by=list(data$var_cluster), FUN=mean)
a1=aggregate(data$DM_1, by=list(data$var_cluster), FUN=mean)
center=inner_join(a1,a2,by="Group.1")
center=inner_join(center,a3,by="Group.1")
colnames(center)=c("var_cluster","x","y","z")
centers=center

p=FetchData(object,append(c('DM_1','DM_2','DM_3'),'curve1'))
p=FetchData(object,append(c('DM_1','DM_2','DM_3'),'curve2'))
p=FetchData(object,append(c('DM_1','DM_2','DM_3'),'curve3'))

pdf("~/Desktop/JZ/curve1.pdf")
scatter3D(data$DM_1,y=data$DM_2,z=data$DM_3, bty="g",phi=0, col = colorRampPalette(brewer.pal(n = 9, 'RdYlBu'))(30),pch=16,colvar=p$curve1,
          xlab = "DM_1",ylab ="DM_2", zlab = "DM_3", NAcol = 'grey')
lines3D(x=curved$DM_1,y=curved$DM_2,z=curved$DM_3, add=TRUE, col="black",lwd=2)
dev.off()

pdf("~/Desktop/JZ/curve2.pdf")
scatter3D(data$DM_1,y=data$DM_2,z=data$DM_3, bty="g",phi=0, col = colorRampPalette(brewer.pal(n = 9, 'RdYlBu'))(30),pch=16,colvar=p$curve2,
          xlab = "DM_1",ylab ="DM_2", zlab = "DM_3", NAcol = 'grey')
lines3D(x=curved$DM_1,y=curved$DM_2,z=curved$DM_3, add=TRUE, col="black",lwd=2)
dev.off()

pdf("~/Desktop/JZ/curve3.pdf")
scatter3D(data$DM_1,y=data$DM_2,z=data$DM_3, bty="g",phi=0, col = colorRampPalette(brewer.pal(n = 9, 'RdYlBu'))(30),pch=16,colvar=p$curve3,
          xlab = "DM_1",ylab ="DM_2", zlab = "DM_3", NAcol = 'grey')
lines3D(x=curved$DM_1,y=curved$DM_2,z=curved$DM_3, add=TRUE, col="black",lwd=2)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
pdf("~/Desktop/test.pdf")
scatter3D(data$DM_1,y=data$DM_2,z=data$DM_3, bty="g", col = gg.col(100),pch=16,
          xlab = "DM_1",ylab ="DM_2", zlab = "DM_3")

text3D(x=center$x,y=center$y,z=center$z,labels=center$var_cluster, add=TRUE,cex=1,adj=5)
lines3D(x=curved$DM_1,y=curved$DM_2,z=curved$DM_3, add=TRUE, col="black",lwd=2)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
pdf("~/Desktop/test2.pdf")
scatter3D(data$DM_1,y=data$DM_2,z=data$DM_3, bty="g", col = data$color,pch=16,col.var=as.integer(data$var_cluster),
          xlab = "DM_1",ylab ="DM_2", zlab = "DM_3")

text3D(x=center$x,y=center$y,z=center$z,labels=center$var_cluster, add=TRUE,cex=1,adj=5)
lines3D(x=curved$DM_1,y=curved$DM_2,z=curved$DM_3, add=TRUE, col="black",lwd=2)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
pdf("~/Desktop/test3.pdf")
scatter3D(data$DM_1,y=data$DM_2,z=data$DM_3, bty="g",phi=0, col = data$color,pch=16,col.var=as.integer(data$var_cluster),
          xlab = "DM_1",ylab ="DM_2", zlab = "DM_3")

text3D(x=center$x,y=center$y,z=center$z,labels=center$var_cluster, add=TRUE,cex=1,adj=5)
lines3D(x=curved$DM_1,y=curved$DM_2,z=curved$DM_3, add=TRUE, col="black",lwd=2)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
pdf("~/Desktop/test4.pdf")
scatter3D(data$DM_1,y=data$DM_2,z=data$DM_3, bty="g",phi=0, col = cpallette,pch=16,colvar=as.integer(data$var_cluster),
          xlab = "DM_1",ylab ="DM_2", zlab = "DM_3")

#text3D(x=center$x,y=center$y,z=center$z,labels=center$var_cluster, add=TRUE,cex=1,adj=5)
lines3D(x=curved$DM_1,y=curved$DM_2,z=curved$DM_3, add=TRUE, col="black")
dev.off()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
pdf("~/Desktop/test5.pdf")
scatter3D(data$DM_1,y=data$DM_2,z=data$DM_3, bty="g",phi=0, col = cpallette,pch=".",colvar=as.integer(data$var_cluster),
          xlab = "DM_1",ylab ="DM_2", zlab = "DM_3")

#text3D(x=center$x,y=center$y,z=center$z,labels=center$var_cluster, add=TRUE,cex=1,adj=5)
lines3D(x=curved$DM_1,y=curved$DM_2,z=curved$DM_3, add=TRUE, col="black",lwd=2)
dev.off()

########################################################################################################################
########################################################################################################################
########################################################################################################################
a <- list()
for (i in 1:nrow(centers)) {
  a[[i]] <- list(x= centers$x[i],y= centers$y[i],z= centers$z[i],text= centers$var_cluster[i],showarrow= T,arrowhead=4,arrowsize=0.5)
}

plot=plot_ly() %>%
  add_trace(x = data$DM_1,y = data$DM_2,z = data$DM_3,colors=cpallette,color=data$var_cluster,size=.5,type = "scatter3d") %>% 
  add_paths(x = curved$DM_1,y = curved$DM_2,z = curved$DM_3, mode="lines",color=I("black"),size=I(7)) %>% 
  layout(
    scene = list(
      aspectratio = list(x = 1,y = 1,z = 1),
      dragmode = "turntable",
      xaxis = list(title = dims[1]),
      yaxis = list(title = dims[2]),
      zaxis = list(title = dims[3]),
      annotations = a
    )
  )

plotly_IMAGE(plot,width = 800,height = 800,out_file = "~/Desktop/image3.png")

########################################################################################################################
########################################################################################################################
########################################################################################################################
map2color<-function(x, limits=NULL){
  inf_vals = is.infinite(x)
  
  if (is.null(limits) == FALSE){
    x[x < limits[1]] = limits[1]
    x[x > limits[2]] = limits[2]
  }
  x[inf_vals] = median(x) 
  ii <- cut(x, breaks = seq(min(x, na.rm=T), max(x, na.rm=T), len = 100), 
            include.lowest = TRUE)
  #colors <- colorRampPalette(c("darkgray", "purple"))(99)[ii]
  colors = viridisLite::viridis(99)[ii]
  colors[inf_vals] = "darkgray"
  return(colors)
}



rgl_init <- function(new.device = FALSE, bg = "white", width = 640) { 
  if( new.device | rgl.cur() == 0 ) {
    rgl.open()
    par3d(windowRect = 50 + c( 0, 0, width, width ) )
    rgl.bg(color = bg )
  }
  rgl.clear(type = c("shapes", "bboxdeco"))
  rgl.viewpoint(theta = 15, phi = 20, zoom = 0.7)
}



rgl_init()
plot3d(data[1:3],
       type='s',
       box=FALSE,
       col=data$color,
       radius = .001
)
text3d(x=center$x,y=center$y,z=center$z,text=center$var_cluster, fontweight = "bold",adj=4.5, color=cpallette)
lines3d(x=curved$DM_1,y=curved$DM_2,z=curved$DM_3, add=TRUE, color="black", lwd=2)

tmp <- data %>% group_by(var_cluster,color) %>% summarise(dm1=mean(DM_1), dm2=mean(DM_2),dm3=mean(DM_3))# %>%
with(tmp,text3d(dm1,dm2,dm3,var_cluster))

movie3d(spin3d(axis = c(0, 1, 1)), duration = 10,
        dir = getwd())

########################################################################################################################
########################################################################################################################
########################################################################################################################
plotPseudoTime = function(object,groupby,reduction='dm',dims=1:2){
  
  curved <- bind_rows(lapply(names(object@misc$sds$data@curves), function(x){c <- slingCurves(object@misc$sds$data)[[x]]
  d <- as.data.frame(c$s[c$ord,seq_len(2)])
  d$curve<-x
  return(d)
  })
  ) 
  
  dims <- paste0(Key(object = object[[reduction]]), dims)
  
  p=FetchData(object = object, vars = c(dims,groupby)) %>%
    ggplot(.,aes_string(x=dims[1],y=dims[2]))+geom_point(aes(color=!!sym(groupby))) + 
    theme(legend.position="top") + 
    guides(col = guide_legend(nrow = 2)) +
    geom_path(aes_string(dims[1], dims[2],linetype="curve"),curved,size=1)
  p
  }
