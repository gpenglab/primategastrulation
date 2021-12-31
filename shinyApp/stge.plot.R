library(ggplot2)
library(plotly)
library(data.table)
library(shiny)
library(plotly)
library(dplyr)

base <- 'data/'

### LOAD data ###
express_vals <- readRDS(paste0(base,'Mfas.exp.meta.rds'))
gene.names<-colnames(express_vals)[-c(1:6)]
gene.co<-colnames(express_vals)[-c(1:18)]

# express_vals[["TBXT"]]

shp.val <-c("Em-Up"=16,
        "Em-Mid"=16,
        "Em-Low"=16,
        "EXMC1"=18,
        "EXMC2"=18,
        "YE1"=17,
        "YE2"=17
)

plot.std.col <- function(df, xname, yname, title){
  
  ggplot(df)+
    geom_point(size=4,aes(x=x, y=y, color=z, shape=shp))+
    xlab(xname)+
    ylab(yname)+
    ggtitle(title)+
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(linetype = "solid", fill = NA),
          panel.background = element_blank(),
          legend.text=element_text(size=8),#size of legend
          legend.title=element_text(size=12),
          axis.text.x = element_blank(),
          axis.text.y = element_blank())  +
    scale_colour_gradient2(low ="darkgreen", mid = "yellow", high = "red",name= "log2(TPM+1)")+ #forestgreen,mid = "palegoldenrod",
    scale_shape_manual(values=shp.val,name="Lineage")
}


plot_2gene<-function(gene.1,gene.2){
  data.1<-unlist(express_vals[[gene.1]])
  data.2<-unlist(express_vals[[gene.2]])
  data.1.norm<-(data.1-min(data.1))/(max(data.1)-min(data.1))
  data.2.norm<-(data.2-min(data.2))/(max(data.2)-min(data.2))
  col.g1<-c(103,39,112)/255#c(30/255,144/255,255/255)
  col.g2<-c(col2rgb("darkorange"))/255   #c(105,154,51)/255#c(1,0,0)
  n<-c(1,1,1)
  f<-c(0,0.961,0.961) #cyan
  A=col.g1[1]-n[1]
  B=col.g2[1]-n[1]
  C=n[1]
  D=col.g1[2]-n[2]
  E=col.g2[2]-n[2]
  Ff=n[2]
  G=col.g1[3]-n[3]
  H=col.g2[3]-n[3]
  I=n[3]
  alpha<-f[1]-A-B-C
  beta<-f[2]-D-E-Ff
  gamma<-f[3]-G-H-I
  data.norm<-cbind(data.1.norm, data.2.norm)
  z1<-apply(data.norm, 1, function(x) max(c(0,A*x[1]+B*x[2]+C+alpha*x[1]*x[2])))
  z2<-apply(data.norm, 1, function(x) max(c(0,D*x[1]+E*x[2]+Ff+beta*x[1]*x[2])))
  z3<-apply(data.norm, 1, function(x) max(c(0,G*x[1]+H*x[2]+I+gamma*x[1]*x[2])))
  df.3<-data.frame(shp=express_vals$cluster,x=express_vals$X, y=express_vals$Y,
                   z1=z1,#0.5*(1+data.1.norm),
                   z2=z2,#0.5-max(c(0.5*(1+data.1.norm), 0.5*(1+data.2.norm))),#rep(0,length(data.1.norm))
                   z3=z3)
  
  p3<-ggplot()+
    xlab("Distal")+
    ylab("")+
    theme(plot.title = element_text(size=2, face="bold"))+
    geom_point(data=df.3,size=4,aes(x=x, y=y, color=rgb(z1,z2, z3), shape=shp))+
    theme(panel.background = element_rect(fill = 'grey80', colour = 'white'),
          plot.title = element_text(hjust = 0.5),
          panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(linetype = "solid", fill = NA),
          legend.text=element_text(size=8),#size of legend
          legend.title=element_text(size=12),
          plot.margin=unit(c(1,1,1.5,1.2),"cm"),
          axis.text.x = element_blank(),
          axis.text.y = element_blank())  +
    scale_color_identity()+
    scale_shape_manual(values=shp.val,name="Lineage")
   
  p3
  
}
plot_2genecolor<-function(gene.1,gene.2){
  data.1<-unlist(express_vals[[gene.1]])
  data.2<-unlist(express_vals[[gene.2]])
  data.1.norm<-(data.1-min(data.1))/(max(data.1)-min(data.1))
  data.2.norm<-(data.2-min(data.2))/(max(data.2)-min(data.2))
  col.g1<-c(103,39,112)/255#c(30/255,144/255,255/255)
  col.g2<-c(col2rgb("darkorange"))/255   #c(105,154,51)/255#c(1,0,0)
  n<-c(1,1,1)
  f<-c(0,0.961,0.961) #cyan
  A=col.g1[1]-n[1]
  B=col.g2[1]-n[1]
  C=n[1]
  D=col.g1[2]-n[2]
  E=col.g2[2]-n[2]
  Ff=n[2]
  G=col.g1[3]-n[3]
  H=col.g2[3]-n[3]
  I=n[3]
  alpha<-f[1]-A-B-C
  beta<-f[2]-D-E-Ff
  gamma<-f[3]-G-H-I
  data.norm<-cbind(data.1.norm, data.2.norm)
  z1<-apply(data.norm, 1, function(x) max(c(0,A*x[1]+B*x[2]+C+alpha*x[1]*x[2])))
  z2<-apply(data.norm, 1, function(x) max(c(0,D*x[1]+E*x[2]+Ff+beta*x[1]*x[2])))
  z3<-apply(data.norm, 1, function(x) max(c(0,G*x[1]+H*x[2]+I+gamma*x[1]*x[2])))
  df.4<-data.frame(x=express_vals$X, y=express_vals$Y, 
                   z1=z1,#0.5*(1+data.1.norm),
                   z2=z2,#0.5-max(c(0.5*(1+data.1.norm), 0.5*(1+data.2.norm))),#rep(0,length(data.1.norm))
                   z3=z3)
  x<-seq(0,1,0.1)
  y<-seq(0,1,0.1)
  m<-expand.grid(x,y)
  res<-t(apply(m, 1, function(x) c(max(c(0,A*x[1]+B*x[2]+C+alpha*x[1]*x[2])), 
                                   max(c(0,D*x[1]+E*x[2]+Ff+beta*x[1]*x[2])),
                                   max(c(0,G*x[1]+H*x[2]+I+gamma*x[1]*x[2]))))
  )
  id<-which(res>1,arr.ind=T)
  res[id]<-1
  df<-data.frame(m=m,col=res)
  
  p4<-ggplot(df, aes(x=m.Var1,y=m.Var2,fill=rgb(col.1,col.2,col.3)))+
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "white"))+
    xlab(gene.1)+
    ylab(gene.2)+
    geom_raster()+
    scale_fill_identity()
  p4
  
}




