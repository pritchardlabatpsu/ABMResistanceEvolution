
## Analyze spatial competition vs tumor size simulations

## Set up
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(RColorBrewer)

## Import and summarize data

csv.files = list.files(pattern="SpatialConstraint_log")
cnames = c("x.pos","y.pos","z.pos","phn") 

# phn = 1 ==> nondividing
# phn = 3 ==> dividing

results.df = as.data.frame(matrix(nrow=length(csv.files),ncol=5))
colnames(results.df) = c("logsize","bdgmax","iter","rel.dst.nondiv","rel.dst.div")

for (i in 1:length(csv.files)) {
  
  # Get meta data
  file.i = csv.files[i]
  file.splt = unlist(strsplit(file.i,"_"))
  
  results.df$logsize[i] = as.numeric(gsub("logsize","",file.splt[2]))
  results.df$bdgmax[i] = as.numeric(gsub("budgmax","",file.splt[3]))
  results.df$iter[i] = as.numeric(gsub("iter","",file.splt[4]))
  
  
  # Import simulation results
  df.i = read.csv(file.i,header=F)
  colnames(df.i) = cnames
  
  # Calculate centroid of tumor
  cntr.i = colMeans(df.i[,1:3])
  df.i$dst = sqrt((df.i$x.pos-cntr.i[1])^2+
                       (df.i$y.pos-cntr.i[2])^2+
                       (df.i$z.pos-cntr.i[3])^2)
  
  mean.dst.nondiv = mean(df.i$dst[df.i$phn==1]) # average distance from centroid of nondividing cells
  mean.dst.div = mean(df.i$dst[df.i$phn==3])    # average distance from centroid of dividing cells
  max.dst = max(df.i$dst)
  
  results.df$rel.dst.nondiv[i] = mean.dst.nondiv/max.dst
  results.df$rel.dst.div[i] = mean.dst.div/max.dst
  
  print(i)
  
}

## Summarize data

smmry.df = aggregate(cbind(rel.dst.nondiv,rel.dst.div)~logsize+bdgmax,data=results.df,FUN=mean)

## Plot results

# Dividing cells
ggplot(smmry.df,aes(x=logsize,y=as.factor(bdgmax),fill=rel.dst.div))+theme_bw()+
  geom_raster()+
  scale_fill_gradientn("Relative Distance\nof Dividing Cells\nfrom Centroid",colors=c("#E88329","#0F8899"),breaks=c(0.7,0.8,0.9))+
  ggtitle("Spatial Competition\nvs Tumor Size")+
  xlab("Tumor Size")+ylab("Max Budging Distance")+
  theme(
    plot.title = element_text(size=26,face="bold",hjust=0.5),
    axis.title = element_text(size=24,face="bold"),
    axis.text = element_text(size=20,color="black",face="bold"),
    legend.title = element_text(size=18,hjust=0.5,face="bold"),
    legend.text = element_text(size=16,face="bold")
  )
ggsave("SpatialCompetition.png",width=8,height=5)

# Nondividing cells
ggplot(smmry.df,aes(x=logsize,y=as.factor(bdgmax),fill=rel.dst.nondiv))+theme_bw()+
  geom_raster()
