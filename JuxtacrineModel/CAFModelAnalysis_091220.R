
### Analyze results from CAF model

## Set up
library(ggplot2)

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Import and compile data

csv.files = list.files(pattern="\\.csv$")

smmry.df = as.data.frame(matrix(nrow=length(csv.files),ncol=6))
colnames(smmry.df) = c("log.size","caf.frac","bdg.max","iter","end.rchd","frac.R","maj.R")

for (i in 1:length(csv.files)) {
  
  file.i = csv.files[i]
  
  # Read file
  df.i = read.csv(file.i,header=F)
  colnames(df.i) = c("j","t","nS","nR","nA")
  
  # Metadata
  df.splt = strsplit(file.i,"_")[[1]]
  smmry.df$caf.frac[i] = as.numeric(gsub("caffrac","",df.splt[2]))
  smmry.df$bdg.max[i] = as.numeric(gsub("bdgmax","",df.splt[3]))
  smmry.df$log.size[i] = as.numeric(gsub("logsize","",df.splt[4]))
  smmry.df$iter[i] = as.numeric(gsub("iter","",df.splt[5]))
  
  # Summarize data
  smmry.df$end.rchd[i] = df.i$t[nrow(df.i)]>=365   # did the simulation get to 1 year?
  
  fin.pop = df.i[nrow(df.i),3:5]
  smmry.df$frac.R[i] = fin.pop$nR/sum(fin.pop) # relative frequency of R pop
  
  smmry.df$maj.R[i] = fin.pop$nR>fin.pop$nA # was genetic resistance more abundant than microenvironmental?
  print(i)
}

## Aggregate data

agg.df = aggregate(smmry.df[,5:ncol(smmry.df)],
                   list(log.size = smmry.df$log.size,
                        bdg.max = smmry.df$bdg.max,
                        caf.frac = smmry.df$caf.frac),
                   FUN = function(x) c(mean=mean(x,na.rm=T),sem=sd(x,na.rm=T)/sqrt(sum(!is.na(x)))) )
agg.df2 = do.call("data.frame",agg.df)

## Plot results

plt.df = agg.df2[agg.df2$log.size==4,]
plt.df$bdg.max = factor(plt.df$bdg.max,levels=c(20,5,2,1))

# Population structure
ggplot(plt.df,aes(x=as.factor(caf.frac),y=bdg.max,fill=frac.R.mean))+theme_bw()+
  geom_raster()+
  ggtitle("Population Structure:\nGenetic Resistance")+
  xlab("CAF Infiltration")+
  ylab("Maximum Budging Distance")+
  scale_fill_gradientn("",colors=viridis::viridis(10),
                       breaks=c(0,0.5,1),labels=c("0%","50%","100%"),
                       limits = c(0,1))+
  scale_x_discrete(labels=c("0%","1%","5%","10%","25%"))+
  theme(
    plot.title = element_text(size=26,face="bold",hjust=0.5),
    axis.title = element_text(size=24,face="bold"),
    axis.text = element_text(size=20,color="black",face="bold"),
    legend.title = element_text(size=18,hjust=0.5,face="bold"),
    legend.text = element_text(size=16,face="bold")
  )

# ggsave("CAFModelHeatMap.png",width=7.5,height=6)

# Population structure (nonresistant)
ggplot(plt.df,aes(x=as.factor(caf.frac),y=bdg.max,fill=1-frac.R.mean))+theme_bw()+
  geom_raster()+
  ggtitle("Population Structure:\nNon-genetic Resistance")+
  xlab("CAF Infiltration")+
  ylab("Maximum Budging Distance")+
  scale_fill_gradientn("",colors=viridis::viridis(10),
                       breaks=c(0,0.5,1),labels=c("0%","50%","100%"),
                       limits = c(0,1))+
  scale_x_discrete(labels=c("0%","1%","5%","10%","25%"))+
  theme(
    plot.title = element_text(size=26,face="bold",hjust=0.5),
    axis.title = element_text(size=24,face="bold"),
    axis.text = element_text(size=20,color="black",face="bold"),
    legend.title = element_text(size=18,hjust=0.5,face="bold"),
    legend.text = element_text(size=16,face="bold")
  )
# ggsave("CAFModelnonRHeatMap.png",width=7.5,height=6)


# Majority resistance?
ggplot(plt.df,aes(x=as.factor(caf.frac),y=bdg.max,fill=maj.R.mean))+theme_bw()+
  geom_raster()+
  ggtitle("Population Structure At One Year")+
  xlab("CAF Infiltration")+
  ylab("Maximum Budging Distance")+
  scale_fill_gradientn("Genetic\nResistance",colors=viridis::viridis(10),
                       breaks=c(0.25,0.5,0.75),labels=c("25%","50%","75%"))+
  scale_x_discrete(labels=c("0%","1%","5%","10%"))+
  theme(
    plot.title = element_text(size=26,face="bold",hjust=0.5),
    axis.title = element_text(size=24,face="bold"),
    axis.text = element_text(size=20,color="black"),
    legend.title = element_text(size=18,hjust=0.5,face="bold"),
    legend.text = element_text(size=16,face="bold")
  )
