### Analyze simulation results for ECM model

## Set up
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(viridis)

## Import and summarize simulations

csv.files = list.files(pattern="\\.csv$")

smmry.df = as.data.frame(matrix(nrow=length(csv.files),ncol=10))
colnames(smmry.df) = c("caf.frac","bdg.max","mgrt","log.size","iter",
                       "S.fin","R.fin","A.fin",
                       "R.frac","A.frac")

for (i in 1:length(csv.files)) {
  
  file.i = csv.files[i]
  
  # Metadata
  df.splt = strsplit(file.i,"_")[[1]]
  smmry.df$caf.frac[i] = as.numeric(gsub("caffrac","",df.splt[2]))
  smmry.df$bdg.max[i] = as.numeric(gsub("bdgmax","",df.splt[3]))
  smmry.df$mgrt[i] = as.numeric(gsub("mgrt","",df.splt[4]))
  smmry.df$log.size[i] = as.numeric(gsub("logsize","",df.splt[5]))
  smmry.df$iter[i] = as.numeric(gsub("iter","",df.splt[6]))
  
  df.i = read.csv(file.i,header=F)
  colnames(df.i) = c("n","t","nS","nR","nA")
  
  # Record final makeup of tumor
  pop.fin = df.i[nrow(df.i),c("nS","nR","nA")]
  smmry.df[i,c("S.fin","R.fin","A.fin")] = pop.fin
  smmry.df$R.frac[i] = pop.fin[2]/sum(pop.fin)
  smmry.df$A.frac[i] = pop.fin[3]/sum(pop.fin)
  
  print(i)
  
}

smmry.df$R.frac = unlist(smmry.df$R.frac)
smmry.df$A.frac = unlist(smmry.df$A.frac)


## Aggregate data

agg.df = do.call(data.frame,
                 aggregate(.~caf.frac+bdg.max+mgrt+log.size, smmry.df, function(x) c(mean = mean(x), sem = sd(x)/sqrt(length(x)))))

agg.df$bdg.max = factor(agg.df$bdg.max,levels=c(20,5,2,1))

## Plot data

# Vary migration rate
gradient.limits = c(min(agg.df$R.fin.mean),max(agg.df$R.fin.mean))

ggplot(agg.df[agg.df$mgrt==0,],aes(x=as.factor(caf.frac),y=as.factor(bdg.max),fill=R.fin.mean))+theme_bw()+
  geom_raster()+
  scale_fill_gradientn(colors=viridis(5),limits=gradient.limits)+
  xlab("CAF Infiltration")+
  ylab("Maximum Budging Distance")+
  ggtitle("Genetic Resistance:\nNo CAF Migration")+
  scale_x_discrete(labels=c("0%","1%","5%","10%","25%"))+
  theme(
    plot.title = element_text(size=26,face="bold",hjust=0.5),
    axis.title = element_text(size=24,face="bold"),
    axis.text = element_text(size=20,color="black",face="bold"),
    legend.title = element_text(size=18,hjust=0.5,face="bold"),
    legend.text = element_text(size=16,face="bold")
  )
# ggsave("NoCAFMigrationHeatMap.png",width=7.5,height=6)

ggplot(agg.df[agg.df$mgrt==0.5,],aes(x=as.factor(caf.frac),y=as.factor(bdg.max),fill=R.fin.mean))+theme_bw()+
  geom_raster()+
  scale_fill_gradientn(colors=viridis(5),limits=gradient.limits)+
  xlab("CAF Infiltration")+
  ylab("Maximum Budging Distance")+
  ggtitle("Genetic Resistance:\nCAF Migration")+
  scale_x_discrete(labels=c("0%","1%","5%","10%","25%"))+
  theme(
    plot.title = element_text(size=26,face="bold",hjust=0.5),
    axis.title = element_text(size=24,face="bold"),
    axis.text = element_text(size=20,color="black",face="bold"),
    legend.title = element_text(size=18,hjust=0.5,face="bold"),
    legend.text = element_text(size=16,face="bold")
  )
# ggsave("CAFMigrationHeatMap.png",width=7.5,height=6)


ggplot(agg.df[agg.df$mgrt==1,],aes(x=as.factor(caf.frac),y=as.factor(bdg.max),fill=R.fin.mean))+theme_bw()+
  geom_raster()+
  scale_fill_gradientn(colors=viridis(5),limits=gradient.limits)

# Alt: plot fraction of final pop

gradient.limits = c(0,100)

ggplot(agg.df[agg.df$mgrt==0,],aes(x=as.factor(caf.frac),y=as.factor(bdg.max),fill=R.frac.mean*100))+theme_bw()+
  geom_raster()+
  scale_fill_gradientn(colors=viridis(5),limits=gradient.limits)+
  xlab("CAF Infiltration")+
  ylab("Maximum Budging Distance")+
  ggtitle("Genetic Resistance:\nNo CAF Migration")+
  scale_x_discrete(labels=c("0%","1%","5%","10%","25%"))+
  theme(
    plot.title = element_text(size=26,face="bold",hjust=0.5),
    axis.title = element_text(size=24,face="bold"),
    axis.text = element_text(size=20,color="black"),
    legend.title = element_text(size=18,hjust=0.5,face="bold"),
    legend.text = element_text(size=16,face="bold")
  )

ggplot(agg.df[agg.df$mgrt==0.5,],aes(x=as.factor(caf.frac),y=as.factor(bdg.max),fill=R.frac.mean*100))+theme_bw()+
  geom_raster()+
  scale_fill_gradientn(colors=viridis(5),limits=gradient.limits)+
  xlab("CAF Infiltration")+
  ylab("Maximum Budging Distance")+
  ggtitle("Genetic Resistance:\nNo CAF Migration")+
  scale_x_discrete(labels=c("0%","1%","5%","10%","25%"))+
  theme(
    plot.title = element_text(size=26,face="bold",hjust=0.5),
    axis.title = element_text(size=24,face="bold"),
    axis.text = element_text(size=20,color="black"),
    legend.title = element_text(size=18,hjust=0.5,face="bold"),
    legend.text = element_text(size=16,face="bold")
  )

# Alt: R.fin as a function of max budging

ggplot(agg.df[agg.df$mgrt==0,],aes(x=bdg.max,color=as.factor(caf.frac)))+theme_bw()+
  geom_point(aes(y=R.fin.mean),size=2.5)+
  geom_line(aes(y=R.fin.mean),size=1.5)+
  geom_errorbar(aes(ymin=(R.fin.mean-R.fin.sd/sqrt(48)),ymax=(R.fin.mean+R.fin.sd/sqrt(48))),size=1)
  

ggplot(agg.df[agg.df$mgrt==0.5,],aes(x=bdg.max,color=as.factor(caf.frac)))+theme_bw()+
  geom_point(aes(y=R.fin.mean),size=2.5)+
  geom_line(aes(y=R.fin.mean),size=1.5)+
  geom_errorbar(aes(ymin=(R.fin.mean-R.fin.sd/sqrt(48)),ymax=(R.fin.mean+R.fin.sd/sqrt(48))),size=1)

