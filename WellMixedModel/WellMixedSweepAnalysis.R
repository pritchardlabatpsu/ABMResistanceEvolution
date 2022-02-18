### Plot outcomes of sweep of recruitment parameter for well-mixed models

## Set up
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(reshape2)

## Summarize files

csv.files = list.files(pattern = "\\.csv$")

smmry.df = as.data.frame(matrix(nrow=length(csv.files),ncol=5))
colnames(smmry.df) = c("nom.adj.frac","pat","fin.adj.frac","adj.dom","res.dom")

for (i in 1:length(csv.files)) {
  
  df.i = read.csv(csv.files[i],header=F)
  colnames(df.i) = c("pat","adj.frac","t","nS","nR","nA","nC")
  
  smmry.df$nom.adj.frac[i] = df.i$adj.frac[1] # Input to simulation
  smmry.df$pat[i] = df.i$pat[1] # patient ID
  smmry.df$fin.adj.frac[i] = df.i$nA[2]/(df.i$nR[2]+df.i$nA[2]) # what fraction of resistant cells at the end of treatment were adj
  smmry.df$fin.res.frac[i] = df.i$nR[2]/(df.i$nR[2]+df.i$nA[2]) # what fraction of resistant cells were res
  
  adj.dom = df.i$nA[2]>df.i$nR[2] # which phenotype is dominant at relapse
  smmry.df$adj.dom[i] = adj.dom
  smmry.df$res.dom[i] = !adj.dom

}

## Aggregate summarized data by model input (nom.adj.frac)

agg.df = do.call(data.frame,
                 aggregate(.~nom.adj.frac, smmry.df, function(x) c(mean = mean(x), sd = sd(x))))

## Plot

# Reshape data

agg.mlt = melt(agg.df,id.vars = "nom.adj.frac")

# Plot mean phenotype frequency
agg.mlt1 = agg.mlt[agg.mlt$variable%in%c("fin.adj.frac.mean","fin.res.frac.mean"),]
agg.mlt1$variable = factor(agg.mlt1$variable,levels = c("fin.res.frac.mean","fin.adj.frac.mean"))

ggplot(agg.mlt1,aes(x=as.factor(log10(nom.adj.frac)),y=value,fill=variable))+
  theme_bw()+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("#EE1D23","#216280"))+
  scale_y_continuous(breaks=c(0,0.5,1))+
  xlab("Percentage of\nCAF-Mediated Resistant Cells")+
  ylab("Final Resistance\nPopulation Structure")+
  ggtitle("Parameter Sweep:\nCAF Recruitment")+
  theme(
    plot.title = element_text(size=20,face="bold",hjust=0.5),
    axis.title = element_text(size=16,face="bold"),
    axis.text = element_text(size=14,face="bold",color="black")
  )+
  guides(fill="none")

# ggsave("ResistanceOutcomes.png",width=5,height=4)

# Alt: Plot dominant phenotypes across cohort

agg.mlt2 = agg.mlt[agg.mlt$variable%in%c("adj.dom.mean","res.dom.mean"),]
ggplot(agg.mlt2,aes(x=as.factor(log10(nom.adj.frac)),y=value,fill=variable))+
  theme_bw()+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("#216280","#EE1D23"))+
  scale_y_continuous(breaks=c(0,0.5,1))+
  xlab("Percentage of\nCAF-Mediated Resistant Cells")+
  ylab("Frequency as Dominant\nResistance Phenotype")+
  ggtitle("Resistance Evolutionary Outcomes")+
  theme(
    plot.title = element_text(size=22,face="bold",hjust=0.5),
    axis.title = element_text(size=20,face="bold"),
    axis.text = element_text(size=18,face="bold",color="black")
  )+
  guides(fill="none")
