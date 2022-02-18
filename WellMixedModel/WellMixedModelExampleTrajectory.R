### Plot median and 5-95 percentile trajectories of well mixed models

## Set up
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(reshape2)

## Import and analyze data

pop.max = 1e10

csv.files = list.files(pattern="\\.csv$")



for (i in 1:length(csv.files)) {
  
  df.i = read.csv(csv.files[i],header = F)
  colnames(df.i) = c("t","nS","nR","nA","nC")
  
  n.trt = which.max(rowSums(df.i[,c("nS","nA")]))
  n.rlps = which(rowSums(df.i[,c("nR","nA")])>=pop.max)[1]
  
  PFS = df.i$t[n.rlps] - df.i$t[n.trt]
  
  smmry.df$pat[i] = i
  smmry.df$PFS[i] = PFS
  
}

## Intermediate steps

# Sort simulation trajectories by PFS
smmry.df = smmry.df[order(smmry.df$PFS),]
rownames(smmry.df) = 1:nrow(smmry.df)

# Note 5th, 50th, and 95th percentile
df.5 = read.csv(csv.files[smmry.df$pat[5]])
colnames(df.5) =  c("t","nS","nR","nA","nC")
df.50 = read.csv(csv.files[smmry.df$pat[50]])
colnames(df.50) =  c("t","nS","nR","nA","nC")
df.95 = read.csv(csv.files[smmry.df$pat[95]])
colnames(df.95) =  c("t","nS","nR","nA","nC")

## Plot

# Tidy data frames

td.5 = melt(df.5,id.vars = "t")
td.5 = td.5[td.5$variable%in%c("nS","nR","nA"),]
td.5 = td.5[td.5$value>0,]

td.50= melt(df.50,id.vars = "t")
td.50 = td.50[td.50$variable%in%c("nS","nR","nA"),]
td.50 = td.50[td.50$value>0,]

td.95 = melt(df.95,id.vars = "t")
td.95 = td.95[td.95$variable%in%c("nS","nR","nA"),]
td.95 = td.95[td.95$value>0,]

# Plot results

ggplot(data=td.50,aes(x=t/30.4,y=log10(value),color=variable))+theme_bw()+
  geom_line(size=2)+
  scale_y_continuous(breaks=c(0,3,6,9))+
  scale_color_manual(values=c("#2ABAFC","#EE1D23","#216280"))+
  xlab("Time (months)")+
  ylab("Population Size")+
  ggtitle("Well-Mixed Model of\nResistance Evolution")+
  theme(
    plot.title = element_text(size=20,face="bold",hjust=0.5),
    axis.title = element_text(size=16,face="bold"),
    axis.text = element_text(size=14,face="bold",color="black")
  )+
  guides(color="none")

ggsave("WellMixedTrajectory.png",width=5,height=4)
