### Analyze Paracrine CAF Model

## Set up
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(viridis)

## Read in files

csv.files = list.files(pattern="logsize4_iter[0-9]{3}_122721\\.csv$")

cnames = c("caf.frac","bdg.max","caf.rad","log.size","iter","S.eval","R.eval","A.eval")
rslts.df = as.data.frame(matrix(nrow=length(csv.files),ncol=length(cnames)))
colnames(rslts.df) = cnames

# t.eval = 360 # choose time point to evaluate

# Load results if available; if not, run for loop below
# load("Results.RData")

for (i in 1:length(csv.files)) {

  file.i = csv.files[i]

  # Identify simulation parameters
  file.splt = strsplit(file.i,"_")[[1]]
  params = as.numeric(gsub("[a-zA-Z]","",file.splt[2:6]))
  names(params) = c("caf.frac","bdg.max","caf.rad","log.size","iter")
  rslts.df[i,1:length(params)] = params

  # Read in simulation results
  df.i = read.csv(file.i,header=F)
  colnames(df.i) = c("step","time","S","R","A")

  # Note population structure after t.eval days
  # step.eval = which(df.i$time>=t.eval)[1]
  # pop.eval = as.vector(df.i[step.eval,c("S","R","A")])
  pop.eval = as.vector(df.i[nrow(df.i),c("S","R","A")])
  rslts.df[i,c("S.eval","R.eval","A.eval")] = pop.eval

  print(i)

}

# save(rslts.df,file="Log4Results.RData")

## Summarize data

smmry.df = aggregate(cbind(S.eval,R.eval,A.eval,iter)~.,rslts.df,FUN=mean)
smmry.df = smmry.df[smmry.df$caf.frac!=0.2,]

# FIX THIS
smmry.df$R.eval[smmry.df$caf.frac==0&smmry.df$bdg.max%in%c(5,20)&smmry.df$caf.rad==2] = smmry.df$R.eval[smmry.df$caf.frac==0&smmry.df$bdg.max%in%c(5,20)&smmry.df$caf.rad==1]
smmry.df$R.eval[smmry.df$caf.frac==0&smmry.df$bdg.max%in%c(5,20)&smmry.df$caf.rad==3] = smmry.df$R.eval[smmry.df$caf.frac==0&smmry.df$bdg.max%in%c(5,20)&smmry.df$caf.rad==1]

## Visualize data

# Heatmaps at different caf.rad
caf.rad = 3
log.size = 4

plot.cndtns = smmry.df$caf.rad==caf.rad & smmry.df$log.size==log.size
gradient.limits = c(min(smmry.df$R.eval[smmry.df$log.size==log.size]),max(smmry.df$R.eval[smmry.df$log.size==log.size]))

ggplot(smmry.df[plot.cndtns,],aes(x=factor(caf.frac),y=factor(bdg.max,levels=c(20,5,2,1)),fill=R.eval))+theme_bw()+
  geom_raster()+
  scale_fill_gradientn(colors=viridis(5),limits=gradient.limits)

# Alternative: Plot resistance outgrowth vs bdg.max for range of caf.rad values

agg.df = do.call(data.frame,aggregate(cbind(S.eval,R.eval,A.eval,iter)~.,rslts.df,function(x) c(mean = mean(x), sd = sd(x))))

caf.frac = 0.05
log.size = 4
plot.cndtns = agg.df$caf.frac==caf.frac & agg.df$log.size==log.size

ggplot(agg.df[plot.cndtns,],aes(x=bdg.max,color=factor(caf.rad)))+theme_bw()+
  geom_point(aes(y=R.eval.mean),size=2.5)+
  geom_line(aes(y=R.eval.mean),size=1.5)+
  geom_errorbar(aes(ymin=(R.eval.mean-R.eval.sd/sqrt(48)),ymax=(R.eval.mean+R.eval.sd/sqrt(48))),size=1)+
  xlab("Max Budging Distance")+
  ylab("Resistance Outgrowth")+
  ggtitle("Effect of Paracrine\nSignaling Range")+
  scale_color_manual(values=c("#E88329","#D2DDA4","#0F8899"))+
  theme(
    plot.title = element_text(size=22,face="bold",hjust=0.5),
    axis.title = element_text(size=20,face="bold"),
    axis.text = element_text(size=18,face="bold",color="black")
  )+
  guides(color="none")

# ggsave("ParacrineModelResults.png",width=5,height=4)
