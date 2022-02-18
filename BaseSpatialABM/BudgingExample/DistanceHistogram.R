### Plot the distances of mitotically active cells from tumor centroid in budging simulations

## Set up
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)

## Import data
cnames = c("x.pos","y.pos","z.pos","phn")

bdg1.df = read.csv("BudgingExample_bdgmax1.csv",header=F)
colnames(bdg1.df) = cnames

bdgInf.df = read.csv("BudgingExample_bdgmaxInf.csv",header=F)
colnames(bdgInf.df) = cnames

## Intermediate calculations
bdg1.cntr = colMeans(bdg1.df[,1:3])
bdg1.df$dst = sqrt((bdg1.df$x.pos-bdg1.cntr[1])^2+
                   (bdg1.df$y.pos-bdg1.cntr[2])^2+
                   (bdg1.df$z.pos-bdg1.cntr[3])^2)

bdgInf.cntr = colMeans(bdgInf.df[,1:3])
bdgInf.df$dst = sqrt((bdgInf.df$x.pos-bdgInf.cntr[1])^2+
                     (bdgInf.df$y.pos-bdgInf.cntr[2])^2+
                     (bdgInf.df$z.pos-bdgInf.cntr[3])^2)

## Plot histograms

stt.clr = "#2ABAFC"
div.clr = "#FFBA08"


# bdgmax = 1
scale.mit = 50

bdg1.stat = bdg1.df[bdg1.df$phn==1,]
bdg1.mit = bdg1.df[bdg1.df$phn==3,]
bdg1.mitr = do.call("rbind", replicate(scale.mit, bdg1.mit, simplify = FALSE))

ggplot()+theme_bw()+
  geom_histogram(data=bdg1.df[bdg1.df$phn==1,],
                 aes(x=dst),binwidth=1,
                 fill=stt.clr,alpha=0.9)+
  geom_histogram(data = bdg1.mitr,
                 aes(x=dst),binwidth=1,
                 fill=div.clr,alpha=0.9)+
  ggtitle("Distribution of Mitotic Activity")+
  scale_y_continuous(
    name = "Static Cells",
    sec.axis = sec_axis( trans=~./scale.mit, name="Dividing Cells")
  )+
  scale_x_continuous(
    name = "Distance from Centroid"
  )+
  theme(
    plot.title = element_text(size=20,face="bold",color="black",hjust=0.5),
    axis.title = element_text(size=16,face="bold"),
    axis.title.y.left = element_text(color=stt.clr),
    axis.title.y.right = element_text(color=div.clr),
    axis.text = element_text(size=14,face="bold"),
    axis.text.y.left = element_text(color=stt.clr),
    axis.text.y.right = element_text(color=div.clr),
  )

ggsave("Histogram_bdgmax1.png",width=5,height=4)

# bdgmax = Inf
scale.mit = 50

bdgInf.stat = bdgInf.df[bdgInf.df$phn==1,]
bdgInf.mit = bdgInf.df[bdgInf.df$phn==3,]
bdgInf.mitr = do.call("rbind", replicate(scale.mit, bdgInf.mit, simplify = FALSE))

ggplot()+theme_bw()+
  geom_histogram(data=bdgInf.df[bdgInf.df$phn==1,],
                 aes(x=dst),binwidth=1,
                 fill=stt.clr,alpha=0.9)+
  geom_histogram(data = bdgInf.mitr,
                 aes(x=dst),binwidth=1,
                 fill=div.clr,alpha=0.9)+
  ggtitle("Distribution of Mitotic Activity")+
  scale_y_continuous(
    name = "Static Cells",
    sec.axis = sec_axis( trans=~./scale.mit, name="Dividing Cells")
  )+
  scale_x_continuous(
    name = "Distance from Centroid"
  )+
  theme(
    plot.title = element_text(size=20,face="bold",color="black",hjust=0.5),
    axis.title = element_text(size=16,face="bold"),
    axis.title.y.left = element_text(color=stt.clr),
    axis.title.y.right = element_text(color=div.clr),
    axis.text = element_text(size=14,face="bold"),
    axis.text.y.left = element_text(color=stt.clr),
    axis.text.y.right = element_text(color=div.clr),
  )

ggsave("Histogram_bdgmaxInf.png",width=5,height=4)

