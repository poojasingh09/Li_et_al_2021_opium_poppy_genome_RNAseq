{\rtf1\ansi\ansicpg1252\cocoartf1671\cocoasubrtf400
{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;\red0\green0\blue0;\red255\green255\blue255;\red36\green38\blue41;
\red235\green236\blue237;\red104\green26\blue29;}
{\*\expandedcolortbl;;\cssrgb\c0\c0\c0;\cssrgb\c100000\c100000\c100000;\cssrgb\c18824\c20000\c21176;
\cssrgb\c93725\c94118\c94510;\cssrgb\c49020\c15294\c15294;}
\paperw11900\paperh16840\margl1440\margr1440\vieww25400\viewh13620\viewkind0
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0

\f0\fs24 \cf0 #### P. Singh Aug 2019 ####\
#### FIgure 3b Opium poppy manuscript (Liu et al) ####\
\
\
\
#######load required libraries #####\
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0
\cf2 \cb3 \CocoaLigature0 options(stringsAsFactors = FALSE)\
setwd(".")\
library(plyr)\
library(ggplot2)\
\
\
\
\
\
\
################# generate plotting data for distances 10 mb or less ##############\
\
#read in pathway allocations, make sure header is path1 path2 dist cor compare\
path10less <- read.table("all.merged.path.10mbpless.rho2", header=T)\
\
\
#read in cor and distance data\
\
decay10less <- read.table("all.merged.dup.NA.sorted.10Mb.rho2.txt", header=T)\
\
\
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\tx14960\pardirnatural\partightenfactor0
\cf2 #prep background data\
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0
\cf2 rng1 <- c(0,10,100,1000,10000,20000,30000,40000,50000,60000,70000,80000,90000,100000,200000,300000,400000,500000)\
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\tx14960\pardirnatural\partightenfactor0
\cf2 rng2 <- seq(600000,10000000, by=200000)\
rng <- append(rng1, rng2)\
labels <- rng[1:65]\
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0
\cf2 labels1 <- labels+1 #to avoid inf on log scale\
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\tx14960\pardirnatural\partightenfactor0
\cf2 \
newdata<-ddply(decay10less,.(cut(value.x,breaks=rng,labels=labels)),summarise,q05=quantile(value.y, 0.05))\
newdata1<-ddply(decay10less,.(cut(value.x,breaks=rng,labels=labels)),summarise,q25=quantile(value.y, 0.25))\
newdata2<-ddply(decay10less,.(cut(value.x,breaks=rng,labels=labels)),summarise,q5=quantile(value.y, 0.50))\
newdata3<-ddply(decay10less,.(cut(value.x,breaks=rng,labels=labels)),summarise,q75=quantile(value.y, 0.75))\
newdata4<-ddply(decay10less,.(cut(value.x,breaks=rng,labels=labels)),summarise,q95=quantile(value.y, 0.95))\
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0
\cf2 newdata5<-ddply(decay10less,.(cut(value.x,breaks=rng,labels=labels)),summarise,mean=mean(value.y))\
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\tx14960\pardirnatural\partightenfactor0
\cf2 \
a <- Reduce(function(x, y) merge(x, y, all=TRUE), list(newdata, newdata1, newdata2, newdata3, newdata4, newdata5))\
colnames(a)[1] <- "bins"\
a <- head(a, -1) \
a <- a[order(a$bins),] \
\
#prep pathway data\
\
same <- path10less[path10less$compare=="same",]\
diff <- path10less[path10less$compare=="different",]\
\
df <- data.frame(x = cut(same$dist, breaks=rng, labels=labels,include.lowest = TRUE), y= same$cor, z=same$path1)\
df1 <- data.frame(x = cut(diff$dist, breaks=rng, labels=labels,include.lowest = TRUE), y= diff$cor)\
\
\
\
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0
\cf2 ################# generate plotting data for distances 10 mb or more ##############\
\
#read in pathway allocations\
\
path10more <- read.table("all.merged.path.10mbpmore.rho2", header=T)\
\
#read in cor and distance data\
\
decay10more <- read.table("all.merged.dup.NA.sorted.10Mbmore.rho2.txt", header=T)\
\
#prep background data\
rng0 <- seq(10000000,280000000, by=5400000)\
labels0 <- rng0[1:49]\
labels10 <- labels0+1 \
\
newdata0<-ddply(decay10more,.(cut(value.x,breaks=rng0,labels=labels0)),summarise,q05=quantile(value.y, 0.05))\
newdata10<-ddply(decay10more,.(cut(value.x,breaks=rng0,labels=labels0)),summarise,q25=quantile(value.y, 0.25))\
newdata20<-ddply(decay10more,.(cut(value.x,breaks=rng0,labels=labels0)),summarise,q5=quantile(value.y, 0.50))\
newdata30<-ddply(decay10more,.(cut(value.x,breaks=rng0,labels=labels0)),summarise,q75=quantile(value.y, 0.75))\
newdata40<-ddply(decay10more,.(cut(value.x,breaks=rng0,labels=labels0)),summarise,q95=quantile(value.y, 0.95))\
newdata50<-ddply(decay10more,.(cut(value.x,breaks=rng0,labels=labels0)),summarise,mean=mean(value.y))\
a0 <- Reduce(function(x, y) merge(x, y, all=TRUE), list(newdata0, newdata10, newdata20, newdata30, newdata40, newdata50))\
colnames(a0)[1] <- "bins"\
a0 <- head(a0, -1) \
a0 <- a0[order(a0$bins),] \
\
#prep pathway data\
\
same0 <- path10more[path10more$compare=="same",]\
diff0 <- path10more[path10more$compare=="different",]\
\
df0 <- data.frame(x = cut(same0$dist, breaks=rng0, labels=labels0,include.lowest = TRUE), y= same0$cor,z=same0$path1)\
df10 <- data.frame(x = cut(diff0$dist, breaks=rng0, labels=labels0,include.lowest = TRUE), y= diff0$cor)\
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\tx14960\pardirnatural\partightenfactor0
\cf2 \
\
\
\
########################################### merge 10mb or less with 10mb or more #########\
\
aboth <- rbind(a,a0)\
sameboth <- rbind(same,same0)\
diffboth <- rbind(diff,diff0)\
labels1both <- append(labels1,labels10)\
\
\
#################################    plot it all ####\
\
\
\
\
\
pdf("10mb.col.rev.final1.pdf",h=7,w=10)\
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0
\cf2 ggplot(aboth, aes(x = as.numeric(log10(labels1both)))) +\
 geom_point(data=diffboth, aes(x=as.numeric(log10(diffboth$dist)), y=diffboth$cor, color=factor(diffboth$compare), fill= factor(diffboth$compare)), size=2, shape=21, alpha = 0.75)+\
 geom_point(data=sameboth, aes(x=as.numeric(log10(sameboth$dist)), y=sameboth$cor, color=factor(sameboth$path1), fill= factor(sameboth$path1)), shape=21, size=2, alpha = 0.75)+\
 \expnd0\expndtw0\kerning0
\CocoaLigature1 scale_color_manual(values = c('noscapine'=\kerning1\expnd0\expndtw0 \CocoaLigature0 "black",'reticuline'=rgb(0.5,0,0), 'sanguinarine'="red", 'thebaine'="purple", 'unknown'=rgb(0,0.5,1), 'different'="black"\expnd0\expndtw0\kerning0
\CocoaLigature1 ))+\kerning1\expnd0\expndtw0 \CocoaLigature0 \
 \expnd0\expndtw0\kerning0
\CocoaLigature1 scale_fill_manual(values = c('noscapine'=\kerning1\expnd0\expndtw0 \CocoaLigature0 "black",'reticuline'=rgb(0.5,0,0), 'sanguinarine'="red", 'thebaine'="purple", 'unknown'=rgb(0,0.5,1), 'different'="white"\expnd0\expndtw0\kerning0
\CocoaLigature1 ))+\kerning1\expnd0\expndtw0 \CocoaLigature0 \
 geom_ribbon(aes(x = as.numeric(log10(labels1both)),ymax=q95, ymin=q75), alpha  = .3)+\
 geom_ribbon(aes(x = as.numeric(log10(labels1both)),ymax=q75, ymin=q25), alpha  = .2)+\
 geom_ribbon(aes(x = as.numeric(log10(labels1both)),ymax=q25, ymin=q05), alpha  = .3)+\
 geom_line(aes(x = as.numeric(log10(labels1both)),y=aboth$mean),size=0.25)+\
 scale_x_continuous(name="Intergenic distance log10[D]", expand = c(0.005, 0)\expnd0\expndtw0\kerning0
\CocoaLigature1 )\kerning1\expnd0\expndtw0 \CocoaLigature0 +\
 scale_y_continuous(name=expression("Pairwise gene expression correlation" (\expnd0\expndtw0\kerning0
\CocoaLigature1 rho\kerning1\expnd0\expndtw0 \CocoaLigature0 ^2)), expand = c(0, 0.07), breaks=c(0,0.25,0.5,0.75,1))+\
 geom_vline(xintercept=log10(100000), linetype="dashed", color = "black", size=0.25)+\
 geom_vline(xintercept=log10(200000), linetype="dashed", color = "black", size=0.25)+\
 geom_vline(xintercept=log10(500000), linetype="dashed", color = "black", size=0.25)+\
 geom_vline(xintercept=log10(1000000), linetype="dashed", color = "black", size=0.25)+\
 geom_vline(xintercept=log10(10000000), linetype="dashed", color = "black", size=0.25)+\
 annotate("text", label="0.1Mb", x=log10(80000), y=1.1, angle=90, size=5)+\
 annotate("text", label="0.2Mb", x=log10(160000), y=1.1, angle=90, size=5)+\
 annotate("text", label="0.5Mb", x=log10(400000), y=1.1, angle=90, size=5)+\
 annotate("text", label="1.0Mb", x=log10(800000), y=1.1, angle=90, size=5)+\
 annotate("text", label="10.0Mb", x=log10(8000000), y=1.1, angle=90, size=5)+\
 #annotate("text", label="<10Mb", x=log10(4500000), y=0.95, size=3)+\
 #annotate("text", label=">10Mb", x=log10(130000000), y=0.95, size=3)+\
 annotate("text", label="", x=log10(130000000), y=0.9, size=3)+\
\
 #geom_rect(xmin = log10(10000000), xmax = log10(297000000), ymin = -0.85, ymax = 1.09, fill=NA, color="grey50", size=0.25)+\
\
 theme_bw()+\
 theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(), \expnd0\expndtw0\kerning0
\CocoaLigature1 legend.title =\kerning1\expnd0\expndtw0 \CocoaLigature0 element_blank(), panel.background = element_blank(), \cf4 \cb5 \expnd0\expndtw0\kerning0
\CocoaLigature1 text = element_text(size=\cf6 22)\cf2 \cb3 \kerning1\expnd0\expndtw0 \CocoaLigature0 )\
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\tx14960\pardirnatural\partightenfactor0
\cf2 dev.off()\
\
}