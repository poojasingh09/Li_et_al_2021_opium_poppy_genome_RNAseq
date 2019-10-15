# PSingh 2019 #
## Opium poppy paper (Li et al)) ##


###############

library(WGCNA)
path <- read.table("SPPP3_Genes_expression_counts_108BIAgenes.txt", header=T, row.names=1)
d <- dist(path[1:13], method = "euclidean") 
fit <- hclust(d, method="ward.D") 
#lab <- labels2colors(path[,c(19,16)], naColor="white", commonColorCode=TRUE)
col= c("orange","pink","cyan","black", "white", rgb(0.75,0.75,0.25), rgb(0.75,0.75,0.25), "white", rgb(0,0.5,0.5), rgb(0.5,0,0), rgb(0.5,0,0), "red", "purple", "deepskyblue2")
lab1 <- labels2colors(path$pathway_new, colorSeq=col)

clust <- path$clust0.5
clust[clust < 0] <- NA

col1 = c("yellow", "green", "lightgreen", "lightblue1", "brown", "cyan4", "burlywood", "antiquewhite3", "deeppink", "aquamarine4", "chocolate", "grey", "aquamarine", "chocolate2", "darkseagreen3", "chartreuse3", "bisque3", "azure2")

lab2 <- labels2colors(clust, colorSeq=col1, naColor="white")

lab <- cbind(lab1, lab2)


colnames(lab) <- c("pathway","clust0.5")
pdf("dendrowithcolors.pdf", w=12,h=6)
plotDendroAndColors(fit, lab, cex.dendroLabels = 0.6,colorHeight=0.1,autoColorHeight=F)
dev.off()
