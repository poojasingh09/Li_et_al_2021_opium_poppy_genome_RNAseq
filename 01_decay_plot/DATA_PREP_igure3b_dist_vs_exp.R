####plot gene co-expression decay ####
####psingh march 2019 ####
####poppy project yeaman lab (Li et al) ####


#0 First split the gff and gene expression matrix into separate scaffolds. and run this script for each scaffold. then merge the output files.

#1 read gtf and calculate distance between gene TSSs as a matrix3, for genes on different chromosomes, set default value of alpha

setwd(".")

library(data.table, lib.loc="/data/home/psingh/tools/R")
library(spaa,lib.loc="/data/home/psingh/tools/R")


genes <- fread("SPPP1A_Psomniferum_HIC_assembly_Li_et_al_s0.gff", header=F, sep="\t")
new_names <- c("chr","source","type","start","end","score","strand","phase","attributes")
setnames(genes, new_names)
genes <- genes[type == "mRNA"]


genes$att <- strsplit(genes$attributes, ";")
genes$geneID <- sapply(genes$att, "[", 3)
genes$geneID <- gsub('Parent=', '', genes$geneID)

for (i in genes){
	if (genes$strand == "+"){
		genes$gene_TSS <- genes$start
		} 
	else {
		genes$gene_TSS <- genes$end
	}
}

genes1 <- data.frame(genes$geneID, genes$chr, genes$gene_TSS)
colnames(genes1) <- c("geneID", "chr", "TSS")
#rownames(genes1) <- genes1$geneID

splitgenes1 <- split(genes1, genes1$chr)
for (split in splitgenes1){
	mdist <- dist(split$TSS)
	mdist <- as.matrix(mdist)
	rownames(mdist) <- split$geneID
	colnames(mdist) <- split$geneID
}

splitgenes1 <- split(genes1, genes1$chr)
for (split in splitgenes1){
	mdist <- dist(split$TSS)
	mdist <- as.matrix(mdist)
	rownames(mdist) <- split$geneID
	colnames(mdist) <- split$geneID
}

mdist <- dist(genes1$TSS)
mdist <- as.matrix(mdist)
rownames(mdist) <- genes1$geneID
colnames(mdist) <- genes1$geneID
distlist <- dist2list(as.dist(mdist))
distmat <- list2dist(mdist)

write.table(distlist, distlist.all.txt, sep="\t", quote=FALSE)
write.table(distmat, distmat.all.txt, sep="\t", quote=FALSE)


#2 read gene expression matrix and get spearman cor between genes across all samples and generate matrix1

exp <- read.table("SPPP2_Gene_expression_counts_all_samples_filt_s0.txt", header=T, sep="\t",  row.names=1)
cor <- cor(t(exp), method="spearman")
corlist <- dist2list(as.dist(cor))

write.table(corlist, corlist.all.txt, sep="\t", quote=FALSE)
write.table(cor, cormat.all.txt, sep="\t", quote=FALSE)

merge(distlist,corlist, by=c("col", "row"))
write.table(merge, mergeall.txt, sep="\t", quote=FALSE)
save.image(file = "decay.all.test.RData")


### the merged files are provided (figure_3b_inputfile_all.merged.dup.NA.sorted.10Mb.rho2; figure_3b_inputfile_all.merged.dup.NA.sorted.10Mbmore.rho2)

