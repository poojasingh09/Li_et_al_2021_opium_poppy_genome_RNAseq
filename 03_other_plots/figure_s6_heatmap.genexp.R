### PSINGH 2019 ####
### OPIUM PAPER (Li et al)


library(pheatmap)


path <- read.table("SPPP3_Genes_expression_counts_108BIAgenes.txt", header=T, row.names=1)

pdf("gene_tpm_all_samples.path.pdf", h=12,w=8)



annrow = data.frame(pathway= as.factor(path$pathway_new))
rownames(annrow) = rownames(path)

ann_colors = list(pathway=c(morphine="cyan", thebaine="purple", noscapine="black", reticuline=rgb(0.5,0,0), reticuline_derivatives=rgb(0.5,0,0), dopamine_4HPAA= "orange", sanguinarine="red", unknown=rgb(0,0.5,1), papaverine=rgb(0.75,0.75,0.25), papaverine_derivatives=rgb(0.75,0.75,0.25), protopine_derivatives=rgb(0,0.5,0.5),magnoflorine="pink", papaverine_sanguinarine="white",  noscapine_sanguinarine="white"))

pheatmap(as.matrix(log1p(path[1:13])),cexRow=0.6, cexCol=0.8, cellwidth=20, cellheight=7, fontsize_row=8, cluster_cols=F,annotation_row = annrow, annotation_colors = ann_colors)
dev.off()
