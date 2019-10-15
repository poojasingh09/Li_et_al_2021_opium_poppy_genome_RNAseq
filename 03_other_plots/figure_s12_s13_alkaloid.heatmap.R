### PSingh 2019####
### Opium poppy paper (Li et al) ####


### first select only known alkaloids on the terminal

grep -v _unknown SPPP4_OrbitrapAlkaloid_Data.txt > s12.txt


#### this part is run in R ####
path <- read.table("s12.txt", header=T, row.names=1)
pdf("s12.pdf", h=40, w=20)
pheatmap(as.matrix(log1p(path)),cexRow=0.6, cexCol=0.8, cellwidth=15, cellheight=30, fontsize_row=8, cluster_cols=F,cluster_rows=F,na_col = "white",border_color="black")
dev.off()


path <- read.table("SPPP4_OrbitrapAlkaloid_Data.txt", header=T, row.names=1)
pdf("s12.pdf", h=40, w=20)
pheatmap(as.matrix(log1p(path)),cexRow=0.6, cexCol=0.8, cellwidth=15, cellheight=30, fontsize_row=8, cluster_cols=F,cluster_rows=F,na_col = "white",border_color="black")
dev.off()

