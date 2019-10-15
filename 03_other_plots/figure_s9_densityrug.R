
### PSingh 2019 ####
## Opium poppy paper (Li et al) ####


############ known BIA pathways genes ####

d <- read.table("figure_s9_inputfile_interscaff_subset_corlist.diff.txt", header=T)
d1<- density(d$cor, adjust=2)
pdf("interscaff_density_rug.pdf")

the_cols <-c("cyan","purple","black","white",rgb(0.5,0,0),"orange","red",rgb(0.75,0.75,0.25),rgb(0.75,0.75,0.25),"white")

leg <- c("morphine","thebaine","noscapine","noscapine_sanguinarine","reticuline","dopamine_4HPAA","sanguinarine","papavarine","papavarine_derivatives","papavarine_sanguinarine")  

plot(d1)
 
rug(d[d$other=="morphine",]$cor, col=the_cols[1], lwd=1.5)
rug(d[d$other=="thebaine",]$cor, col=the_cols[2], lwd=1.5)
rug(d[d$other=="noscapine",]$cor, col=the_cols[3], lwd=1.5)
rug(d[d$other=="noscapine_sanguinarine",]$cor, col=the_cols[4], lwd=1.5)
rug(d[d$other=="reticuline",]$cor, col=the_cols[5], lwd=1.5)
rug(d[d$other=="dopamine_4HPAA",]$cor, col=the_cols[6], lwd=1.5)
rug(d[d$other=="sanguinarine",]$cor, col=the_cols[7], lwd=1.5)
#rug(d[d$other=="unknown",]$cor, col=the_cols[8], lwd=1.5)
rug(d[d$other=="papavarine",]$cor, col=the_cols[8], lwd=1.5)
rug(d[d$other=="papavarine_derivatives",]$cor, col=the_cols[9], lwd=1.5)
rug(d[d$other=="papavarine_sanguinarine",]$cor, col=the_cols[10], lwd=1.5)
#rug(d[d$other=="other",]$cor, col="grey50", lwd=1.5)

legend("topright", legend=leg, fill = the_cols, text.col="black", cex=0.75)
 
dev.off()

####known ##### rho2 ####


d1<- density((d$cor)^2, adjust=2)
pdf("interscaff_density_rug_rho2.pdf")


the_cols <-c("cyan","purple","black",rgb(0.5,0,0),"orange","red",rgb(0.75,0.75,0.25))

leg <- c("morphine","thebaine","noscapine","reticuline","dopamine_4HPAA","sanguinarine","papavarine")  


plot(d1, xlab=expression("Pairwise gene expression correlation" (rho^2)), main=NA)
rug((d[d$other=="morphine",]$cor)^2, col=the_cols[1], lwd=1.5)
rug((d[d$other=="thebaine",]$cor)^2, col=the_cols[2], lwd=1.5)
rug((d[d$other=="noscapine",]$cor)^2, col=the_cols[3], lwd=1.5)
#rug((d[d$other=="noscapine_sanguinarine",]$cor)^2, col=the_cols[4], lwd=1.5)
rug((d[d$other=="reticuline",]$cor)^2, col=the_cols[4], lwd=1.5)
rug((d[d$other=="dopamine_4HPAA",]$cor)^2, col=the_cols[5], lwd=1.5)
rug((d[d$other=="sanguinarine",]$cor)^2, col=the_cols[6], lwd=1.5)
rug((d[d$other=="papavarine",]$cor)^2, col=the_cols[7], lwd=1.5)
rug((d[d$other=="papavarine_derivatives",]$cor)^2, col=the_cols[7], lwd=1.5)
#rug((d[d$other=="papavarine_sanguinarine",]$cor)^2, col=the_cols[10], lwd=1.5)


q1 <- quantile((d$cor)^2, 0.95)
q2 <- quantile((d$cor)^2, 0.75)

abline(v=q1, lty=2)
abline(v=q2, lty=2)

text(x=0.48, y=6,labels="95%", cex=0.5)
text(x=0.18, y=6,labels="75%", cex=0.5)

legend("topright", legend=leg, fill = the_cols, text.col="black", cex=0.75)
 
dev.off()




##### unknown pathway genes######


d1<- density((d$cor)^2, adjust=2)
pdf("interscaff_density_rug_rho2_unknown.pdf")


the_cols <-c("blue")

leg <- c("unknown")  


plot(d1, xlab=expression("Pairwise gene expression correlation" (rho^2)), main=NA)
rug((d[d$other=="unknown",]$cor)^2, col=the_cols[1], lwd=1.5)



q1 <- quantile((d$cor)^2, 0.95)
q2 <- quantile((d$cor)^2, 0.75)

abline(v=q1, lty=2)
abline(v=q2, lty=2)

text(x=0.48, y=6,labels="95%", cex=0.5)
text(x=0.18, y=6,labels="75%", cex=0.5)

legend("topright", legend=leg, fill = the_cols, text.col="black", cex=0.75)
 
dev.off()


