.libPaths(c("/home/samir/R/x86_64-pc-linux-gnu-library/4.0/", .libPaths()))
library("ComplexHeatmap")
library("circlize")
library("pheatmap")
library("fastcluster")
args=as.character(commandArgs(trailingOnly=T))
M=read.delim(args[1],header=T,row.names=1)
Compare=as.character(args[2])
#M=M[,-1]

#a=read.delim(args[1])
#b=read.delim(args[2])
#c=read.delim(args[3])
#d=read.delim(args[4])

#a[which(a[,4] > 3),4]=3
#a[which(a[,4] < -3),4]=-3
#b[which(b[,4] > 3),4]=3
#b[which(b[,4] < -3),4]=-3
#c[which(c[,4] > 3),4]=3
#c[which(c[,4] < -3),4]=-3
#d[which(d[,4] > 3),4]=3
#d[which(d[,4] < -3),4]=-3


length(M[,1])
#M=na.omit(-M[index,])
#N=M[,-grep("SALvTX.GC", colnames(M))]
N=M
head(N)
ordr=hclust(dist(N))$order

H=Heatmap(M[ordr,], name = paste("Log2FC", sep=""), col = colorRamp2(c(-.75, -.1, .1, .75), c("blue", "white", "white", "red")), show_row_names = T, cluster_rows = FALSE, show_heatmap_legend = TRUE, cluster_columns = FALSE, row_names_gp = gpar(fontsize = 20))

png(paste(Compare, "_heatmap_Vmao.png", sep=""), height = 2000, width = 2000, res = 300)
draw(H)
dev.off()
