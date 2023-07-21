.libPaths()
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


#M=as.matrix(cbind(as.numeric(a[,4]),as.numeric(b[,4]),as.numeric(c[,4]),as.numeric(d[,4])))
#colnames(M)=c("FCortex","Cortex","CA","DG")
head(M)
tail(M)
length(M[,1])
#index=((a[,8]>1.3 & abs(a[,4])>1) | (b[,8]>1.3 & abs(b[,4])>1) | (c[,8]>1.3 & abs(c[,4])>1) | (d[,8]>1.3 & abs(d[,4])>1))
#M=na.omit(-M[index,])
ordr=hclust(dist(M))$order

myCol <- c(colorRampPalette(c("blue", "white", "red"))(303))
myBreaks <- c(-3,seq(-.75,-.1,length=150), 0, seq(.1,.75,length=150),3)
#png(paste(Compare, "_heatmap_Vmao.png", sep=""), height = 2000, width = 2000, res = 300)
#pheatmap(M[ordr,], cluster_rows=F, cluster_cols=F, scale="none", legend=T, fontsize_col=20,color=myCol,breaks=myBreaks, show_rownames = F)
#dev.off()

write.table(M[ordr,], file="Hclust_Log2Table.txt", sep="\t", quote=FALSE)
