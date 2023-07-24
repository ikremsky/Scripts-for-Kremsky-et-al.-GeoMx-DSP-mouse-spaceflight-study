.libPaths(c("/home/samir/R/x86_64-pc-linux-gnu-library/4.0/", .libPaths()))
library("BSDA")
args=commandArgs(trailingOnly=T)

x=read.delim(args[1],header=T)
name=as.character(args[2])

Gene=x[,1]
Log2=x[,2]
head(Gene)
head(Log2)

png(paste(name,"_Log2_barplot.png",sep=""), height = 2000, width = 2000, res = 300)
barplot(Log2, names=Gene, col="red", ylab="Log2", cex.names=1.3)
dev.off()
