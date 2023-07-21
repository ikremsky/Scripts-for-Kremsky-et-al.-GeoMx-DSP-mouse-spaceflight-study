.libPaths(c("/home/samir/R/x86_64-pc-linux-gnu-library/4.0/", .libPaths()))
library("FactoMineR")
library(fastcluster)
library(gplots)
library("pheatmap")
library("qvalue")
args=as.character(commandArgs(trailingOnly=T))

x=read.delim("/data2/samir/Vmao_Nanopore/tables/GeneExpressionTable.txt", sep="\t", header=T)
sums=apply(x[,-1], 1, max)
y=x[,-1]
genes=x[,1]
names=colnames(y)
####
samples=as.character(args[1:length(args)])
#head(samples)

M=NULL
se=NULL
i=1
for (sample in samples)
{
	assign(sample, y[,grep(sample, names)])
	
	se=c(se, rep(i,3))
	i=i+1
}

M=get(samples[1])

for (j in 2:length(samples))
{
	M=cbind(M,get(samples[j]))
}

y=t(y)
M=as.matrix(M)
#head(M)

P=NULL

ses=factor(se)
levels(ses)=samples

for (i in 1:length(x[,1]))
{
        a=aov(M[i,] ~ ses)
        P=c(P, unlist(summary(a))["Pr(>F)1"])
}
OPmax=P[order(P, decreasing=T)][1]
OPmin=P[order(P, decreasing=F)][1]
q=qvalue(P)$qvalues
#q=qvalue(P, lambda = seq(OPmin, OPmax, 0.025), smooth.log.pi0="TRUE")$qvalues

signif=which(q<.05)

x=NULL
for (sample in samples)
{
	assign(sample, get(sample)[signif,])
}

x=get(samples[1])
#head(x)
for (j in 2:length(samples))
{
        x=cbind(x,get(samples[j]))
}
head(x)
y=t(x)

png(paste("PCA_Vmoa_All_Samples", ".png", sep=""), height = 2000, width = 2000, res = 300)
P=PCA(y, graph = F, scale.unit = T)

p=plot(P,choix="ind",habillage ="none", col.ind=rep("white",20), graph.type="classic",label="none")
points(P$ind$coord[,1],P$ind$coord[,2],col=c(rep("red",6),rep("cyan3",6),rep("gold4",6),rep("purple",6)), pch=c(rep(16,3),rep(1,3),rep(16,3),rep(1,3),rep(16,3),rep(1,3),rep(16,3),rep(1,3)), cex=2)
legend("top", ncol=2, pch=c(rep(16,4),16,1), legend= c("CA","DG","FCortex","Cortex","GC", "FLT"), col= c("red","cyan3","gold4","purple",rep("black",2)), bty="o", pt.cex=2)
dev.off()
