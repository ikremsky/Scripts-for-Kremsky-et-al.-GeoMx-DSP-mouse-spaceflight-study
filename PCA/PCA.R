library("FactoMineR")
library(fastcluster)
library(gplots)
library("pheatmap")
library("qvalue")
args=as.character(commandArgs(trailingOnly=T))

x=read.delim("/home/ikremsky/VmaoData/GeneExpressionTable.txt", sep="\t", header=T)
sums=apply(x[,-1], 1, max)
y=x[,-1]
genes=x[,1]
names=colnames(y)
####
celltype=as.character(args[1])
samples=as.character(args[2:length(args)])

head(samples)

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
P=NULL

ses=factor(se)
levels(ses)=samples

for (i in 1:length(x[,1]))
{
        a=aov(M[i,] ~ ses)
        P=c(P, unlist(summary(a))["Pr(>F)1"])
}

#q=qvalue(P)$qvalues
q=qvalue(P, lambda = seq(min(P), max(P), 0.025), smooth.log.pi0="TRUE")$qvalues

signif=which(q<.05)
head(signif)
x=NULL
#for (sample in samples)
#{
#	assign(sample, get(sample)[signif,])
#}

x=get(samples[1])

for (j in 2:length(samples))
{
        x=cbind(x,get(samples[j]))
}

#colnames(x)=samples
y=t(x)
#rownames(x)=genes[signif]

png(paste("PCA_", celltype, "_RNAseq.png", sep=""), height = 2000, width = 2000, res = 300)
P=PCA(y, graph = F, scale.unit = F)
#p=plot(P,choix="ind",habillage="ind",graph.type="classic", cex=2.5)
p=plot(P,choix="ind",habillage ="none", col.ind=rep("white",20), graph.type="classic",label="none")
points(P$ind$coord[,1],P$ind$coord[,2],col=c("black","black","black","black"),pch=c(rep(15,3),rep(19,3)), cex=1.5)
legend("top", ncol=2, pch=c(15,19), legend=samples, col=rep("black",3), bty="o", cex=1.5)
dev.off()
