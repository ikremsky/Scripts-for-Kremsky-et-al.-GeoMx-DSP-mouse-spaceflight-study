library(FactoMineR)
library("fastcluster")
library("gplots")
library("pheatmap")
library("qvalue")

x=read.delim("/data2/samir/Vmao_Nanopore/tables/GeneExpressionTable.txt", sep="\t", header=T)
sums=apply(x[,-1], 1, max)
y=x[,-1]
genes=x[,1]
names=colnames(y)

CA_GC_SAL=y[,grep("CA_GC_SAL", names)]
CA_FLT_SAL=y[,grep("CA_FLT_SAL", names)]

DG_GC_SAL=y[,grep("DG_GC_SAL", names)]
DG_FLT_SAL=y[,grep("DG_FLT_SAL", names)]

FCortex_GC_SAL=y[,grep("FCortex_GC_SAL", names)]
FCortex_FLT_SAL=y[,grep("FCortex_FLT_SAL", names)]

RCortex_GC_SAL=y[,grep("RCortex_GC_SAL", names)]
RCortex_FLT_SAL=y[,grep("RCortex_FLT_SAL", names)]

y=t(y)
M=as.matrix(cbind(CA_GC_SAL, CA_FLT_SAL, DG_GC_SAL, DG_FLT_SAL, FCortex_GC_SAL, FCortex_FLT_SAL, RCortex_GC_SAL, RCortex_FLT_SAL))
#head(M)

P=NULL
ses=factor(c(rep(1,3), rep(2,3),rep(3,3), rep(4,3), rep(5,3), rep(6,3), rep(7,3), rep(8,3)))
levels(ses)=c("CA_GC_SAL", "CA_FLT_SAL", "DG_GC_SAL", "DG_FLT_SAL", "FCortex_GC_SAL", "FCortex_FLT_SAL", "RCortex_GC_SAL", "RCortex_FLT_SAL")
for(i in 1:length(x[,1]))
{
        a=aov(M[i,] ~ ses)
        P=c(P, unlist(summary(a))["Pr(>F)1"])
}
q=qvalue(P)$qvalues
#q=qvalue(P, lambda = seq(min(P), max(P), 0.025), smooth.log.pi0="TRUE")$qvalues

signif=which(q<.05)

CA_GC_SAL=CA_GC_SAL[signif,]
CA_FLT_SAL=CA_FLT_SAL[signif,]

DG_GC_SAL=DG_GC_SAL[signif,]
DG_FLT_SAL=DG_FLT_SAL[signif,]

FCortex_GC_SAL=FCortex_GC_SAL[signif,]
FCortex_FLT_SAL=FCortex_FLT_SAL[signif,]

RCortex_GC_SAL=RCortex_GC_SAL[signif,]
RCortex_FLT_SAL=RCortex_FLT_SAL[signif,]

#length(signif)

CA_GC_SAL=apply(CA_GC_SAL, 1, mean)
CA_FLT_SAL=apply(CA_FLT_SAL, 1, mean)

DG_GC_SAL=apply(DG_GC_SAL, 1, mean)
DG_FLT_SAL=apply(DG_FLT_SAL, 1, mean)

FCortex_GC_SAL=apply(FCortex_GC_SAL, 1, mean)
FCortex_FLT_SAL=apply(FCortex_FLT_SAL, 1, mean)

RCortex_GC_SAL=apply(RCortex_GC_SAL, 1, mean)
RCortex_FLT_SAL=apply(RCortex_FLT_SAL, 1, mean)


x = cbind(CA_GC_SAL, CA_FLT_SAL, DG_GC_SAL, DG_FLT_SAL, FCortex_GC_SAL, FCortex_FLT_SAL, RCortex_GC_SAL, RCortex_FLT_SAL)
colnames(x)=c("CA_GC", "CA_FLT", "DG_GC", "DG_FLT", "FCortex_GC", "FCortex_FLT", "Cortex_GC", "Cortex_FLT")
head(x)
y=t(x)

png(paste("PCA_FLTvsGC-SAL_Samples_Vmao.png", sep=""), height = 2000, width = 2000, res = 300)
P=PCA(y, graph = F, scale.unit = T)
p=plot(P,choix="ind",habillage="ind",graph.type="classic")
dev.off()
