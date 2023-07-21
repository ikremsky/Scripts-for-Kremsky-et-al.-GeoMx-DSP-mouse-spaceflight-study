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

CA_FLT_TX=y[,grep("CA_FLT_TX", names)]
CA_GC_SAL=y[,grep("CA_GC_SAL", names)]
CA_GC_TX=y[,grep("CA_GC_TX", names)]
CA_FLT_SAL=y[,grep("CA_FLT_SAL", names)]

DG_FLT_TX=y[,grep("DG_FLT_TX", names)]
DG_GC_SAL=y[,grep("DG_GC_SAL", names)]
DG_GC_TX=y[,grep("DG_GC_TX", names)]
DG_FLT_SAL=y[,grep("DG_FLT_SAL", names)]

FCortex_FLT_TX=y[,grep("FCortex_FLT_TX", names)]
FCortex_GC_SAL=y[,grep("FCortex_GC_SAL", names)]
FCortex_GC_TX=y[,grep("FCortex_GC_TX", names)]
FCortex_FLT_SAL=y[,grep("FCortex_FLT_SAL", names)]

RCortex_FLT_TX=y[,grep("RCortex_FLT_TX", names)]
RCortex_GC_SAL=y[,grep("RCortex_GC_SAL", names)]
RCortex_GC_TX=y[,grep("RCortex_GC_TX", names)]
RCortex_FLT_SAL=y[,grep("RCortex_FLT_SAL", names)]

y=t(y)
M=as.matrix(cbind(CA_FLT_TX, CA_GC_SAL, CA_GC_TX, CA_FLT_SAL, DG_FLT_TX, DG_GC_SAL, DG_GC_TX, DG_FLT_SAL, FCortex_FLT_TX, FCortex_GC_SAL, FCortex_GC_TX, FCortex_FLT_SAL, RCortex_FLT_TX, RCortex_GC_SAL, RCortex_GC_TX, RCortex_FLT_SAL))
print(names)
#head(M)
#head(x)
P=NULL
ses=factor(c(rep(1,3), rep(2,3),rep(3,3), rep(4,3), rep(5,3), rep(6,3), rep(7,3), rep(8,3), rep(9,3), rep(10,3), rep(11,3), rep(12,3), rep(13,3), rep(14,3), rep(15,3), rep(16,3)))
levels(ses)=c("CA_FLT_TX", "CA_GC_SAL", "CA_GC_TX", "CA_FLT_SAL", "DG_FLT_TX", "DG_GC_SAL", "DG_GC_TX", "DG_FLT_SAL", "FCT_FLT_TX", "FCT_GC_SAL", "FCT_GC_TX", "FCT_FLT_SAL", "CT_FLT_TX", "CT_GC_SAL", "CT_GC_TX", "CT_FLT_SAL")
for(i in 1:length(x[,1]))
{
        a=aov(M[i,] ~ ses)
        P=c(P, unlist(summary(a))["Pr(>F)1"])
}
q=qvalue(P)$qvalues
#q=qvalue(P, lambda = seq(min(P), max(P), 0.025), smooth.log.pi0="TRUE")$qvalues

signif=which(q<.05)
#head(signif)
min(q)
CA_FLT_TX=CA_FLT_TX[signif,]
CA_GC_SAL=CA_GC_SAL[signif,]
CA_GC_TX=CA_GC_TX[signif,]
CA_FLT_SAL=CA_FLT_SAL[signif,]

DG_FLT_TX=DG_FLT_TX[signif,]
DG_GC_SAL=DG_GC_SAL[signif,]
DG_GC_TX=DG_GC_TX[signif,]
DG_FLT_SAL=DG_FLT_SAL[signif,]

FCortex_FLT_TX=FCortex_FLT_TX[signif,]
FCortex_GC_SAL=FCortex_GC_SAL[signif,]
FCortex_GC_TX=FCortex_GC_TX[signif,]
FCortex_FLT_SAL=FCortex_FLT_SAL[signif,]

RCortex_FLT_TX=RCortex_FLT_TX[signif,]
RCortex_GC_SAL=RCortex_GC_SAL[signif,]
RCortex_GC_TX=RCortex_GC_TX[signif,]
RCortex_FLT_SAL=RCortex_FLT_SAL[signif,]

length(signif)
#head(y)

x = cbind(CA_FLT_TX,CA_GC_SAL,CA_GC_TX,CA_FLT_SAL,DG_FLT_TX,DG_GC_SAL,DG_GC_TX,DG_FLT_SAL,FCortex_FLT_TX,FCortex_GC_SAL,FCortex_GC_TX,FCortex_FLT_SAL,RCortex_FLT_TX,RCortex_GC_SAL,RCortex_GC_TX,RCortex_FLT_SAL)
colnames(x)=c(rep("CA.FLT.TX", 3), rep("CA.GC.SAL", 3), rep("CA.GC.TX", 3), rep("CA.FLT.SAL", 3), rep("DG.FLT.TX", 3), rep("DG.GC.SAL", 3), rep("DG.GC.TX", 3), rep("DG.FLT.SAL", 3), rep("FCT.FLT.TX", 3), rep("FCT.GC.SAL", 3), rep("FCT.GC.TX", 3), rep("FCT.FLT.SAL", 3), rep("CT.FLT.TX", 3), rep("CT.GC.SAL", 3), rep("CT.GC.TX", 3), rep("CT.FLT.SAL", 3))
y=t(x)

png(paste("PCA_AllSamples_AllReps_Vmao", ".png", sep=""), height = 2000, width = 2000, res = 300)
P=PCA(y, graph = F, scale.unit = T)
p=plot(P,choix="ind",habillage ="none", col.ind=rep("white",20), graph.type="classic",label="none")
legend("top", ncol=2, pch=c(rep(16,4),0,16,15,1), legend= c("CA","DG","FCT","CT","FLT-BuOE", "GC-SAL", "GC-BuOE", "FLT-SAL"), col= c("red","cyan3","gold4","purple",rep("black",4)), bty="o")
points(P$ind$coord[,1],P$ind$coord[,2],col=c(rep("red",12),rep("cyan3",12),rep("gold4",12),rep("purple",12)), pch=rep(c(rep(0, 3),rep(16, 3),rep(15, 3),rep(1, 3)), 4))
dev.off()

CA_FLT_TX=apply(CA_FLT_TX, 1, mean)
CA_GC_SAL=apply(CA_GC_SAL, 1, mean)
CA_GC_TX=apply(CA_GC_TX, 1, mean)
CA_FLT_SAL=apply(CA_FLT_SAL, 1, mean)

DG_FLT_TX=apply(DG_FLT_TX, 1, mean)
DG_GC_SAL=apply(DG_GC_SAL, 1, mean)
DG_GC_TX=apply(DG_GC_TX, 1, mean)
DG_FLT_SAL=apply(DG_FLT_SAL, 1, mean)

FCortex_FLT_TX=apply(FCortex_FLT_TX, 1, mean)
FCortex_GC_SAL=apply(FCortex_GC_SAL, 1, mean)
FCortex_GC_TX=apply(FCortex_GC_TX, 1, mean)
FCortex_FLT_SAL=apply(FCortex_FLT_SAL, 1, mean)

RCortex_FLT_TX=apply(RCortex_FLT_TX, 1, mean)
RCortex_GC_SAL=apply(RCortex_GC_SAL, 1, mean)
RCortex_GC_TX=apply(RCortex_GC_TX, 1, mean)
RCortex_FLT_SAL=apply(RCortex_FLT_SAL, 1, mean)

x = cbind(CA_FLT_TX,CA_GC_SAL,CA_GC_TX,CA_FLT_SAL,DG_FLT_TX,DG_GC_SAL,DG_GC_TX,DG_FLT_SAL,FCortex_FLT_TX,FCortex_GC_SAL,FCortex_GC_TX,FCortex_FLT_SAL,RCortex_FLT_TX,RCortex_GC_SAL,RCortex_GC_TX,RCortex_FLT_SAL)
colnames(x)=c("CA.FLT.TX", "CA.GC.SAL", "CA.GC.TX", "CA.FLT.SAL", "DG.FLT.TX", "DG.GC.SAL", "DG.GC.TX", "DG.FLT.SAL", "FCT.FLT.TX", "FCT.GC.SAL", "FCT.GC.TX", "FCT.FLT.SAL", "CT.FLT.TX", "CT.GC.SAL", "CT.GC.TX", "CT.FLT.SAL")
y=t(x)

png(paste("PCA_AllSamples_Vmao", ".png", sep=""), height = 2000, width = 2000, res = 300)
P=PCA(y, graph = F, scale.unit = T)
#p=plot(P,choix="ind",habillage="ind",graph.type="classic")
p=plot(P,choix="ind",habillage ="none", col.ind=rep("white",20), graph.type="classic",label="none")
points(P$ind$coord[,1],P$ind$coord[,2],col=c(rep("red",4),rep("cyan3",4),rep("gold4",4),rep("purple",4)), pch=rep(c(0,16,15,1), 4), cex=2)
legend("top", ncol=2, pch=c(rep(16,4),0,16,15,1), legend= c("CA","DG","FCT","CT","FLT-BuOE", "GC-SAL", "GC-BuOE", "FLT-SAL"), col= c("red","cyan3","gold4","purple",rep("black",4)), bty="o", pt.cex=2)
dev.off()
