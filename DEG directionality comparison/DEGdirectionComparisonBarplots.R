args=as.character(commandArgs(trailingOnly=T))

TF <- args[1]

#overlaps, ordered downdown, downup, updown, upup

over_downdown <- length(readLines(args[2]))
over_downup <- length(readLines(args[3]))
over_updown <- length(readLines(args[4]))
over_upup <- length(readLines(args[5]))

#DEGenes per treatment, ordered GCvsFLT down, GCvsFLT up, SALvsTX down, SALvsTX up

GCvsFLTdown <- length(readLines(args[6]))
GCvsFLTup <- length(readLines(args[7]))
SALvsTXdown <- length(readLines(args[8]))
SALvsTXup <- length(readLines(args[9]))

#Total DEGenes for each overlap combination

total_downdown <- GCvsFLTdown+SALvsTXdown-over_downdown
total_downup <- GCvsFLTdown+SALvsTXup-over_downup
total_updown <- GCvsFLTup+SALvsTXdown-over_updown
total_upup <- GCvsFLTup+SALvsTXup-over_upup

#Combined values for same direction (downdown+upup) and opposite (downup+updown)

noncountsame <- total_downdown+total_upup-over_downdown-over_upup
noncountopposite <- total_downup+total_updown-over_downup-over_updown
countsame <- over_downdown+over_upup
countopposite <- over_downup+over_updown

plotPvals2 = function(p, maxH, asteriskDist, x1, x2)
{
	magnify=2
        if(p < .01 & p >= .001)
        {
                text(mean(c(x1,x2)), maxH+asteriskDist, labels="*", cex=magnify)
                arrows(x1,maxH,x2,maxH,code=3,lwd=magnify,angle=90,length=0.05,col="purple")
        }
        if(p < .001 & p >= .00001)
        {
                text(mean(c(x1,x2)), maxH+asteriskDist, labels="**", cex=magnify)
                arrows(x1,maxH,x2,maxH,code=3,lwd=magnify,angle=90,length=0.05,col="purple")
        }
        if(p < .00001 & p >= .0000000001)
        {
                text(mean(c(x1,x2)), maxH+asteriskDist, labels="***", cex=magnify)
                arrows(x1,maxH,x2,maxH,code=3,lwd=magnify,angle=90,length=0.05,col="purple")
        }
        if(p < .0000000001)
        {
                text(mean(c(x1,x2)), maxH+asteriskDist, labels="****", cex=magnify)
                arrows(x1,maxH,x2,maxH,code=3,lwd=magnify,angle=90,length=0.05,col="purple")
        }
}

#percent1 <- 100*countsame/(countsame + noncountsame)
#percent2 <- 100*countopposite/(countopposite + noncountopposite)

M <- as.matrix(rbind(c(countsame, noncountsame), c(countopposite, noncountopposite)))
maxH <- max(countsame, countopposite)+1
head(maxH)

P <- fisher.test(M)
png(paste(TF, "_barplot_TEPIC.png", sep=""), height = 2000, width = 2000, res=300)
mp=barplot(c(countopposite, countsame), names.arg=c("opposite Dir", "same Dir"), cex.main=2, cex.axis=1.4, cex.lab=2, cex.names=1.7, ylab="", main=TF, ylim=c(0,maxH+11))
plotPvals2(P$p.value,maxH+5,2,mp[1],mp[2])
title(ylab="Overlapping DEGs", line=2.5, cex.lab=1.8)
dev.off()


