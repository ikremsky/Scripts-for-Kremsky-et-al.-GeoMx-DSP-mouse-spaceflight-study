args=as.character(commandArgs(trailingOnly=T))

region <- args[1]

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

#Input file of overlaps from 100 random shufflings
shuf=unlist(read.delim(as.character(args[10]), header=F))

#Total DEGenes for each overlap combination

total_downdown <- GCvsFLTdown+SALvsTXdown-over_downdown
total_downup <- GCvsFLTdown+SALvsTXup-over_downup
total_updown <- GCvsFLTup+SALvsTXdown-over_updown
total_upup <- GCvsFLTup+SALvsTXup-over_upup

#Combined values for shuf direction (downdown+upup) and opposite (downup+updown)

noncountshuf <- total_downdown+total_upup-over_downdown-over_upup
noncountopposite <- total_downup+total_updown-over_downup-over_updown
countshuf <- median(shuf)
countopposite <- over_downup+over_updown

png(paste(region, "_barplot_shuffled.png", sep=""), height = 2000, width = 2000, res=300)
#mp=barplot(c(countopposite, countshuf), names.arg=c("opposite Dir", "random"), cex.main=2, cex.axis=1.4, cex.lab=2, cex.names=1.7, ylab="", main=region, ylim=c(0,maxH+11))
#plotPvals2(P$p.value,maxH+5,2,mp[1],mp[2])
hist_info <- hist(shuf, plot=F)
max_index <- which.max(hist_info$counts)
x_max <- hist_info$mids[max_index]
y_max <- max(hist_info$counts)
#shift value for CA and FCT
shift=4
#shift value for CT and DG
shift=5
hist(shuf, xlim=c(0,countopposite+5), ylim=c(0, y_max+7), xlab="# of overlapping DEGs (opposite dir)", col="grey", plot=T, main="")

params <- par()
abline(v=countopposite)
text(x=countopposite, y=y_max+10, "actual", cex=1.25)
text(x=x_max+shift, y=y_max+6, paste("permutation distribution \n (N=", length(shuf), ")", sep=""), col="grey", cex=1.25)
dev.off()


