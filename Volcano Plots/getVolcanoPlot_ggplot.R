library(ggrepel)
library(dplyr)
library(ggplot2)
args=commandArgs(trailingOnly=T)
i=1
tmp=read.delim(as.character(args[i]), sep="\t", header=T); i=i+1
name=as.character(args[i]); i=i+1
y=read.delim(as.character(args[i]),header=F); i=i+1

de <- tmp[complete.cases(tmp), ]

png(paste(name, "_ggplot_volcanoplot", ".png", sep=""), height = 2000, width = 2000, res=300)

de$diffexpressed <- "NO"
signif=de$Adjusted_pvalue < .05 & abs(de$log2FoldChange) > .585

de$diffexpressed[signif] <- "Significant"
de$diffexpressed[!signif] <- "NS"

length(which(signif))

top5degs <- print(y[,1])
de$delabel <- ifelse(de$gene_symbol %in% top5degs, de$gene_symbol, NA)

ggplot(data=de, aes(x=log2FoldChange, y=NeglogPadjValue, col=diffexpressed, label=delabel)) + geom_point(size = 3) + geom_text_repel(data=filter(de), aes(label = delabel), size = 5, nudge_y = -0.30, nudge_x = -0.5, max.overlaps = Inf, show.legend = F) + scale_color_manual(values=c("black", "red")) + labs(color = 'DiffExpressed', x = "Log2 Fold Change", y = "-Log10 padj-value") + ggtitle(paste(name, "_GCvsFLT-SAL", sep="")) + theme_classic() + coord_cartesian(ylim = c(0, 2.2), xlim = c(-2, 2), clip = "off")
dev.off()
