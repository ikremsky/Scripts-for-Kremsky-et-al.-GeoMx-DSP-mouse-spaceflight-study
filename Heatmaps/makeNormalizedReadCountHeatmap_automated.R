.libPaths(c("/home/samir/R/x86_64-pc-linux-gnu-library/4.0/", .libPaths()))
# Load required packages
library("ComplexHeatmap")
library("circlize")

# Load data
args=as.character(commandArgs(trailingOnly=T))
celltype=as.character(args[1:length(args)])
NRC_mat=read.delim(file="/data2/samir/Vmao_Nanopore/heatmaps/GeneExpressionTable.sorted.mod", header=T,row.names=1)
Log2_mat=read.delim(file="/data2/samir/Vmao_Nanopore/heatmaps/Hclust_Log2Table.txt",header=T,row.names=1)

# Match row order
NRC_mat=NRC_mat[rownames(Log2_mat),]

# Identify the columns to exclude
exclude_cols <- grep("TX", colnames(NRC_mat))
new_NRC_mat <- NRC_mat[, -exclude_cols]

# Select the appropriate column from the Normalized Read Count (NRC) matrix
for (sample in celltype)
{
	# Normalized ReadCount Matrix is created below
        assign(paste("NRC_", sample, sep=""), new_NRC_mat[,grep(sample, colnames(new_NRC_mat))])

        assign(paste("NRC_", sample, "_GC", sep=""), get(paste("NRC_", sample,sep=""))[,grep("GC", colnames(get(paste("NRC_", sample,sep=""))))])
        assign(paste("NRC_", sample, "_FLT", sep=""), get(paste("NRC_", sample,sep=""))[,grep("FLT", colnames(get(paste("NRC_", sample,sep=""))))])

	assign(paste("NRC_", sample, "_GC_mean", sep=""), apply(get(paste("NRC_", sample, "_GC", sep="")), 1, mean))
	assign(paste("NRC_", sample, "_FLT_mean", sep=""), apply(get(paste("NRC_", sample, "_FLT", sep="")), 1, mean))

	assign(paste("NRC_", sample, "_Final", sep=""), cbind(get(paste("NRC_", sample, "_GC_mean", sep="")), get(paste("NRC_", sample, "_FLT_mean", sep=""))))
	new_var_NRC <- paste("NRC_", sample, "_Final_named", sep="")
	new_val_NRC <- setNames(data.frame(get(paste("NRC_", sample, "_Final", sep=""))), c(paste(sample, "_GC", sep=""), paste(sample, "_FLT", sep="")))
	new_var_NRC <- log(new_val_NRC,2)

	# Log2 Matrix is created below
	assign(paste("Log2_", sample, sep=""), Log2_mat[,grep(sample, colnames(Log2_mat))])
	new_var_Log2 <- paste("Log2_", sample, "_named", sep="")
        new_val_Log2 <- list(setNames(data.frame(get(paste("Log2_", sample, sep=""))), paste(sample, "_Log2", sep="")))
        new_var_Log2 <- new_val_Log2
	
	# Heatmaps are made below
	assign(paste("NRC_", sample, sep=""), Heatmap(as.matrix(new_var_NRC), name = paste("NRC_", sample, sep=""), col = colorRamp2(c(2, 6), c("white", "black")), column_order = c(paste(sample, "_FLT", sep=""), paste(sample, "_GC", sep="")), show_row_names = FALSE, cluster_rows = FALSE, show_heatmap_legend = FALSE))
	assign(paste("Log2_", sample, sep=""), Heatmap(data.frame(new_var_Log2), name = paste("Log2_", sample, sep=""), col = colorRamp2(c(-.75, -.25, .25, .75), c("blue", "white", "white", "red")), show_row_names = FALSE, cluster_rows = FALSE, show_heatmap_legend = TRUE))
}
CA=get(paste("Log2_", args[1], sep="")) + get(paste("NRC_", args[1], sep=""))
DG=get(paste("Log2_", args[2], sep="")) + get(paste("NRC_", args[2], sep=""))
FCortex=get(paste("Log2_", args[3], sep="")) + get(paste("NRC_", args[3], sep=""))
RCortex=get(paste("Log2_", args[4], sep="")) + get(paste("NRC_", args[4], sep=""))

png(paste("heatmap_automated_Combined.png", sep=""), height = 2000, width = 2000, res = 300)
draw(CA + DG + FCortex + RCortex)
dev.off()
