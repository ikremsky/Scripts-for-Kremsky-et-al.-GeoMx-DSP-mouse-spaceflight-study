# Load required packages
library("ComplexHeatmap")
library("circlize")

# Load data
args=as.character(commandArgs(trailingOnly=T))
NRC_mat=read.delim(args[1],header=T,row.names=1)
Log2_mat=read.delim(args[2],header=T,row.names=1)

# Match row order
NRC_mat=NRC_mat[rownames(Log2_mat),]

# Identify the columns to exclude
exclude_cols <- grep("TX", colnames(NRC_mat))

# Create a new matrix without the excluded columns
new_NRC_mat <- NRC_mat[, -exclude_cols]

# Select the appropriate column from the Normalized Read Count (NRC) matrix
NRC_CA <- new_NRC_mat[,grep("CA", colnames(new_NRC_mat))]
NRC_DG <- new_NRC_mat[,grep("DG", colnames(new_NRC_mat))]
NRC_FCortex <- new_NRC_mat[,grep("FCortex", colnames(new_NRC_mat))]
NRC_RCortex <- new_NRC_mat[,grep("RCortex", colnames(new_NRC_mat))]

# Create the inputs for the NRC heatmaps
NRC_CA_GC <- log(apply(NRC_CA[ ,grep("GC", colnames(NRC_CA))], 1, mean),2)
NRC_CA_FLT <- log(apply(NRC_CA[ ,grep("FLT", colnames(NRC_CA))], 1, mean),2)
NRC_CA_Final <- cbind(NRC_CA_GC, NRC_CA_FLT)

NRC_DG_GC <- log(apply(NRC_DG[ ,grep("GC", colnames(NRC_DG))], 1, mean),2)
NRC_DG_FLT <- log(apply(NRC_DG[ ,grep("FLT", colnames(NRC_DG))], 1, mean),2)
NRC_DG_Final <- cbind(NRC_DG_GC, NRC_DG_FLT)

NRC_FCortex_GC <- log(apply(NRC_FCortex[ ,grep("GC", colnames(NRC_FCortex))], 1, mean),2)
NRC_FCortex_FLT <- log(apply(NRC_FCortex[ ,grep("FLT", colnames(NRC_FCortex))], 1, mean),2)
NRC_FCortex_Final <- cbind(NRC_FCortex_GC, NRC_FCortex_FLT)

NRC_RCortex_GC <- log(apply(NRC_RCortex[ ,grep("GC", colnames(NRC_RCortex))], 1, mean),2)
NRC_RCortex_FLT <- log(apply(NRC_RCortex[ ,grep("FLT", colnames(NRC_RCortex))], 1, mean),2)
NRC_RCortex_Final <- cbind(NRC_RCortex_GC, NRC_RCortex_FLT)

# Select the appropriate column input for the Log2 heatmap
Log2_CA <- Log2_mat[,grep("CA", colnames(Log2_mat))]
Log2_DG <- Log2_mat[,grep("DG", colnames(Log2_mat))]
Log2_FCortex <- Log2_mat[,grep("FCortex", colnames(Log2_mat))]
Log2_RCortex <- Log2_mat[,grep("RCortex", colnames(Log2_mat))]

# Create heatmaps for Log2 values
heatmap_Log2_CA <- Heatmap(Log2_CA, name = "Log2_CA", col = colorRamp2(c(-.75, 0, .75), c("blue", "white", "red")), show_row_names = FALSE, cluster_rows = FALSE)
heatmap_Log2_DG <- Heatmap(Log2_DG, name = "Log2_DG", col = colorRamp2(c(-.75, 0, .75), c("blue", "white", "red")), show_row_names = FALSE, cluster_rows = FALSE)
heatmap_Log2_FCortex <- Heatmap(Log2_FCortex, name = "Log2_FCortex", col = colorRamp2(c(-.75, 0, .75), c("blue", "white", "red")), show_row_names = FALSE, cluster_rows = FALSE)
heatmap_Log2_RCortex <- Heatmap(Log2_RCortex, name = "Log2_RCortex", col = colorRamp2(c(-.75, 0, .75), c("blue", "white", "red")), show_row_names = FALSE, cluster_rows = FALSE)

# Create heatmaps for Normalized Readcount Gene Expression values
heatmap_NRC_CA <- Heatmap(NRC_CA_Final, name = "NRC_CA_Final", col = colorRamp2(c(2, 6), c("white", "black")), show_row_names = FALSE, cluster_rows = FALSE, column_order = c("NRC_CA_FLT","NRC_CA_GC"))
heatmap_NRC_DG <- Heatmap(NRC_DG_Final, name = "NRC_DG_Final", col = colorRamp2(c(2, 6), c("white", "black")), show_row_names = FALSE, cluster_rows = FALSE, column_order = c("NRC_DG_FLT","NRC_DG_GC"))
heatmap_NRC_FCortex <- Heatmap(NRC_FCortex_Final, name = "NRC_FCortex_Final", col = colorRamp2(c(2, 6), c("white", "black")), show_row_names = FALSE, cluster_rows = FALSE, column_order = c("NRC_FCortex_FLT","NRC_FCortex_GC"))
heatmap_NRC_RCortex <- Heatmap(NRC_RCortex_Final, name = "NRC_RCortex_Final", col = colorRamp2(c(2, 6), c("white", "black")), show_row_names = FALSE, cluster_rows = FALSE, column_order = c("NRC_RCortex_FLT","NRC_RCortex_GC"))

# Draw the combined heatmap
png("heatmap_Combined.png", height = 2000, width = 2000, res = 300)
draw(heatmap_Log2_CA + heatmap_NRC_CA + heatmap_Log2_DG + heatmap_NRC_DG + heatmap_Log2_FCortex + heatmap_NRC_FCortex + heatmap_Log2_RCortex + heatmap_NRC_RCortex)
dev.off()

# make histogram
Hist_mat=cbind(NRC_CA_Final, NRC_DG_Final, NRC_FCortex_Final, NRC_RCortex_Final)
png("Histogram_Combined.png", height = 2000, width = 2000, res = 300)
hist(c(Hist_mat)) 
dev.off()
