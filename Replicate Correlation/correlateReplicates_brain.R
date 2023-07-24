# Input data
args=commandArgs(trailingOnly=T)
#data=read.delim(args[1])
name=as.character(args[2])

# Read the matrix from a file
mat <- as.matrix(read.table(args[1], header = TRUE))

# Take the log of each element in the matrix, adding 1 to each element first
dataset <- log(mat[,seq.int(1,3)] + 1)

cormat <- round(cor(dataset),3)
library(reshape2)
melted_cormat <- melt(cormat)
library(ggplot2)

# Make scatterplot
png(paste(name,"_FLT_TX.png",sep = ""))
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "yellow", high = "blue", mid = "white",
   midpoint = 0.70, limit = c(0.5,1), space = "Lab",
   name="Pearson\nCorrelation") +
  theme_minimal()+
 theme(axis.text.x = element_text(angle = 45, vjust = 1,
    size = 12, hjust = 1))+
 coord_fixed()
dev.off()



# Take the log of each element in the matrix, adding 1 to each element first
dataset <- log(mat[,seq.int(4,6)] + 1)

cormat <- round(cor(dataset),3)
library(reshape2)
melted_cormat <- melt(cormat)
library(ggplot2)

# Make scatterplot
png(paste(name,"_CG-SAL.png",sep = ""))
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "yellow", high = "blue", mid = "white",
   midpoint = 0.70, limit = c(0.5,1), space = "Lab",
   name="Pearson\nCorrelation") +
  theme_minimal()+
 theme(axis.text.x = element_text(angle = 45, vjust = 1,
    size = 12, hjust = 1))+
 coord_fixed()
dev.off()

# Take the log of each element in the matrix, adding 1 to each element first
dataset <- log(mat[,seq.int(7,9)] + 1)

cormat <- round(cor(dataset),3)
library(reshape2)
melted_cormat <- melt(cormat)
library(ggplot2)

# Make scatterplot
png(paste(name,"_GC-TX.png",sep = ""))
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "yellow", high = "blue", mid = "white",
   midpoint = 0.70, limit = c(0.5,1), space = "Lab",
   name="Pearson\nCorrelation") +
  theme_minimal()+
 theme(axis.text.x = element_text(angle = 45, vjust = 1,
    size = 12, hjust = 1))+
 coord_fixed()
dev.off()

# Take the log of each element in the matrix, adding 1 to each element first
dataset <- log(mat[,seq.int(10,12)] + 1)

cormat <- round(cor(dataset),3)
library(reshape2)
melted_cormat <- melt(cormat)
library(ggplot2)

# Make scatterplot
png(paste(name,"_FLT-SAL.png",sep = ""))
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "yellow", high = "blue", mid = "white",
   midpoint = 0.70, limit = c(0.5,1), space = "Lab",
   name="Pearson\nCorrelation") +
  theme_minimal()+
 theme(axis.text.x = element_text(angle = 45, vjust = 1,
    size = 12, hjust = 1))+
 coord_fixed()
dev.off()

#FLT-TX, CG-SAL, GC-TX, FLT-SAL
