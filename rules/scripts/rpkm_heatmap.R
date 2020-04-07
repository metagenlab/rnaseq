
library("edgeR")
library("rtracklayer")
library("reshape2")
library("ggplot2")
library("pheatmap")
library("stringr")
library("svglite")

# READ custom_annotation table 
counts <- read.csv(snakemake@input[[1]], header=T, as.is=T, stringsAsFactors = FALSE, sep="\t")

subset <- counts[,2:length(counts)]
row.names(subset) <- counts[,1]

svglite(snakemake@output[[1]], height=23, width=9)
pheatmap(as.matrix(log2(subset)), 
         cellwidth = 12, 
         cellheight = 10, 
         main = "Heatmap") # , annotation_col=my_sample_col
dev.off()

