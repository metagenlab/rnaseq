
library("edgeR")
library("rtracklayer")
library("reshape2")
library("ggplot2")
library("pheatmap")
library("stringr")
library("svglite")
library("ComplexHeatmap")

# READ custom_annotation table 
counts <- read.csv(snakemake@input[[1]], header=T, as.is=T, stringsAsFactors = FALSE, sep="\t")

subset <- counts[,c(4:length(counts))]
#row.names(subset) <- counts[,3]
print(head(subset))
print(dim(subset))
print(counts$label)
ha <- rowAnnotation(foo = counts$label)
print("ok")
svglite(snakemake@output[[1]], height=70, width=8)
Heatmap(as.matrix(subset), right_annotation = ha, split =counts$label, show_row_dend = FALSE, rect_gp = gpar(col = "white", lwd = 2))
dev.off()

