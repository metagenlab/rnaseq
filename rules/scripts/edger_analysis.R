
library("edgeR")
library("rtracklayer")
library("reshape2")
library("ggplot2")

split_path <- function(path) {
  if (dirname(path) %in% c(".", path)) return(basename(path))
  return(c(basename(path), split_path(dirname(path))))
}

# READ COUNT TABLE
read_counts <- read.table(snakemake@input$count_file, header=T)


# extract conditions from output path
path <- split_path(snakemake@output[[1]])
cond1 <- strsplit(path[[2]], '_vs_')[[1]][[1]]
cond2 <- strsplit(path[[2]], '_vs_')[[1]][[2]]

# READ SAMPLE TABLE
sample_info_table <- read.table(snakemake@params$sample_table, header=T)
sample_info_table[,"total_reads"] <- NA

# read reference gff
annotations <- readGFF(snakemake@input$reference_gff)
annotations <- annotations[annotations$type == 'gene',]
rownames(annotations) <- annotations$locus_tag
print(head(annotations))
dim(annotations)

# get DGEList
y <- DGEList(counts=read_counts, 
             samples=sample_info_table, 
             genes=annotations,
             group=sample_info_table$condition)


keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y$samples

# PLOT histotgramm rawcounts 
#pdf("counts_distrib.pdf")
#ggplot(read_counts, aes(x = "AOZQ-2")) + geom_histogram(fill = "#525252", binwidth = 2000) + ggtitle("Histogram on Raw Counts for WL 1 sample")
#dev.off()

pdf(snakemake@output[[1]])
plotMDS(y, top = 1000, labels = y$samples$SampleName, col = as.numeric(y$samples$group), 
        pch = as.numeric(y$samples$group), cex = 2)
dev.off()

d <- estimateCommonDisp(y) 
d <- estimateTagwiseDisp(d)

de <- exactTest(d, pair = c(cond1, cond2))

#hist(de$table$PValue, breaks = 50, xlab = 'p-value (without correction)')

# compute FDR all genes, save a sorted table of differentially expressed genes.
# gathering differential expressed genes

tT <- topTags(de, n = nrow(d))
# tabular form of differentially expressed genes 
deg.list <- tT$table

locus_tags <- rownames(deg.list)
# select genes that have 1% false discovery rate
top.deg <- locus_tags[deg.list$FDR < .01]

## For students - VOLCANO PLOT
pdf(snakemake@output[[2]])
plot(deg.list$logFC, -log10(deg.list$PValue), 
     pch=20, main=paste(cond1, "vs", cond2, "comparison"), xlim=c(-4,4),
     xlab = "Log2 Fold Change", ylab = "-log10(pvalue)")
with(subset(deg.list, FDR<.01 & abs(logFC)>2), 
     points(logFC, -log10(PValue), pch=20, col="lightblue"))
dev.off()

pdf(snakemake@output[[3]])
plotSmear(de, de.tags = top.deg, main = paste(cond1, "vs", cond2, "comparison"))
dev.off()

# WRITE table
print(tT)
write.table(tT, snakemake@output[[4]], sep="\t")