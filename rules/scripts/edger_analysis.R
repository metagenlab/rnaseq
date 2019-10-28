
library("edgeR")
library("rtracklayer")
library("reshape2")
library("ggplot2")
library("pheatmap")

split_path <- function(path) {
  if (dirname(path) %in% c(".", path)) return(basename(path))
  return(c(basename(path), split_path(dirname(path))))
}

# READ COUNT TABLE
read_counts <- read.table(snakemake@input$count_file, header=T)
head(read_counts )
# extract conditions from output path
path <- split_path(snakemake@output[[1]])
cond1 <- strsplit(path[[2]], '_vs_')[[1]][[1]]
cond2 <- strsplit(path[[2]], '_vs_')[[1]][[2]]

# READ custom_annotation table 
custom_annotation_table <- read.csv(snakemake@params$annotation, header=T, as.is=T, stringsAsFactors = FALSE, sep="\t")
head(custom_annotation_table)

# READ SAMPLE TABLE
sample_info_table <- read.table(snakemake@params$sample_table, header=T, as.is=T, stringsAsFactors = FALSE, sep="\t")
sample_info_table[,"total_reads"] <- NA

# read reference gff
annotations <- readGFF(snakemake@input$reference_gff)
annotations <- annotations[annotations$type == 'gene',]
rownames(annotations) <- annotations$locus_tag

m <- match(annotations$locus_tag, custom_annotation_table$locus_tag)

print(m)
annotations <- cbind(annotations, custom_annotation_table[m,])

print(head(annotations))
dim(annotations)

# get DGEList
y <- DGEList(counts=read_counts, 
             samples=sample_info_table, 
             genes=annotations,
             group=sample_info_table$condition)


keep <- filterByExpr(y)
print("KEEP:")
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y$samples

pdf(snakemake@output[[1]])
plotMDS(y, top = 1000, labels = y$samples$SampleName, col = as.numeric(y$samples$group), 
        pch = as.numeric(y$samples$group), cex = 2)
dev.off()

d <- estimateCommonDisp(y) 
d <- estimateTagwiseDisp(d)

de <- exactTest(d, pair = c(cond2, cond1))

#hist(de$table$PValue, breaks = 50, xlab = 'p-value (without correction)')

# compute FDR all genes, save a sorted table of differentially expressed genes.
# gathering differential expressed genes

tT <- topTags(de, n = nrow(d))
# tabular form of differentially expressed genes 
deg.list <- tT$table

locus_tags <- rownames(deg.list)
# select genes that have 10% false discovery rate
top.deg <- locus_tags[deg.list$FDR < .01]
print("deg.list$")
head(deg.list)
## For students - VOLCANO PLOT
pdf(snakemake@output[[2]])
plot(deg.list$logFC, -log10(deg.list$PValue), 
     pch=20, main=paste(cond1, "vs", cond2, "comparison"),
     xlab = "Log2 Fold Change", ylab = "-log10(pvalue)")
with(subset(deg.list, FDR<.01 & abs(logFC)>2), 
     points(logFC, -log10(PValue), pch=20, col="red"))
s <- subset(deg.list, FDR<.01 & abs(logFC)>2)
text(s$logFC, -log10(s$PValue), labels = paste(s$locus_tag, s$KO))
abline(v=2, lty=2, lwd=1, col="lightblue")
abline(v=-2, lty=2, lwd=1, col="lightblue")
dev.off()

pdf(snakemake@output[[3]])
plotSmear(de, de.tags = top.deg, main = paste(cond1, "vs", cond2, "comparison"))
dev.off()


# Normalize counts
pseudo_counts <- log2(y$counts + 1)
pseudocountsFilter.ggplot <- as.data.frame(pseudo_counts)
dfFilter <- melt(pseudocountsFilter.ggplot)
dfFilter <- data.frame(dfFilter, sample = dfFilter$variable)
m <- match(dfFilter$sample, sample_info_table$OldSampleName)
dfFilter$condition <- sample_info_table$condition[m]

pdf("psuedocounts.pdf")
ggplot(dfFilter, aes(x = value, colour = variable, fill = variable)) +
  geom_histogram(binwidth = 0.6) + facet_wrap(~ condition) +
  theme(legend.position = "top") + xlab("log2(counts)") + ggtitle("log2(counts) distribution")
dev.off()

pdf("psuedocounts_boxplots.pdf")
ggplot(dfFilter, aes(x=sample, y=value)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4)+ facet_wrap(~ condition)
dev.off()


annotations$gene_length <- (annotations$end - annotations$start) + 1
m <- match(rownames(de$genes), rownames(annotations))
y$genes$length <- annotations$gene_length[m]
head(y$genes)
print("rpkm")
genes_rpkm <- rpkm(y)
print("ok")
head(genes_rpkm)


considered_samples <- sample_info_table$OldSampleName[sample_info_table$condition %in% c(cond1, cond2)]
print("considered_samples")
print(considered_samples)
print("topgenes")
log_rpkm.topgenes <- log2(genes_rpkm[rownames(deg.list[1:100,]), considered_samples]+1)
m <- match(rownames(log_rpkm.topgenes), custom_annotation_table$locus_tag)
gene <- custom_annotation_table$gene[m]
product <- custom_annotation_table$product[m]
annot <- paste0(rownames(log_rpkm.topgenes), '/', gene, '/', product)
rownames(log_rpkm.topgenes) <- annot

pdf(paste0("pheatmap_", cond1, "_vs_", cond2, ".pdf"), height=17, width=10)
#my_sample_col <- data.frame(sample = sample_info_table$condition)
#row.names(my_sample_col) <- colnames(log_rpkm.topgenes)
#print(log_rpkm.topgenes)
pheatmap(log_rpkm.topgenes, main = paste('Heatmap of top 100 genes:', cond1, " vs ", cond2)) # , annotation_col=my_sample_col
dev.off()
# WRITE table
write.table(tT, snakemake@output[[4]], sep="\t")
