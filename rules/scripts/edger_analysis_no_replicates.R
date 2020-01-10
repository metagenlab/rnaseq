
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

annotations <- cbind(annotations, custom_annotation_table[m,])

# COUNT table
y <- DGEList(counts=read_counts, 
             samples=sample_info_table, 
             genes=annotations,
             group=sample_info_table$condition)

# filter genes with low counts
keep <- filterByExpr(y)

# plot MDS based on the top 1000 genes y$samples$SampleName
pdf(snakemake@output[[1]])
plotMDS(y, top = 1000, labels = NULL, col = as.numeric(y$samples$group), 
        pch = as.numeric(y$samples$group), cex = 1, xlim=c(-5,2), ylim=c(-2,1.5), legend = "auto")
legend("topleft", legend=y$samples$SampleName, pch=as.numeric(y$samples$group), col=as.numeric(y$samples$group))
dev.off()

# add gene length to table
annotations$gene_length <- (annotations$end - annotations$start) + 1
m <- match(rownames(y$genes), rownames(annotations))
y$genes$length <- annotations$gene_length[m]

# claculated RPKM
genes_rpkm_raw <- as.data.frame(rpkm(y))

genes_rpkm <- genes_rpkm_raw[,c(cond1, cond2)]

print("adding fold change")
genes_rpkm$LogFoldChange <- log2(genes_rpkm[,cond1]/genes_rpkm[,cond2])
#genes_rpkm$median <- apply(genes_rpkm_raw, 1, median, na.rm = TRUE)
genes_rpkm$mean_all <- rowMeans(genes_rpkm_raw)
genes_rpkm$log2_mean_all <- log2(genes_rpkm$mean_all + 1)
w <- which(genes_rpkm[,cond2] == 0 & genes_rpkm[,cond1] != 0)
genes_rpkm[, paste("only", cond1, sep="_")] <- FALSE
genes_rpkm[, paste("only", cond1, sep="_")][w] <- TRUE 

w <- which(genes_rpkm[,cond1] == 0 & genes_rpkm[,cond2] != 0)
genes_rpkm[, paste("only", cond2, sep="_")] <- FALSE
genes_rpkm[, paste("only", cond2, sep="_")][w] <- TRUE 

# calculated fold change based on rpkm

# sort by expression and fold change


# plot median RPKM vs fold change
pdf(snakemake@output[[2]], width=12, height=6)
plot(genes_rpkm$log2_mean_all, genes_rpkm$LogFoldChange, main = paste(cond1, "vs", cond2, "comparison"), pch = 16)
abline(h=2, col="red")
abline(h=-2, col="red")
abline(v=2, col="blue")

w1 <- which((genes_rpkm$LogFoldChange >= 2 & genes_rpkm$log2_mean_all > 2) | (genes_rpkm$LogFoldChange <= -2 & genes_rpkm$log2_mean_all > 2)) 

points(genes_rpkm$log2_mean_all[w1], genes_rpkm$LogFoldChange[w1], col="green", pch = 16)

dev.off()

###################################################
# plot lomarlized counts distribution and boxplots 
###################################################

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

genes_rpkm_filtered <- genes_rpkm[w1,]
# order by the absolute value of fold change
genes_rpkm_filtered_ordered <- genes_rpkm_filtered[order(-abs(genes_rpkm_filtered$LogFoldChange)),]


m <- match(rownames(genes_rpkm_filtered_ordered), custom_annotation_table$locus_tag)
gene <- custom_annotation_table$gene[m]
product <- custom_annotation_table$product[m]
annot <-paste0(rownames(genes_rpkm_filtered_ordered), '/', gene, '/', product) # strtrim
rownames(genes_rpkm_filtered_ordered) <- annot

w1 <- which(genes_rpkm_filtered_ordered$LogFoldChange > 0)
w2 <- which(genes_rpkm_filtered_ordered$LogFoldChange < 0)

metadata_gene <- data.frame(up_or_down=rep("gg", nrow(genes_rpkm_filtered_ordered)),
                            row.names=rownames(genes_rpkm_filtered_ordered), stringsAsFactors=FALSE
                            )

metadata_gene$up_or_down[w1] <- "up"
metadata_gene$up_or_down[w2] <- "down"

metadata_gene$up_or_down <- as.factor(metadata_gene$up_or_down)

# top genes merged
pdf(snakemake@output[[4]], height=23, width=11)
pheatmap(log2(genes_rpkm_filtered_ordered[1:100,c(cond2, cond1)] + 1), 
         cellwidth = 20, 
         cellheight = 14, 
         annotation_row=metadata_gene, main = paste(cond2, " vs ", cond1))

dev.off()

# add annotations and column tagging significant genes  
m <- match(rownames(genes_rpkm), custom_annotation_table$locus_tag)
genes_rpkm <- cbind(genes_rpkm, custom_annotation_table[m,])
genes_rpkm$signif <- FALSE 
genes_rpkm$signif[w1] <- TRUE

# write table 
write.table(genes_rpkm, snakemake@output[[3]], sep="\t", col.names=NA)


###  top 100 downregulated  ###
# order by LogFoldChange
genes_rpkm_filtered_ordered <- genes_rpkm_filtered_ordered[order(genes_rpkm_filtered_ordered$LogFoldChange),]
# keep only downregulated genes 
genes_rpkm_filtered_ordered_down <- genes_rpkm_filtered_ordered[genes_rpkm_filtered_ordered$LogFoldChange < 0,]
if (length(genes_rpkm_filtered_ordered_down[,1]) > 100) {
  genes_rpkm_filtered_ordered_down <- genes_rpkm_filtered_ordered_down[1:100,]
}

# calculate plot height and width according to labels size and number of rows 
h <- 23 * length(genes_rpkm_filtered_ordered_down[,1])/100
w <- 4 + nchar(max(rownames(genes_rpkm_filtered_ordered_down))) / 10

pdf(snakemake@output[[5]], height=h, width=11)
pheatmap(log2(genes_rpkm_filtered_ordered_down[,c(cond2,cond1)] + 1), 
         cellwidth = 20, 
         cellheight = 14, 
         cluster_rows=TRUE, 
         cluster_cols=FALSE, 
         main = paste('Top', length(genes_rpkm_filtered_ordered_down[,1]) ,' down:', cond2, " vs ", cond1)) # , annotation_col=my_sample_col
dev.off()


###  top 100 upregulated  ###
# order by LogFoldChange
genes_rpkm_filtered_ordered <- genes_rpkm_filtered_ordered[order(-genes_rpkm_filtered_ordered$LogFoldChange),]
# keep only upregulated genes
genes_rpkm_filtered_ordered_up <- genes_rpkm_filtered_ordered[genes_rpkm_filtered_ordered$LogFoldChange > 0,]
if (length(genes_rpkm_filtered_ordered_up[,1]) > 100) {
  genes_rpkm_filtered_ordered_up <- genes_rpkm_filtered_ordered_up[1:100,]
}
h <- 23 * length(genes_rpkm_filtered_ordered_up[,1])/100
w <- 4 + nchar(max(rownames(genes_rpkm_filtered_ordered_up))) / 10
pdf(snakemake@output[[6]], height=h, width=w)
pheatmap(log2(genes_rpkm_filtered_ordered_up[,c(cond2,cond1)] + 1), 
         cellwidth = 20, 
         cellheight = 14, 
         cluster_rows=TRUE, 
         cluster_cols=FALSE, 
         main = paste('Top ', length(genes_rpkm_filtered_ordered_up[,1]) ,' up:', cond2, " vs ", cond1)) # , annotation_col=my_sample_col
dev.off()

