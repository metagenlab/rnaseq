
library("edgeR")
library("rtracklayer")
library("reshape2")
library("ggplot2")
library("pheatmap")
library("stringr")
library("svglite")

split_path <- function(path) {
  if (dirname(path) %in% c(".", path)) return(basename(path))
  return(c(basename(path), split_path(dirname(path))))
}

# cutoffs for filtering
FDR_cutoff = snakemake@params$FDR_cutoff
logFC_cutoff = snakemake@params$logFC_cutoff

print("FRD_cutoff", FDR_cutoff)
print(FDR_cutoff/2)

# READ COUNT TABLE
read_counts <- read.table(snakemake@input$count_file, header=T)

# extract conditions from output path
path <- split_path(snakemake@output[[1]])
cond1 <- strsplit(path[[2]], '_vs_')[[1]][[1]]
cond2 <- strsplit(path[[2]], '_vs_')[[1]][[2]]

# READ custom_annotation table 
custom_annotation_table <- read.csv(snakemake@params$annotation, header=T, as.is=T, stringsAsFactors = FALSE, sep="\t")

# READ SAMPLE TABLE
sample_info_table <- read.table(snakemake@params$sample_table, header=T, as.is=T, stringsAsFactors = FALSE, sep="\t")
sample_info_table[,"total_reads"] <- NA

# read reference gff
annotations <- readGFF(snakemake@input$reference_gff)
annotations <- annotations[annotations$type == 'gene',]
rownames(annotations) <- annotations$locus_tag

m <- match(annotations$locus_tag, custom_annotation_table$locus_tag)

annotations <- cbind(annotations, custom_annotation_table[m,])

# get DGEList
y <- DGEList(counts=read_counts, 
             samples=sample_info_table, 
             genes=annotations,
             group=sample_info_table$condition)


keep <- filterByExpr(y)

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
# cutoffs
w <- which(tT$table$FDR < FDR_cutoff & abs(tT$table$logFC) > logFC_cutoff)
tT$table$signif <- FALSE
tT$table$signif[w] <- TRUE

locus_tags <- rownames(deg.list)
# select genes that have 10% false discovery rate
top.deg <- locus_tags[deg.list$FDR < FDR_cutoff]

## VOLCANO PLOT: cutoffs of 0.01 FDR and 2 logFC
pdf(snakemake@output[[2]])
  plot(deg.list$logFC, -log10(deg.list$FDR), # PValue
      pch=20, main=paste(cond1, "vs", cond2, "comparison"),
      xlab = "Log2 Fold Change", ylab = "-log10(FDR)",
      xlim=c(-max(abs(deg.list$logFC))*1.1, max(abs(deg.list$logFC))*1.1) )
  with(subset(deg.list, FDR < .01 & abs(logFC) > 2), 
      points(logFC, -log10(FDR), pch=20, col="red"))
  s <- subset(deg.list, FDR < .01 & abs(logFC) > 2)
  text(s$logFC, -log10(s$FDR), labels = paste(s$locus_tag, s$KO))
  abline(v=2, lty=2, lwd=1, col="lightblue")
  abline(v=-2, lty=2, lwd=1, col="lightblue")
dev.off()

## VOLCANO PLOT: cutoffs of FRD and logFC as defined in config file
pdf(snakemake@output[[3]])
  plot(deg.list$logFC, -log10(deg.list$FDR), 
  pch=20, main=paste(cond1, "vs", cond2, "comparison"),
  xlab = "Log2 Fold Change", ylab = "-log10(FDR)", 
  xlim=c(-max(abs(deg.list$logFC))*1.1, max(abs(deg.list$logFC))*1.1))
  with(subset(deg.list, FDR < FRD_cutoff & abs(deg.list$logFC) > logFC_cutoff), 
  points(logFC, -log10(FDR), pch=20, col="red"))
  s <- subset(deg.list, FDR < FRD_cutoff & abs(deg.list$logFC) > logFC_cutoff)
  text(s$logFC, -log10(s$FDR), labels = paste(s$locus_tag, s$KO))
  abline(v=1, lty=2, lwd=1, col="lightblue")
  abline(v=-1, lty=2, lwd=1, col="lightblue")
dev.off()

## plot SMEAR
pdf(snakemake@output[[4]])
plotSmear(de, de.tags = top.deg, main = paste(cond1, "vs", cond2, "comparison"))
dev.off()

# Normalize counts
pseudo_counts <- log2(y$counts + 1)
pseudocountsFilter.ggplot <- as.data.frame(pseudo_counts)
dfFilter <- melt(pseudocountsFilter.ggplot)
dfFilter <- data.frame(dfFilter, sample = dfFilter$variable)
m <- match(dfFilter$sample, sample_info_table$OldSampleName)
dfFilter$condition <- sample_info_table$condition[m]

pdf("pseudocounts.pdf")
ggplot(dfFilter, aes(x = value, colour = variable, fill = variable)) +
  geom_histogram(binwidth = 0.6) + facet_wrap(~ condition) +
  theme(legend.position = "top") + xlab("log2(counts)") + ggtitle("log2(counts) distribution")
dev.off()

pdf("pseudocounts_boxplots.pdf")
ggplot(dfFilter, aes(x=sample, y=value)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4)+ facet_wrap(~ condition)
dev.off()

annotations$gene_length <- (annotations$end - annotations$start) + 1
m <- match(rownames(de$genes), rownames(annotations))
y$genes$length <- annotations$gene_length[m]

genes_rpkm <- rpkm(y)

considered_samples <- sample_info_table$SampleName[sample_info_table$condition %in% c(cond1, cond2)]
considered_samples <- str_replace(considered_samples, "-", ".")

log_rpkm.topgenes <- log2(genes_rpkm[rownames(deg.list[1:100,]), considered_samples]+1)
m <- match(rownames(log_rpkm.topgenes), custom_annotation_table$locus_tag)

gene <- custom_annotation_table$gene[m]
product <- custom_annotation_table$product[m]
annot <- paste0(rownames(log_rpkm.topgenes), '/', gene, '/', product)
rownames(log_rpkm.topgenes) <- annot

log_rpkm.topgenes_FC <- tT[custom_annotation_table$locus_tag[m],"logFC"]

metadata_gene <- data.frame(up_or_down=rep("gg", nrow(log_rpkm.topgenes_FC)),
                            row.names=rownames(log_rpkm.topgenes), stringsAsFactors=FALSE
                            )

metadata_gene$up_or_down[which(log_rpkm.topgenes_FC[[1]]$logFC > 0)] <- "up"
metadata_gene$up_or_down[which(log_rpkm.topgenes_FC[[1]]$logFC < 0)] <- "down"
metadata_gene$up_or_down <- as.factor(metadata_gene$up_or_down)

w <- length(considered_samples)*2.5 + nchar(max(rownames(log_rpkm.topgenes))) / 10


#################################
###  top 100 genes heatmap ######
#################################


svglite(snakemake@output[[5]], height=23, width=w)
pheatmap(log_rpkm.topgenes[,considered_samples], 
         cellwidth = 12, 
         cellheight = 10, 
         annotation_row=metadata_gene,
         main = paste('Heatmap of top 100 genes:', cond1, " vs ", cond2)) # , annotation_col=my_sample_col
dev.off()

#################################
###  Write table ################
#################################

write.table(tT, snakemake@output[[6]], sep="\t")

# filter according to config file cutoffs for FRD and logFC
w <- which(tT$table$FDR < FRD_cutoff & abs(tT$table$logFC) > logFC_cutoff)

locus_list <- tT$table$locus_tag[w] 
genes_rpkm_filtered <- as.data.frame(genes_rpkm[locus_list,])
genes_rpkm_filtered$logFC <- tT$table$logFC[w] 

# prepare consensus annotation for each gene
m <- match(rownames(genes_rpkm_filtered), custom_annotation_table$locus_tag)
gene <- custom_annotation_table$gene[m]
product <- custom_annotation_table$product[m]
annot <- paste0(rownames(genes_rpkm_filtered), '/', gene, '/', product)
rownames(genes_rpkm_filtered) <- annot

# order by the absolute value of fold change
genes_rpkm_filtered_ordered <- genes_rpkm_filtered[order(-abs(genes_rpkm_filtered$logFC)),]

#################################
###  top 100 downregulated  #####
#################################

# order by logFC
genes_rpkm_filtered_ordered <- genes_rpkm_filtered_ordered[order(genes_rpkm_filtered_ordered$logFC),]
# keep only downregulated genes 
genes_rpkm_filtered_ordered_down <- genes_rpkm_filtered_ordered[genes_rpkm_filtered_ordered$logFC < 0,]
if (length(genes_rpkm_filtered_ordered_down[,1]) > 100) {
  genes_rpkm_filtered_ordered_down <- genes_rpkm_filtered_ordered_down[1:100,]
}

# calculate plot height and width according to labels size and number of rows 
h <- 23 * length(genes_rpkm_filtered_ordered_down[,1])/100
w <- 5 + nchar(max(rownames(genes_rpkm_filtered_ordered_down))) / 10

svglite(snakemake@output[[7]], height=h, width=11)
pheatmap(log2(genes_rpkm_filtered_ordered_down[,considered_samples] + 1), 
         cellwidth = 20, 
         cellheight = 14, 
         cluster_rows=TRUE, 
         cluster_cols=FALSE, 
         main = paste('Top', length(genes_rpkm_filtered_ordered_down[,1]) ,' down:', cond2, " vs ", cond1)) # , annotation_col=my_sample_col
dev.off()


#################################
###  top 100 upregulated  #######
#################################

# order by logFC
genes_rpkm_filtered_ordered <- genes_rpkm_filtered_ordered[order(-genes_rpkm_filtered_ordered$logFC),]
# keep only upregulated genes
genes_rpkm_filtered_ordered_up <- genes_rpkm_filtered_ordered[genes_rpkm_filtered_ordered$logFC > 0,]
if (length(genes_rpkm_filtered_ordered_up[,1]) > 100) {
  genes_rpkm_filtered_ordered_up <- genes_rpkm_filtered_ordered_up[1:100,]
}
h <- 23 * length(genes_rpkm_filtered_ordered_up[,1])/100
w <- 5 + nchar(max(rownames(genes_rpkm_filtered_ordered_up))) / 10
svglite(snakemake@output[[8]], height=h, width=w)
pheatmap(log2(genes_rpkm_filtered_ordered_up[,considered_samples] + 1), 
         cellwidth = 20, 
         cellheight = 14, 
         cluster_rows=TRUE, 
         cluster_cols=FALSE, 
         main = paste('Top ', length(genes_rpkm_filtered_ordered_up[,1]) ,' up:', cond2, " vs ", cond1)) # , annotation_col=my_sample_col
dev.off()
