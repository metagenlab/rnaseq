
library("edgeR")
library("rtracklayer")
library("reshape2")
library("ggplot2")
library("pheatmap")
library("stringr")
library("svglite")
library("ggfortify")

# cutoffs for filtering
FDR_cutoff = as.numeric(snakemake@params$FDR_cutoff)
logFC_cutoff =  as.numeric(snakemake@params$logFC_cutoff)
cond1 <- snakemake@params$cond1
cond2 <- snakemake@params$cond2

#print("FDR_cutoff", FDR_cutoff)
#print(FDR_cutoff/2)

# READ COUNT TABLE
read_counts <- read.table(snakemake@input$count_file, header=T) # , check.names=FALSE

# READ custom_annotation table 
custom_annotation_table <- read.csv(snakemake@params$annotation, header=T, as.is=T, stringsAsFactors = FALSE, sep="\t")

# READ SAMPLE TABLE
sample_info_table <- read.table(snakemake@params$sample_table, header=T, as.is=T, stringsAsFactors = FALSE, sep="\t")
sample_info_table[,"total_reads"] <- NA

# read reference gff
annotations <- readGFF(snakemake@input$reference_gff)
annotations <- annotations[annotations$type == 'gene',]
print("adding rownames")
rownames(annotations) <- annotations$locus_tag
print("ok")
m <- match(annotations$locus_tag, custom_annotation_table$locus_tag)
annotations <- cbind(annotations, custom_annotation_table[m,])

m <- match(rownames(read_counts), annotations$locus_tag)
annotations <- annotations[m,]

print("dim read_counts")
dim(read_counts)
print("dim genes")
dim(annotations)
# get DGEList
y <- DGEList(counts=read_counts, 
             samples=sample_info_table, 
             genes=annotations,
             group=sample_info_table$condition)


keep <- filterByExpr(y)

y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)


print(y$samples)


pdf("test_norm_qc1.pdf")
plotMD(cpm(y, log=TRUE), column=1)
abline(h=0, col="red", lty=2, lwd=2)
dev.off()
pdf("test_norm_qc2.pdf")
plotMD(cpm(y, log=TRUE), column=2)
abline(h=0, col="red", lty=2, lwd=2)
dev.off()

my_data <- cpm(y, log=T)

#pca <- prcomp(t(my_data), scale = F)
#expl.var <- round(pca$sdev^2/sum(pca$sdev^2)*100) # percent explained variance
unscaled = prcomp(t(my_data))
scaled = prcomp(t(my_data), scale = TRUE)
expl.var <- round(unscaled$sdev^2/sum(unscaled$sdev^2)*100)
# Let's plot the results as a screeplot
print("OK")
print(rownames(t(my_data)))
pdf(snakemake@output[[1]])
  #par(mfrow=c(2,2))
  #plot(unscaled)
  #plot(scaled)
  pc <- c(1,2)
  plot(unscaled$x[,pc], 
  xlab=paste0("PC ", pc[1], " (", expl.var[pc[1]], "%)"), 
  ylab=paste0("PC ", pc[2], " (", expl.var[pc[2]], "%)"),
  pch=20,
  col=factor(sample_info_table$condition)
  )
  text(unscaled$x[,pc], labels=rownames(t(my_data)), cex = 1.5)
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
#print("number of match after foilter")
#print(length(w))
tT$table$signif <- FALSE
tT$table$signif[w] <- TRUE
print("locus tags")
locus_tags <- rownames(deg.list)
print("locus tags ok")
# select genes that have 10% false discovery rate
top.deg <- locus_tags[deg.list$FDR < FDR_cutoff]

## VOLCANO PLOT: cutoffs of FRD and logFC as defined in config file
pdf(snakemake@output[[2]])
  plot(deg.list$logFC, -log10(deg.list$FDR), 
  pch=20, main=paste(cond1, "vs", cond2, "comparison"),
  xlab = "Log2 Fold Change", ylab = "-log10(FDR)", 
  xlim=c(-max(abs(deg.list$logFC))*1.1, max(abs(deg.list$logFC))*1.1))
  with(subset(deg.list, FDR < FDR_cutoff & abs(deg.list$logFC) > logFC_cutoff), 
  points(logFC, -log10(FDR), pch=20, col="red"))
  # cutoff for label hard coded to 2xfoldchange
  s <- subset(deg.list, FDR < FDR_cutoff & abs(deg.list$logFC) > 2)
  #text(s$logFC, -log10(s$FDR), labels = paste(s$locus_tag, s$KO))
  abline(v=1, lty=2, lwd=1, col="lightblue")
  abline(v=-1, lty=2, lwd=1, col="lightblue")
dev.off()

## plot SMEAR
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

pdf("pseudocounts.pdf")
ggplot(dfFilter, aes(x = value, colour = variable, fill = variable)) +
  geom_histogram(binwidth = 0.6) + facet_wrap(~ condition) +
  theme(legend.position = "top") + xlab("log2(counts)") + ggtitle("log2(counts) distribution")
dev.off()


heatmap.plotting.replicates <- function(x, name){
	# calculate the spearman correlation on your samples
  spearman.tissue <- melt(cor(x, method = "pearson")) # spearman
	colnames(spearman.tissue)[3] <- "spear.corr"
  p <- ggplot(spearman.tissue, aes(Var1, Var2)) +
  geom_tile(aes(fill = spear.corr, color = "gray", width=0.95, height=0.95)) + 
  geom_text(aes(label = round(spear.corr, 4)), size = 4.3) +
  scale_fill_gradient(low = "white", high = "red") +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle=30,hjust=1,vjust=1.0,size = 14),
        axis.text.y = element_text(size = 14))
	return(p)
}

pdf("pseudocounts_pearson.pdf")
#  pearson
heatmap.plotting.replicates(pseudo_counts, "heatmap-pearson")
dev.off()

pdf("pseudocounts_boxplots.pdf")
ggplot(dfFilter, aes(x=sample, y=value)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4)+ facet_wrap(~ condition)
dev.off()

annotations$gene_length <- (annotations$end - annotations$start) + 1
print("match")
m <- match(rownames(de$genes), rownames(annotations))
print("match ok")
y$genes$length <- annotations$gene_length[m]

genes_rpkm <- rpkm(y)

considered_samples <- sample_info_table$SampleName[sample_info_table$condition %in% c(cond1, cond2)]
considered_samples <- str_replace(considered_samples, "-", ".")
print("------dim rpkm-----")
dim(genes_rpkm)
head(genes_rpkm)
if (length(deg.list[,1])>=100){
log_rpkm.topgenes <- log2(genes_rpkm[rownames(deg.list[1:100,]), considered_samples]+1)
} else {
log_rpkm.topgenes <- log2(genes_rpkm[rownames(deg.list[1:length(deg.list[,1]),]), considered_samples]+1)
}

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

#print("heatmap top 100")

#print(head(log_rpkm.topgenes[,considered_samples]))

svglite(snakemake@output[[4]], height=23, width=w)
pheatmap(log_rpkm.topgenes[,considered_samples], 
         cellwidth = 12, 
         cellheight = 10, 
         annotation_row=metadata_gene,
         main = paste('Heatmap of top 100 genes:', cond1, " vs ", cond2)) # , annotation_col=my_sample_col
dev.off()

#################################
###  Write table ################
#################################



write.table(tT, snakemake@output[[5]], sep="\t")

# filter according to wildcards cutoffs for FRD and logFC
w <- which(tT$table$FDR < FDR_cutoff & abs(tT$table$logFC) > logFC_cutoff)

locus_list <- tT$table$locus_tag[w] 
w_up <- which(tT$table$FDR < FDR_cutoff & abs(tT$table$logFC) > logFC_cutoff & tT$table$logFC > 0 )
w_down <- which(tT$table$FDR < FDR_cutoff & abs(tT$table$logFC) > logFC_cutoff & tT$table$logFC < 0 )
locus_list_up <- tT$table$locus_tag[w_up]
locus_list_down <- tT$table$locus_tag[w_down]
write.table(locus_list_down, snakemake@output[[8]], sep="\t", col.names =F, row.names =F, quote = F)
write.table(locus_list_up, snakemake@output[[9]], sep="\t", col.names =F, row.names =F, quote = F)

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

svglite(snakemake@output[[6]], height=h, width=11)
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
svglite(snakemake@output[[7]], height=h, width=w)
pheatmap(log2(genes_rpkm_filtered_ordered_up[,considered_samples] + 1), 
         cellwidth = 20, 
         cellheight = 14, 
         cluster_rows=TRUE, 
         cluster_cols=FALSE, 
         main = paste('Top ', length(genes_rpkm_filtered_ordered_up[,1]) ,' up:', cond2, " vs ", cond1)) # , annotation_col=my_sample_col
dev.off()
