
library("edgeR")
library("rtracklayer")
library("stringr")


# READ COUNT TABLE
read_counts <- read.table(snakemake@input$count_file, header=T)

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

# COUNT table
y <- DGEList(counts=read_counts, 
             samples=sample_info_table, 
             genes=annotations,
             group=sample_info_table$condition)

# add gene length to table
annotations$gene_length <- (annotations$end - annotations$start) + 1
m <- match(rownames(y$genes), rownames(annotations))
y$genes$length <- annotations$gene_length[m]

# claculated RPKM
genes_rpkm_raw <- as.data.frame(rpkm(y))

genes_rpkm_raw <- cbind(locus_tag=rownames(genes_rpkm_raw), genes_rpkm_raw)

write.table(genes_rpkm_raw, snakemake@output[[1]], sep="\t", row.names=F)
