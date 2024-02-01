
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
annotations <- annotations[annotations$type %in% c('gene', 'pseudogene'),]
annotations <- annotations[annotations$gene_biotype != 'rRNA' ,]
annotations <- annotations[annotations$gene_biotype != 'tRNA' ,]
annotations <- as.data.frame(annotations)

start_lst <- aggregate(annotations$start, by = list(annotations$locus_tag), min)
end_lst <- aggregate(annotations$end, by = list(annotations$locus_tag), max)

annotations <- annotations[!duplicated(annotations$locus_tag), ]
annotations$start <- start_lst$x
annotations$end <- end_lst$x

m <- match(annotations$locus_tag, custom_annotation_table$locus_tag)

annotations <- cbind(annotations, custom_annotation_table[m,])

rownames(annotations) <- annotations$locus_tag


# COUNT table
y <- DGEList(counts=read_counts, 
             samples=sample_info_table, 
             genes=annotations,
             group=sample_info_table$condition)

# add gene length to table
annotations$gene_length <- (annotations$end - annotations$start) + 1
print(rownames(annotations))
m <- match(rownames(y$genes), rownames(annotations))

y$genes$length <- annotations$gene_length[m]

# claculated RPKM
genes_rpkm_raw <- rpkm(y)

genes_rpkm_raw <- as.data.frame(genes_rpkm_raw)

genes_rpkm_raw <- cbind(locus_tag=rownames(genes_rpkm_raw), genes_rpkm_raw)

write.table(genes_rpkm_raw, snakemake@output[[1]], sep="\t", row.names=F)
