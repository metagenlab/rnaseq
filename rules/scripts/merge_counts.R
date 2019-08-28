library("reshape2")
library("ggplot2")

split_path <- function(path) {
  if (dirname(path) %in% c(".", path)) return(basename(path))
  return(c(basename(path), split_path(dirname(path))))
}

# LOAD AND COMBINE COUNTS INTO A SINGLE DATAFRAME
htseq_count_file_list <- snakemake@input$htseq_files

first_file <- htseq_count_file_list[[1]]
sample_name <- split_path(first_file)[[4]]
read_counts <- read.table(first_file, header=F, as.is=T)
rownames(read_counts) <- read_counts[,1]

colnames(read_counts) <- c("locus_tag", sample_name)

for (path in htseq_count_file_list[c(2:length(htseq_count_file_list))]) {
    sample_name <- split_path(path)[[4]]
    d = read.table(path, header=F, as.is=T)
    read_counts[,sample_name] <- d[,2]
}

# remove locus tag raw
read_counts$locus_tag <- NULL
read_counts <- subset(read_counts, !(rownames(read_counts) %in% c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")))

write.table(read_counts, snakemake@output[[1]], sep="\t")


heatmap.plotting.replicates <- function(x, name){
	# calculate the spearman correlation on your samples
  spearman.tissue <- melt(cor(x, method = "spearman"))
	colnames(spearman.tissue)[3] <- "spear.corr"

	qp <- qplot(x=Var1, y=Var2, 
	            data=spearman.tissue, fill=spear.corr, 
	            geom="tile", xlab = "", ylab = "") + 
	  theme(panel.grid.major = element_blank(), 
	        panel.grid.minor = element_blank(), 
	        panel.background = element_blank(), 
	        axis.text.x = element_text(angle = 45, hjust = 1)) + 
	  labs(title = name)

	return(qp)
}


pdf(snakemake@output[[2]])
heatmap.plotting.replicates(read_counts, "Heatmap samples similarity")
dev.off()