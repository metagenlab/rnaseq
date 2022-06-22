

library(svglite)
require(gridExtra)
library(ggplot2)


print(snakemake@input[[1]])
t <- read.table(snakemake@input[[1]], header=T)

p1 <- ggplot(t, aes(x=n_reads, y=larger_than_0, color=sample)) + geom_point() + geom_line() + geom_hline(yintercept=894) + ylab("Number of genes") + ggtitle(">0 reads")
p2 <- ggplot(t, aes(x=n_reads, y=larger_than_10, color=sample)) + geom_point() + geom_line() + geom_hline(yintercept=894) + ylab("Number of genes") + ggtitle(">10 reads")
p3 <- ggplot(t, aes(x=n_reads, y=larger_than_100, color=sample)) + geom_point() + geom_line() + geom_hline(yintercept=894) + ylab("Number of genes") + ggtitle(">100 reads")
p4 <- ggplot(t, aes(x=n_reads, y=larger_than_1000, color=sample)) + geom_point() + geom_line() + geom_hline(yintercept=894) + ylab("Number of genes") + ggtitle(">1000 reads")

svglite(snakemake@output[[1]],height=12, width=12)
grid.arrange(p1, p2, p3, p4, ncol=2)
dev.off()

pdf(snakemake@output[[2]],height=12, width=12)
grid.arrange(p1, p2, p3, p4, ncol=2)
dev.off()
