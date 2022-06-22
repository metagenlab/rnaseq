
library(ggplot2)
library(reshape2)

t <- read.csv(snakemake@input[[1]]   , sep='\t', header=1, stringsAsFactors=FALSE)
# COG_category	COG_category_description	n_genes	n_genes_up	n_genes_down	n_genes_up_and_down	fraction_all	fraction_up	fraction_down	fraction_up_and_down
# exclude "NO cog" category from the plot (generally very large, distort the plot)
w <- which(t$COG_category_description != "No COG annotation")
sub_t <- t[w,c("COG_category_description", "fraction_all", "fraction_up", "fraction_down", "fraction_up_and_down")]

long <- melt(sub_t, id.vars = c("COG_category_description"))

long$variable <- as.factor(long$variable)

p <- ggplot(long, aes(y=value, x=COG_category_description, fill=variable)) + geom_bar(position = "dodge2", stat="identity")
p <- p + coord_flip() + ylab("COG fraction")
p <- p + scale_x_discrete(limits = rev(levels(factor(long$COG_category_description))))
#p <- p + scale_fill_brewer(palette="Set1")

ggsave(
  snakemake@output[[1]],
  plot = last_plot(),
  width = 30,
  height = 20,
  units = c("cm"),
  dpi = 300)
