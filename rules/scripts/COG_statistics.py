
import pandas
import re

FDR_cutoff = snakemake.params["FDR_cutoff"]
logFC_cutoff = snakemake.params["logFC_cutoff"]

locus_tag2COG = pandas.read_csv(snakemake.input["locus_tag2COG"], sep="\t", names=["locus_tag", "COG_corresp"]).set_index("locus_tag")

COG2annotation = pandas.read_csv(snakemake.input["COG2annotation"], sep="\t", names=["COG_corresp","COG_category", "COG_description2"]).set_index("COG_corresp")

COG_category2description = pandas.read_csv(snakemake.input["COG_category2description"], sep="\t", comment='#', names=["COG_category", "COG_category_description"]).set_index("COG_category").to_dict()["COG_category_description"]

# count table
count_table = pandas.read_csv(snakemake.input["count_file"], sep="\t", index_col=0)

df_COG = count_table.join(locus_tag2COG, on='locus_tag')

df_cog_with_annotations = df_COG.join(COG2annotation, on='COG_corresp')

COG_category2n_up = df_cog_with_annotations[['COG_category', 'locus_tag', 'FDR', 'logFC']].query(f'FDR < {FDR_cutoff} and logFC > {logFC_cutoff}')[['COG_category']].groupby("COG_category").size().to_dict()
COG_category2n_down = df_cog_with_annotations[['COG_category', 'locus_tag', 'FDR', 'logFC']].query(f'FDR < {FDR_cutoff} and logFC < {logFC_cutoff}')[['COG_category']].groupby("COG_category").size().to_dict()

print("updown", COG_category2n_up)

COG_category2n_genes = df_cog_with_annotations.groupby("COG_category").size().to_dict()

locus_tag2COG = locus_tag2COG.to_dict()["COG_corresp"]
COG2annotation = COG2annotation.to_dict()["COG_description2"]

print("dico", COG2annotation)

o = open(snakemake.output[0], "w")
p = open(snakemake.output[1], "w")


h = ["COG_category",
     "COG_category_description",
     "locus_tag",
     "COG",
     "COG_description",
     "logFC",
     "FDR",
     "signif"]

p.write("\t".join(h) + '\n')

h = ["COG_category",
     "COG_category_description",
     "n_genes",
     "n_genes_up",
     "n_genes_down",
     "n_genes_up_and_down",
     "fraction_all",
     "fraction_up",
     "fraction_down",
     "fraction_up_and_down"]

o.write("\t".join(h) + '\n')


for COG_category in COG_category2description:
    COG_category_description = COG_category2description[COG_category]
    
    n_genes = float(sum(COG_category2n_genes.values()))
    n_genes_up = float(sum(COG_category2n_up.values()))
    n_genes_down = float(sum(COG_category2n_down.values()))
    n_genes_up_and_down = n_genes_up + n_genes_down
    
    
    try:
        n_genes_category = COG_category2n_genes[COG_category]
        fracton_category = float(n_genes_category)/float(n_genes)
    except KeyError:
        # pathway absent from the considered genome
        continue
    try:
        n_genes_up_category = COG_category2n_up[COG_category]
    except KeyError:
        n_genes_up_category = 0
    try:
        n_genes_down_category = COG_category2n_down[COG_category]
    except KeyError:
        n_genes_down_category = 0

    n_up_and_down_category = n_genes_up_category + n_genes_down_category
    
    # calculate fraction if not null 
    if n_genes_up == 0:
        fraction_up = 0
    else:
        fraction_up = float(n_genes_up_category)/n_genes_up
    if n_genes_down == 0:
        fraction_down = 0
    else:
        fraction_down = float(n_genes_down_category)/n_genes_down
    if n_genes_up_and_down == 0:
        fraction_diff = 0
    else:
        fraction_diff = float(n_up_and_down_category)/n_genes_up_and_down

    o.write(f"{COG_category}\t{COG_category_description}\t{n_genes_category}\t{n_genes_up_category}\t{n_genes_down_category}\t{n_up_and_down_category}\t{fracton_category}\t{fraction_up}\t{fraction_down}\t{fraction_diff}\n")

    locus_list = df_cog_with_annotations[df_cog_with_annotations["COG_category"] == COG_category]["locus_tag"]


    for locus in locus_list:
        COG = locus_tag2COG[locus]
        COG_description = COG2annotation[COG]
        try:
            logFC  = count_table[count_table["locus_tag"] == locus]["logFC"][0]
            FDR  = count_table[count_table["locus_tag"] == locus]["FDR"][0]
            product  = ' '.join(str(count_table[count_table["locus_tag"] == locus]["product.1"][0]).split())
            signif  = count_table[count_table["locus_tag"] == locus]["signif"][0]
        except:
            # case very lowly expressed genes removed from the table
            logFC, FDR, product, signif = None, None, None, None
        p.write(f"{COG_category}\t{COG_category_description}\t{locus}\t{product}\t{COG}\t{COG_description}\t{logFC}\t{FDR}\t{signif}\n")


