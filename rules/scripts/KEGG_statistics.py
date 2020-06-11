
import pandas
import re

FDR_cutoff = snakemake.params["FDR_cutoff"]
logFC_cutoff = snakemake.params["logFC_cutoff"]

# b0114	map00010
locus_tag2pathway = pandas.read_csv(snakemake.input["locus_tag2pathway"], sep="\t", names=["locus_tag","pathway"])

# b2388	M00001
locus_tag2module = pandas.read_csv(snakemake.input["locus_tag2module"], sep="\t", names=["locus_tag","module"])

# M00001	Glycolysis (Embden-Meyerhof pathway), glucose => pyruvate
module2description = pandas.read_csv(snakemake.input["module2description"], sep="\t", names=["module","module_description"]).set_index("module").to_dict()["module_description"]

# map00010	Glycolysis / Gluconeogenesis
pathway2description = pandas.read_csv(snakemake.input["pathway2description"], sep="\t", names=["pathway","pathway_description"]).set_index("pathway").to_dict()["pathway_description"]


locus_tag2ko = pandas.read_csv(snakemake.input["locus_tag2ko"], sep="\t", names=["locus_tag","KO"]).set_index("locus_tag").to_dict()["KO"]

ko2description = pandas.read_csv(snakemake.input["ko2description"], sep="\t", names=["KO","description"]).set_index("KO").to_dict()["description"]

# count table
count_table = pandas.read_csv(snakemake.input["count_file"], sep="\t", index_col=0)

df_module = count_table.join(locus_tag2module.set_index('locus_tag'), on='locus_tag')
df_pathway = count_table.join(locus_tag2pathway.set_index('locus_tag'), on='locus_tag')

module2n_up = df_module[['module', 'locus_tag', 'FDR', 'logFC']].query(f'FDR < {FDR_cutoff} and logFC > {logFC_cutoff}')[['module']].groupby("module").size().to_dict()
module2n_down = df_module[['module', 'locus_tag', 'FDR', 'logFC']].query(f'FDR < {FDR_cutoff} and logFC < {logFC_cutoff}')[['module']].groupby("module").size().to_dict()

pathway2n_up = df_pathway[["pathway", 'locus_tag', 'FDR', 'logFC']].query(f'FDR < {FDR_cutoff} and logFC > {logFC_cutoff}')[["pathway"]].groupby("pathway").size().to_dict()
pathway2n_down = df_pathway[["pathway", 'locus_tag', 'FDR', 'logFC']].query(f'FDR < {FDR_cutoff} and logFC < {logFC_cutoff}')[["pathway"]].groupby("pathway").size().to_dict()

module2n_genes = locus_tag2module.groupby("module").size().to_dict()
pathway2n_genes = locus_tag2pathway.groupby("pathway").size().to_dict()

o = open(snakemake.output[0], "w")
p = open(snakemake.output[1], "w")
q = open(snakemake.output[2], "w")
r = open(snakemake.output[3], "w") 


h = ["pathway",
     "pathway_description",
     "locus_tag",
     "product",
     "ko",
     "ko_description",
     "logFC",
     "FDR",
     "signif"]

q.write("\t".join(h) + '\n')

h = ["module",
     "module_description",
     "locus_tag",
     "product",
     "ko",
     "ko_description",
     "logFC",
     "FDR",
     "signif"]

r.write("\t".join(h) + '\n')

h = ["pathway",
     "pathway_description",
     "n_genes",
     "n_genes_up",
     "n_genes_down",
     "n_genes_up_and_down",
     "fraction_up",
     "fraction_down",
     "fraction_up_and_down"]

o.write("\t".join(h) + '\n')

for pathway in pathway2description:
    pathway_description = pathway2description[pathway]
    try:
        n_genes = pathway2n_genes[pathway]
    except KeyError:
        # pathway absent from the considered genome
        continue
    try:
        n_genes_up = pathway2n_up[pathway]
    except KeyError:
        n_genes_up = 0
    try:
        n_genes_down = pathway2n_down[pathway]
    except KeyError:
        n_genes_down = 0

    n_all = n_genes_up + n_genes_down   
    if n_all == 0:
        fraction_up = 0
        fraction_down = 0
        fraction_diff = 0
    else:
        fraction_up = float(n_genes_up)/float(n_genes)
        fraction_down = float(n_genes_down)/float(n_genes)
        fraction_diff = float(n_genes_up + n_genes_down)/float(n_genes)

    o.write(f"{pathway}\t{pathway_description}\t{n_genes}\t{n_genes_up}\t{n_genes_down}\t{n_all}\t{fraction_up}\t{fraction_down}\t{fraction_diff}\n")

    locus_list = locus_tag2pathway[locus_tag2pathway["pathway"] == pathway]["locus_tag"]

    for locus in locus_list:
        
        ko = locus_tag2ko[locus]
        ko_description = ko2description[ko]

        try:
            logFC  = count_table[count_table["locus_tag"] == locus]["logFC"][0]
            FDR  = count_table[count_table["locus_tag"] == locus]["FDR"][0]
            product  = ' '.join(str(count_table[count_table["locus_tag"] == locus]["product.1"][0]).split())
            signif  = count_table[count_table["locus_tag"] == locus]["signif"][0]
        except:
            # case very lowly expressed genes removed from the table
            logFC, FDR, product, signif = None, None, None, None
        q.write(f"{pathway}\t{pathway_description}\t{locus}\t{product}\t{ko}\t{ko_description}\t{logFC}\t{FDR}\t{signif}\n")

h = ["module",
     "module_description",
     "n_genes",
     "n_genes_up",
     "n_genes_down",
     "n_genes_up_and_down",
     "fraction_up",
     "fraction_down",
     "fraction_up_and_down"]

p.write("\t".join(h) + '\n')

for module in module2description:

    module_description = module2description[module]
    try:
        n_genes = module2n_genes[module]
    except KeyError:
        # module absent from the considered genome
        continue
    try:
        n_genes_up = module2n_up[module]
    except KeyError:
        n_genes_up = 0
    try:
        n_genes_down = module2n_down[module]
    except KeyError:
        n_genes_down = 0

    n_all = n_genes_up + n_genes_down   

    if n_all == 0:
        fraction_up = 0
        fraction_down = 0
        fraction_diff = 0
    else:
        fraction_up = float(n_genes_up)/float(n_genes)
        fraction_down = float(n_genes_down)/float(n_genes)
        fraction_diff = float(n_genes_up + n_genes_down)/float(n_genes)

    p.write(f"{module}\t{module_description}\t{n_genes}\t{n_genes_up}\t{n_genes_down}\t{n_all}\t{fraction_up}\t{fraction_down}\t{fraction_diff}\n")

    locus_list = locus_tag2module[locus_tag2module["module"] == module]["locus_tag"]

    for locus in locus_list:

        ko = locus_tag2ko[locus]
        ko_description = ko2description[ko]
        try:
            logFC  = count_table[count_table["locus_tag"] == locus]["logFC"][0]
            FDR  = count_table[count_table["locus_tag"] == locus]["FDR"][0]
            product  = ' '.join(str(count_table[count_table["locus_tag"] == locus]["product.1"][0]).split())
            signif  = count_table[count_table["locus_tag"] == locus]["signif"][0]
        except:
            logFC, FDR, product, signif = None, None, None, None

        r.write(f"{module}\t{module_description}\t{locus}\t{product}\t{ko}\t{ko_description}\t{logFC}\t{FDR}\t{signif}\n")




