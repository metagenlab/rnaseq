
import pandas
import re

locus_tag2pathway = pandas.read_csv(snakemake.input["locus_tag2pathway"], sep="\t", names=["locus_tag","pathway"])

pathway2n_genes = locus_tag2pathway.groupby("pathway").size().to_dict()

pathway2description = pandas.read_csv(snakemake.input["pathway2description"], sep="\t", names=["pathway","pathway_description"]).set_index("pathway").to_dict()["pathway_description"]

# rpkm_table 
rpkm_table = pandas.read_csv(snakemake.input["rpkm_table"], sep="\t", header=0).set_index('locus_tag')

sample_table = pandas.read_csv(snakemake.params["sample_table"], sep="\t", header=0).set_index("SampleName")

# extract condition/sample correspondance
condition2samples = {}
for n, row in sample_table.iterrows():
    n = re.sub("-", ".", n)
    if row["condition"] not in condition2samples:
        condition2samples[row["condition"]] = [n] 
    else:
        condition2samples[row["condition"]].append(n)

cond1, cond2 = snakemake.input["count_file"].split("/")[3].split("_vs_")

# count table
count_table = pandas.read_csv(snakemake.input["count_file"], sep="\t", index_col=0)

df_pathway = count_table.join(locus_tag2pathway.set_index('locus_tag'), on='locus_tag').join(rpkm_table, on='locus_tag')

rpkm_column_list = condition2samples[cond1] + condition2samples[cond2]

column_list = ['pathway', 'locus_tag'] + rpkm_column_list

pathway2rpkm_sum = df_pathway[column_list].groupby(["pathway"]).agg('sum')

median_rpkm = []

for pathway, row in pathway2rpkm_sum.iterrows():
    n_genes = pathway2n_genes[pathway]
    description = pathway2description[pathway]
    label = f"{pathway}/{description}"
    m = [str(i/float(n_genes)) for i in row[rpkm_column_list]]
    median_rpkm.append([label] + m)

with open(snakemake.output[0], 'w') as f:
    h = ["pathway"] + rpkm_column_list
    f.write("\t".join(h) + '\n')
    for row in median_rpkm:
        f.write("\t".join(row) + '\n')

with open(snakemake.output[1], 'w') as f:
    h = ["pathway"] + rpkm_column_list
    f.write("\t".join(h) + '\n')
    for pathway, row in pathway2rpkm_sum.iterrows():
        r = [pathway] + [str(i) for i in row[rpkm_column_list]]
        f.write("\t".join(r) + '\n')

