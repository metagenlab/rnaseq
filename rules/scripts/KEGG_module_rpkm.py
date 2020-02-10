
import pandas
import re

# b2388	M00001
locus_tag2module = pandas.read_csv(snakemake.input["locus_tag2module"], sep="\t", names=["locus_tag","module"])

module2n_genes = locus_tag2module.groupby("module").size().to_dict()

# M00001	Glycolysis (Embden-Meyerhof pathway), glucose => pyruvate
module2description = pandas.read_csv(snakemake.input["module2description"], sep="\t", names=["module","module_description"]).set_index("module").to_dict()["module_description"]

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

df_module = count_table.join(locus_tag2module.set_index('locus_tag'), on='locus_tag').join(rpkm_table, on='locus_tag')

rpkm_column_list = condition2samples[cond1] + condition2samples[cond2]

column_list = ['module', 'locus_tag'] + rpkm_column_list

module2rpkm_sum = df_module[column_list].groupby(["module"]).agg('sum')

median_rpkm = []

for module, row in module2rpkm_sum.iterrows():
    n_genes = module2n_genes[module]
    description = module2description[module]
    label = f"{module}/{description}"
    m = [str(i/float(n_genes)) for i in row[rpkm_column_list]]
    median_rpkm.append([label] + m)

with open(snakemake.output[0], 'w') as f:
    h = ["module"] + rpkm_column_list
    f.write("\t".join(h) + '\n')
    for row in median_rpkm:
        f.write("\t".join(row) + '\n')
