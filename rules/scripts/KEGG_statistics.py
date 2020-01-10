
import pandas
import re
# 
count_file = snakemake.input["count_file"]

# eco:b0114	path:eco00010
locus_tag2pathway = snakemake.input["locus_tag2pathway"]

# eco:b2388	md:eco_M00001
locus_tag2module = snakemake.input["locus_tag2module"]

# md:M00001	Glycolysis (Embden-Meyerhof pathway), glucose => pyruvate
module2description = snakemake.input["module2description"]

# path:map00010	Glycolysis / Gluconeogenesis
pathway2description = snakemake.input["pathway2description"]

def parse_simple(file_name):
    col1_2_col2 = {}
    with open(file_name, "r") as f:
        for row in f:
            data = row.rstrip().split("\t")
            data[0] = data[0].split(":")[1]
            col1_2_col2[data[0]] = data[1]
    return col1_2_col2
    

def parse_kegg_table(file_name, split_string=False, ):
    col1_to_col2 = {}
    col2_to_col1 = {}
    with open(file_name, "r") as f:
        for row in f:
            data = row.rstrip().split("\t")
            data[0] = data[0].split(':')[1]
            if split_string:
                data[1] = data[1].split(split_string)[1]
            if data[0] not in col1_to_col2:
                col1_to_col2[data[0]] = [data[1]]
            else:
                col1_to_col2[data[0]].append(data[1])

            if data[1] not in col2_to_col1:
                col2_to_col1[data[1]] = [data[0]]
            else:
                col2_to_col1[data[1]].append(data[0])
    return col1_to_col2, col2_to_col1

locus_tag2pathway_list, pathway2locus_list = parse_kegg_table(locus_tag2pathway, split_string=":")        
locus_tag2module_list, module2locus_list = parse_kegg_table(locus_tag2module, split_string="_") 
modle2description_dico = parse_simple(module2description)
pathway2description_dico = parse_simple(pathway2description)

pathway2n_significative_change = {}
module2n_significative_change = {}

count_table = pandas.read_csv(count_file, sep="\t", index_col=0)

count_table.head()

for n, row in count_table.iterrows():
    locus = row["locus_tag"]
    if row["signif"] is True:
        if locus in locus_tag2module_list:
            for module in locus_tag2module_list[locus]:
                if module not in module2n_significative_change:
                    module2n_significative_change[module] = 1
                else:
                    module2n_significative_change[module] += 1
        if locus in locus_tag2pathway_list:
            for pathway in locus_tag2pathway_list[locus]:
                if pathway not in pathway2n_significative_change:
                    pathway2n_significative_change[pathway] = 1
                else:
                    pathway2n_significative_change[pathway] += 1

o = open(snakemake.output[0], "w")
p = open(snakemake.output[1], "w")

print("module2n_significative_change", module2n_significative_change)
print("pathway2n_significative_change", pathway2n_significative_change)

o.write("pathway\tpathway_description\tn_genes\tn_genes_diff\tfraction\n")
for pathway in pathway2n_significative_change:
    pathway_id = re.search("\d+", pathway)[0]
    path_edit = "map%s" % pathway_id
    pathway_description = pathway2description_dico[path_edit]
    n_genes = len(pathway2locus_list[pathway])
    n_genes_diff = pathway2n_significative_change[pathway]
    fraction = float(n_genes_diff)/float(n_genes)
    o.write("%s\t%s\t%s\t%s\t%s\n" % (path_edit,
                                  pathway_description,
                                  n_genes,
                                  n_genes_diff,
                                  fraction))

p.write("module\tmodule_description\tn_genes\tn_genes_diff\tfraction\n")
for module in module2n_significative_change:
    module_description = modle2description_dico[module]
    n_genes = len(module2locus_list[module])
    n_genes_diff = module2n_significative_change[module]
    fraction = float(n_genes_diff)/float(n_genes)
    p.write("%s\t%s\t%s\t%s\t%s\n" % (module,
                                  module_description,
                                  n_genes,
                                  n_genes_diff,
                                  fraction))