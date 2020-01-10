
import pandas
from Bio import SeqIO

reference_gbk = snakemake.input[0] # reference_genome.gff
reference_genome_uniprotKB_ko = snakemake.input[1] # reference/reference_genome_uniprotKB_ko.txt
ko_description = snakemake.input[2] # reference/ko_description.txt
reference_genome_uniprotKB_COG = snakemake.input[3] # reference/reference_genome_uniprotKB_COG.txt
COG_annotations = snakemake.input[4] # reference/COG_annotations.tsv"
o = open(snakemake.output[0], 'w')

header = ["locus_tag",	
          "gene",	
          "product",
          "COG",
          "COG_description",
          "KO",
          "KO_description",
          "orthogroup",
          "orthogroup_freq"]

locus_tag2ko = pandas.read_csv(reference_genome_uniprotKB_ko,
                           delimiter='\t',
                           names=["locus_tag","ko"]).set_index("locus_tag").to_dict()["ko"]

locus_tag2cog = pandas.read_csv(reference_genome_uniprotKB_COG,
                           delimiter='\t',
                           names=["locus_tag","cog"]).set_index("locus_tag").to_dict()["cog"]

ko2description = pandas.read_csv(ko_description,
                           delimiter='\t',
                           names=["ko","description"]).set_index("ko").to_dict()["description"]
for ko in list(ko2description.keys()):
    ko_edit = ko.split(":")[1]
    if ko != ko_edit:
        ko2description[ko_edit] = ko2description[ko]

cog2description = pandas.read_csv(COG_annotations,
                           delimiter='\t',
                           names=["cog","category","description"]).set_index("cog").to_dict()["description"]

o.write("\t".join(header)+'\n')

for rec in SeqIO.parse(reference_gbk, "genbank"):
    for feature in rec.features:
        if feature.type == 'CDS':
            locus_tag = feature.qualifiers["locus_tag"][0]
            product = feature.qualifiers["product"][0]
            gene = feature.qualifiers["gene"][0]
            try:
                ko = locus_tag2ko[locus_tag]
                ko_description = ko2description[ko]
            except KeyError:
                ko = '-'
                ko_description = '-'
            try:
                cog = locus_tag2cog[locus_tag]
                cog_description = cog2description[cog]
            except:
                cog = ''
                cog_description = '-'
            o.write(f'{locus_tag}\t{gene}\t{product}\t{cog_description}\t{ko}\t{ko_description}\t-\t-\n')
        if feature.type in ["tRNA", "rRNA", "ncRNA"]:
            locus_tag = feature.qualifiers["locus_tag"][0]
            product = feature.qualifiers["product"][0]
            gene = feature.qualifiers["gene"][0]
            o.write(f'{locus_tag}\t{gene}\t{product}\t-\t-\t-\t-\t-\n')



        
