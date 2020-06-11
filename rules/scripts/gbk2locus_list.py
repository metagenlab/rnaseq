

from Bio import SeqIO
records = SeqIO.parse(snakemake.input[0], "genbank")
with open(snakemake.output[0], "w") as f:
    for record in records:
        for feature in record.features:
            if "locus_tag" in feature.qualifiers:
                f.write("%s\n" % feature.qualifiers["locus_tag"][0])