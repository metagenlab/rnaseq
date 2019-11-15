


import pandas
import itertools

include: "rules/making_sample_dataset.rules"
include: "rules/logging.rules"
include: "rules/fastqc.rules"
include: "rules/qualimap.rules"
include: "rules/trimmomatic.rules"
include: "rules/prepare_fastq.rules"
include: "rules/bwa.rules"
include: "rules/counts.rules"
include: "rules/reads_sample.rules"
include: "rules/edgeR.rules"
include: "rules/reference_genome.rules"

rule rarefaction:
    input:
        expand("samples/{sample}/reads/trimmed/fastqc/single_fastqc.html", sample=read_naming.keys()),
        "report/rnaseq/htseq/summary_all.tab",
        expand("report/qualimap/{sample}/bwa/GCA_000092785/genome_results.txt", sample=read_naming.keys()),
        "aggregated/rarefaction.pdf",

def condition_combination(wildcards):
    file_list = []
    condition_list = list(set(all_samples["condition"]))
    for pair in itertools.combinations(condition_list, 2):
        cons1, cond2 = pair
        file_list.append(f'report/rnaseq/edger/{cons1}_vs_{cond2}/plotMDS.pdf')
    print(file_list)
    return file_list


rule edger:
    input:
        condition_combination