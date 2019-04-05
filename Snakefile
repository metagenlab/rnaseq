

# fastqc
# timming


import pandas

include: "rules/making_sample_dataset.rules"
include: "rules/logging.rules"
include: "rules/fastqc.rules"
include: "rules/qualimap.rules"
include: "rules/trimmomatic.rules"
include: "rules/prepare_fastq.rules"
include: "rules/bwa.rules"
include: "rules/counts.rules"
include: "rules/reads_sample.rules"

rule all:
    input:
        expand("samples/{sample}/reads/trimmed/fastqc/single_fastqc.html", sample=read_naming.keys()),
        "report/rnaseq/htseq/summary_all.tab",
        expand("report/qualimap/{sample}/bwa/GCA_000092785/genome_results.txt", sample=read_naming.keys()),
        expand("samples/{sample}/mapping/bwa/GCA_000092785_mapped_sample_1000000.fastq", sample=read_naming.keys()),
        expand("aggregated/{sample}.txt", sample=read_naming.keys()),
