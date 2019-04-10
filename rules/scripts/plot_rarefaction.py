
import re
htseq_files = snakemake.input["gene_counts"]
output_file = snakemake.output[0]


sample2data = {}
for htseq_file in htseq_files:
    m = re.findall("samples/(.*)/rnaseq/htseq/sample_([0-9]+)_counts.txt", htseq_file)[0]
    print(m)
    sample, n_reads = m
    if sample not in sample2data:
        sample2data[sample] = {}
    sample2data[sample][n_reads] = {}
    sample2data[sample][n_reads]["equal_0"] = 0
    sample2data[sample][n_reads]["larger_than_0"] = 0
    sample2data[sample][n_reads]["larger_than_10"] = 0
    sample2data[sample][n_reads]["larger_than_100"] = 0
    sample2data[sample][n_reads]["larger_than_1000"] = 0
    sample2data[sample][n_reads]["larger_than_10000"] = 0
    sample2data[sample][n_reads]["total"] = 0
    with open(htseq_file, 'r') as f:
        for row in f:
            data = row.rstrip().split("\t")
            if data[0] == '__no_feature':
                sample2data[sample][n_reads]["__no_feature"] = data[1]
            elif data[0] == '__ambiguous':
                sample2data[sample][n_reads]["__ambiguous"] = data[1]
            elif data[0] == '__too_low_aQual':
                sample2data[sample][n_reads]["__too_low_aQual"] = data[1]
            elif data[0] == '__not_aligned':
                sample2data[sample][n_reads]["__not_aligned"] = data[1]
            elif data[0] == '__alignment_not_unique':
                sample2data[sample][n_reads]["__alignment_not_unique"] = data[1]
            else:
                if float(data[1]) == 0:
                    sample2data[sample][n_reads]["equal_0"]+=1
                if float(data[1]) > 0:
                    sample2data[sample][n_reads]["larger_than_0"]+=1
                if float(data[1]) > 10:
                    sample2data[sample][n_reads]["larger_than_10"]+=1
                if float(data[1]) > 100:
                    sample2data[sample][n_reads]["larger_than_100"]+=1
                if float(data[1]) > 1000:
                    sample2data[sample][n_reads]["larger_than_1000"]+=1
                if float(data[1]) > 10000:
                    sample2data[sample][n_reads]["larger_than_10000"]+=1
                sample2data[sample][n_reads]["total"] += int(data[1])

with open(output_file, 'w') as o:
    for sample in sample2data:
        for n_reads in sample2data[sample]:
                    o.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (sample,n_reads,
                                                      sample2data[sample][n_reads]["__not_aligned"],
                                                      sample2data[sample][n_reads]["__no_feature"],
                                                      sample2data[sample][n_reads]["total"],
                                                      sample2data[sample][n_reads]["equal_0"],
                                                      sample2data[sample][n_reads]["larger_than_0"],
                                                      sample2data[sample][n_reads]["larger_than_10"],
                                                      sample2data[sample][n_reads]["larger_than_100"],
                                                      sample2data[sample][n_reads]["larger_than_1000"],
                                                      sample2data[sample][n_reads]["larger_than_10000"]))
