
import re
htseq_files = snakemake.input["htseq_files"]
output_file = snakemake.output[0]


'''

['__no_feature', '978816']
['__ambiguous', '294']
['__too_low_aQual', '217114']
['__not_aligned', '18978990']
['__alignment_not_unique', '0']

'''


sample2data = {}
for htseq_file in htseq_files:
    sample = re.findall("samples/(.*)/rnaseq/htseq/counts.txt", htseq_file)[0]
    sample2data[sample] = {}
    sample2data[sample]["equal_0"] = 0
    sample2data[sample]["larger_than_0"] = 0
    sample2data[sample]["larger_than_10"] = 0
    sample2data[sample]["larger_than_100"] = 0
    sample2data[sample]["larger_than_1000"] = 0
    sample2data[sample]["larger_than_10000"] = 0
    sample2data[sample]["total"] = 0
    with open(htseq_file, 'r') as f:
        for row in f:
            data = row.rstrip().split("\t")
            if data[0] == '__no_feature':
                sample2data[sample]["__no_feature"] = data[1]
            elif data[0] == '__ambiguous':
                sample2data[sample]["__ambiguous"] = data[1]
            elif data[0] == '__too_low_aQual':
                sample2data[sample]["__too_low_aQual"] = data[1]
            elif data[0] == '__not_aligned':
                sample2data[sample]["__not_aligned"] = data[1]
            elif data[0] == '__alignment_not_unique':
                sample2data[sample]["__alignment_not_unique"] = data[1]
            else:
                if float(data[1]) == 0:
                    sample2data[sample]["equal_0"]+=1
                if float(data[1]) > 0:
                    sample2data[sample]["larger_than_0"]+=1
                if float(data[1]) > 10:
                    sample2data[sample]["larger_than_10"]+=1
                if float(data[1]) > 100:
                    sample2data[sample]["larger_than_100"]+=1
                if float(data[1]) > 1000:
                    sample2data[sample]["larger_than_1000"]+=1
                if float(data[1]) > 10000:
                    sample2data[sample]["larger_than_10000"]+=1
                sample2data[sample]["total"] += int(data[1])
with open(output_file, 'w') as o:
    for n, sample in enumerate(sample2data.keys()):
        if n == 0:
            o.write("sample\tnot_aligned\tno_feature\taligned\t0\t>0\t>10\t>100\t>1000\t>10000\\n")
        total_non_ambiguous_reads = int(sample2data[sample]["__not_aligned"]) + sample2data[sample]["total"]
        o.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (sample,
                                          sample2data[sample]["__not_aligned"],
                                          sample2data[sample]["__no_feature"],
                                          sample2data[sample]["total"],
                                          sample2data[sample]["equal_0"],
                                          sample2data[sample]["larger_than_0"],
                                          sample2data[sample]["larger_than_10"],
                                          sample2data[sample]["larger_than_100"],
                                          sample2data[sample]["larger_than_1000"],
                                          sample2data[sample]["larger_than_10000"]))
