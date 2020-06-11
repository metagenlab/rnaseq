import sqlite3
from pronto import Ontology
from Bio.SeqUtils import CheckSum
from Bio import SeqIO

faa_path = snakemake.input["faa_file"]
goa_path = snakemake.input["goa"]
go = Ontology(snakemake.input["go_obo"])
go_annotations = open(snakemake.output["go_annotations"], 'w')
uniparcdb = snakemake.input["uniparcdb"]

conn = sqlite3.connect(goa_path)
cursor = conn.cursor()

sqlatt = f'attach database "{uniparcdb}" as uniparc;'
cursor.execute(sqlatt,)

# 1. retrieve uniprot accession from exact match (hash)
# 2. retrieve GO annotations

for record in SeqIO.parse(faa_path, "fasta"):
    
    checksum = CheckSum.seguid(record.seq)
    

      
    sqlq = 'select * from uniparc.uniparc_accession where sequence_hash="%s"' % checksum
    
    uniparc_id = cursor.execute(sqlq,).fetchall()[0][0]
    
    print("uid", uniparc_id)
    
    sqlq2 = 'select distinct accession from uniparc_cross_references t1 ' \
      ' inner join crossref_databases t2 on t1.db_id=t2.db_id ' \
      ' where t1.uniparc_id=%s and db_name in ("UniProtKB/Swiss-Prot", "UniProtKB/TrEMBL");' % uniparc_id
    print(sqlq2)
    uniprotkb_acc_list = [i[0].split(".")[0] for i in cursor.execute(sqlq2,).fetchall()]
    print("hits:", uniprotkb_acc_list)
    
    sql = 'select distinct GO_id from goa_table where uniprotkb_accession in ("%s");' % ('","'.join(uniprotkb_acc_list))
    go_terms_data = cursor.execute(sql,).fetchall()
    
    go_list = ';'.join([i[0] for i in go_terms_data])
    
    go_annotations.write(f"{record.id}\t{go_list}\n")

go_annotations.close()

