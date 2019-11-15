

from Bio import Entrez
import ftplib
import logging
import os

Entrez.email = snakemake.params['NCBI_email']

ftp_login = ftplib.FTP(host='ftp.ncbi.nih.gov',
                       user='anonymous', 
                       passwd=snakemake.params['NCBI_email'])

handle1 = Entrez.esearch(db="assembly", term=snakemake.params['assembly_accession'])
record1 = Entrez.read(handle1)

ncbi_id = record1['IdList'][-1]

try:
    handle_assembly = Entrez.esummary(db="assembly", id=ncbi_id)
except:
    raise IOError('link to assembly could not be found')

assembly_record = Entrez.read(handle_assembly, validate=False)

try:
    ftp_path_refseq = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
except IndexError:
    raise IOError('refseq ftp path could not be found')

split_path=ftp_path_refseq.split('/')
redundant_name=split_path[-1]
path=ftp_path_refseq.split('ftp://ftp.ncbi.nlm.nih.gov/')[1]
ftp_login.cwd(path)
logging.info('Downloading from : {0}'.format(path))
logging.info('Downloading to {0}'.format(snakemake.output[0]))
ftp_login.retrbinary("RETR "+redundant_name+'_genomic.fna.gz',
                        open(os.path.join(os.getcwd(),snakemake.output[0]),"wb").write)
logging.info('Downloading to {0}'.format(snakemake.output[0]))
ftp_login.retrbinary("RETR "+redundant_name+'_genomic.gff.gz',
                        open(os.path.join(os.getcwd(),snakemake.output[1]),"wb").write)
logging.info('Done')







