import os
import argparse
from Bio import SeqIO
from Bio import Entrez
from Bio import SearchIO
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
os.getcwd()

def InptFiles(SRR):
    """
    This will retrieve the transcriptomes of given SRR numbers and convert them to
    paired-end fastq files. Files from wget download as SRR#######.1 get renamed to SRR#######  

    """
    #getFiles = 'wget https://sra-download.st-va.ncbi.nlm.nih.gov/sos1/sra-pub-run-12/SRR' + SRR + '/' + SRR + '.1'
    #wget https://sra-download.st-va.ncbi.nlm.nih.gov/sos1/sra-pub-run/12/SRR5660030/SRR5660030.1
    renameFile = 'mv '+ str(SRR)+'.1 ' + SRR
    splitFiles = 'fastq-dump -I --split-files '+ str(SRR)
    #os.system(getFiles)
    os.system(renameFile)
    os.system(splitFiles)


def main():
    """ Takes in SRR id number arguments in command line and runs through"""
    
    parser = argparse.ArgumentParser(description= 'Process SRR (RunID) numbers')
    parser.add_argument('SRR', metavar= 'N', type=str, nargs = '+', help = 'SRRs')
    parser.add_argument('--downloadfiles', help='Retrieve SRA from database and convert to paired-end fastq files')
    args = parser.parse_args()

    
    
    if args.downloadfiles:
        for i in args.SRR:
            InptFiles(i)
    with open('miniproject.log', 'a') as f_out:
        f_out.write('SRA files download done.')
        f_out.close()


    
if __name__ == '__main__':
    main()
