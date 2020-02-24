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

def getTranscriptomeIndex():
    outFasta = open("EF999921.fasta", "w")
    outFile = open("EF999921_CDS.fasta","w")
    Entrez.email = 'kdelany@luc.edu'
    handle = Entrez.efetch(db = 'nucleotide', id= 'EF999921', rettype= 'fasta')
    records = list(SeqIO.parse(handle, "fasta"))
    outFasta.write('>' + str(records[0].description)+ '\n' + str(records[0].seq))
    outFasta.close()
    #SeqIO.write(records[0], 'EF999921.fasta', 'fasta')
    GBhandle = Entrez.efetch(db = 'nucleotide', id = 'EF999921', rettype= 'gb', retmode='text')
    count = 0
    for record in SeqIO.parse(GBhandle, 'genbank'):
            for feature in record.features:
                if feature.type == "CDS":
                    count +=1
                    outFile.write('>' + str(feature.qualifiers['protein_id']).replace('[', '').replace(']', '').replace("'","") + '\n' + str(feature.location.extract(record).seq) +'\n')                
    outFile.close()
    return(count)


def Kallisto(SRR):
    """
    Takes in SRR number and creates/runs the Kallisto command line commands.
    Uses CDS fasta file created in the transcriptome index

    """
    kallisto_cmd = 'time kallisto index -i HCMVindex.idx EF999921_CDS.fasta'
    os.system(kallisto_cmd)
    kallisto_run = 'time kallisto quant -i HCMVindex.idx -o ./' + str(SRR) +' -b 30 -t 4 '+ str(SRR) + '_1.fastq '+ str(SRR)+ '_2.fastq'
    os.system(kallisto_run)


def main():
    """ Takes in SRR id number arguments in command line and runs through"""
    
    parser = argparse.ArgumentParser(description= 'Process SRR (RunID) numbers')
    parser.add_argument('SRR', metavar= 'N', type=str, nargs = '+', help = 'SRRs')
    parser.add_argument('--downloadfiles', help='Retrieve SRA from database and convert to paired-end fastq files')
    args = parser.parse_args()

    
    
#    if args.downloadfiles:
#        for i in args.SRR:
#            InptFiles(i)
#    with open('miniproject.log', 'a') as f_out:
#        f_out.write('SRA files download done.')
#        f_out.close()
    result =  getTranscriptomeIndex()
    with open('miniproject.log', 'a') as f_out:
        f_out.write('The HCMV genome (EF999921) has ' + str(result) + ' CDS.')
        f_out.close()
    
#    for i in args.SRR:
#        Kallisto(i)
    
if __name__ == '__main__':
    main()
