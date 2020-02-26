import os
import logging
import argparse
from Bio import SeqIO
from Bio import Entrez
from Bio import SearchIO
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
curr =  os.getcwd()
os.chdir(curr)


def arg_parser():
    parser = argparse.ArgumentParser(description='Evaluate sra files')
    parser.add_argument('srrfiles', metavar= '.srr', type=str, nargs = '*', help = 'SRRs')

    return parser.parse_args()




def InptFiles(SRR):
    """
    This will retrieve the transcriptomes of given SRR numbers and convert them to
    paired-end fastq files. Files from wget download as SRR#######.1 get renamed to SRR#######  

    """
    getFiles = 'wget https://sra-download.st-va.ncbi.nlm.nih.gov/sos1/sra-pub-run-12/SRR' + SRR + '/' + SRR + '.1'
    renameFile = 'mv '+ str(SRR)+'.1 ' + SRR
    splitFiles = 'fastq-dump -I --split-files '+ str(SRR)
    os.system(getFiles)
    os.system(renameFile)
    os.system(splitFiles)


def getTranscriptomeIndex():
    outFasta = open("EF999921.fasta", "w")
    outFile = open("EF999921_CDS.fasta","w")
    Entrez.email = 'kdelany@luc.edu'
    #retrieve record
    handle = Entrez.efetch(db = 'nucleotide', id= 'EF999921', rettype= 'fasta')
    records = list(SeqIO.parse(handle, "fasta"))
    #use .seq object to write out fasta file of whole sequence
    outFasta.write('>' + str(records[0].description)+ '\n' + str(records[0].seq))
    outFasta.close()
    #Fetch genbank format object
    GBhandle = Entrez.efetch(db = 'nucleotide', id = 'EF999921', rettype= 'gb', retmode='text')
    count = 0
    #this writes out the CDS sequences with a >identifier to a new file 
    #also adds a count for recording number of CDS sequences
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



def SleuthInput(SRRs):
    """
    Create file for input into sleuth
    
    """
    #input file for sleuth
    covFile = open('cov.txt','w')
    condition1 = "2dpi"
    condition2 = "6dpi"
    #initial line in file
    covFile.write('sample'+ '\t' + 'condition' + '\t' + 'path' + '\n')
    #based on SRR number, write condition and path to input file
    for i in SRRs:
        path = '/data/kdelany/compBio_miniProject/'+i
        if int(i[3:])%2==0:
            covFile.write(str(i)+ '\t' + condition1 + '\t'+ str(path)+ '\n')
        else:
            covFile.write(str(i)+ '\t' + condition2 + '\t'+ str(path)+ '\n')
    covFile.close()


def Sleuth():
    """
    Runs sleuth in R
    Reads sleuth output and adds to log file

    """
    runSleuth = 'Rscript sleuth.R'
    os.system(runSleuth)
    output = open('topten.txt','r')
    listoflines= output.readlines()
    for line in listoflines:
        logging.info(line)



def bowtie2build(SRR):
    """ Builds a Bowtie index for HCMV """
    build_cmd = 'bowtie2-build ./EF999921_CDS.fasta EF999921'
    os.system(build_cmd)
    bowtie_cmd = 'bowtie2 --quiet --no-unal --al-conc BOW_'+SRR+'.fastq -x EF999921 -1 '+SRR+ '_1.fastq -2'+SRR+'_2.fastq -S '+SRR+ '.sam'
    os.system(bowtie_cmd)



def Sam2Fastq(SRR):
    """ Bash script to convert sam files to fastq files """
    runConvert = 'bash samtofastq.sh ' + str(SRR)
    os.system(runConvert)



def getNumReads(SRR):
    """ Returns the number of reads before bowtie2 and after """
    name = ''
    #tells which donor to write into log file
    if SRR == 'SRR5660030':
        name = 'Donor 1 (2dpi)'
    elif SRR == 'SRR5660033':
        name = 'Donor 1 (6dpi)'
    elif SRR == 'SRR5660044':
        name = 'Donor 3 (2dpi)'
    else:
        name = 'Donor 3 (6dpi)'
    #before running Bowtie2
    SRRfile = open(str(SRR)+'_1.fastq')
    SRRfile1 = open(str(SRR)+'_2.fastq')
    count1 = 0
    count = 0
    #count reads in SRR file before bowtie
    for line in SRRfile:
        count1+=1
    for line in SRRfile1:
        count+=1               
    beforeCount =(count+ count1)/4
    AfterFile1 = open('BOW_'+SRR+'.1.fastq')
    AfterFile2 = open('BOW_'+SRR+'.2.fastq')
    #count reads after bowtie mapping
    count2 = 0
    count3 = 0
    for line in AfterFile1:
        count2 +=1
    for line in AfterFile2:
        count3 +=1
    afterCount =(count3+ count2)/4

    logging.info(str(name)+' had ' +str(beforeCount) + ' read pairs before Bowtie2 filtering and '+ str(afterCount)+' pairs after.')



def SPAdes(SRRs):
    """
    Sets SRRs to separate variables
    """
    SRR1= SRRs[0]
    SRR2= SRRs[1]
    SRR3= SRRs[2]
    SRR4= SRRs[3]
   # spades_cmd = 'spades -k 55,77,99,127 -t 2 --only-assembler -s BOW_'+SRR1+'.1.fastq -s '+SRR2+'_bow.fastq -s '+SRR3+ '_bow.fastq -s '+SRR4 + '_bow.fastq -o SpadesAssembly/'
    spades_cmd ='spades -k 55,77,99,127 -t 2 --only-assembler --pe1-1 BOW_'+SRR1+'.1.fastq --pe1-2 BOW_'+SRR1+'.2.fastq --pe2-1 BOW_'+SRR2+'.1.fastq --pe2-2 BOW_'+SRR2+'.2.fastq --pe3-1 BOW_'+SRR3+'.1.fastq --pe3-2 BOW_'+SRR3+'.2.fastq --pe4-1 BOW_'+SRR4+'.1.fastq --pe4-2 BOW_'+SRR4+'.2.fastq -o SpadesAssembly/'
    os.system(spades_cmd)
    logging.info(str(spades_cmd))



def numContigs():
    """
    Counts number of contigs greater than 100

    """
    newFile = open('LargeContigs.txt', 'w')
    count = 0
    #initialize count to 0 and parse SPAdes output as fasta
    handle = SeqIO.parse('./SpadesAssembly/contigs.fasta','fasta')
    #if the sequence len  is greater than 1000 add to count
    #Add sequences greater than 1000 to outfile
    for record in handle:
        m = len(record.seq)
        if m > 1000:
            count +=1
            newFile.write('> '+str(record.id) + '\n' + str(record.seq) + '\n')
    newFile.close()

    logging.info('There are '+str(count)+' contigs > 1000 bp in the assembly.')



def countContigs():
    """
    This sums the total amount of base pairs in the assembly

    """
    newFile = open('LargeContigs.txt', 'r')
    #parse file with the contigs > 1000 as a fasta
    handle = SeqIO.parse('LargeContigs.txt', 'fasta')
    lenList = []
    #add each sequence len to a list
    for record in handle:
        m = len(record.seq)
        lenList.append(int(m))
    #sum the list
    total = sum(lenList) 
   
    logging.info('There are '+str(total)+' bp in the assembly.')


def assembleContigs():
    """
    This concatenates the contigs > 1000 into one fasta seq separated by 50 Ns

    """
    assemblyFile = open('Assemble.fasta','w')
    inFile = open('LargeContigs.txt','r')
    handle = SeqIO.parse(inFile, 'fasta')
    concat = ''
    #parse large contigs file as fasta, build a string of record.seq objects
    #separate by 50 N's
    for record in handle:
        seq = str(record.seq)
        concat+= seq+ 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'
    assemblyFile.write(concat)
    assemblyFile.close()


def blast():
    inFile = open('Assemble.fasta').read()
    #blast command with added entrez_query to limit search
    blast1 = NCBIWWW.qblast('blastn','nr',inFile, entrez_query='"Herpesviridae"[organism]')
    with open("my_blast.xml", "w") as outhandle:
        outhandle.write(blast1.read())
    blast1.close()
    outhandle.close()
    #parse output
    qresult = SearchIO.read('my_blast.xml','blast-xml')
    #adds first line to log file of desired information
    logging.info('seq_title\talign_len\tnumber_HSPs\ttopHSP_ident\ttopHSP_gaps\ttopHSP_bits\ttopHSP_expect')
    #length of results
    m = len(qresult)
    #num for top ten (limit to first 10 results)
    top_ten = 9
    #was not sure if there would be 10 results, so if there arent, set topten to m
    if top_ten > m:
        top_ten = m
    #for the first 10 resutls
    #store hit and hsp information as variables. res would be the hit. 
    #add variables to log separated by tab
    for i in range(0, top_ten):
        res = qresult[i]
        hsp1 = qresult[i][0]
        seq_title = str(res.description)
        align_len = str(res.seq_len)
        number_HSPs = str(len(res.hsps))
        topHSP_ident = str(hsp1.ident_num)
        topHSP_gaps = str(hsp1.gap_num)
        topHSP_bits = str(hsp1.bitscore)
        topHSP_expect = str(hsp1.evalue)
        logging.info(seq_title + '\t' + align_len+ '\t'+ number_HSPs+ '\t'+topHSP_ident+ '\t'+ topHSP_gaps+ '\t'+ topHSP_bits+ '\t'+ topHSP_expect)



def main():
    """ Takes in SRR id number arguments in command line and runs through"""
    
    args = arg_parser()
    
    logging.info('SRA values: %s', args.srrfiles)
    
    print(args.srrfiles)
  
    #This is where I download and split the SRR data files. Commented out to save time.
    #Use test data from cloned git hub
#    for i in args.srrfiles:
#       InptFiles(i)

    result =  getTranscriptomeIndex()
    logging.info('The HCMV genome (EF999921) has '+str(result)+' CDS.')


    for i in args.srrfiles:
        Kallisto(i)

    SleuthInput(args.srrfiles)
    Sleuth()

    
    for i in args.srrfiles:
        bowtie2build(i)
       # Sam2Fastq(i)
        getNumReads(i)
        
    SPAdes(args.srrfiles)
    numContigs()
    countContigs()
    assembleContigs()
    blast()
    

if __name__ == '__main__':
    logger = logging.getLogger(__name__)
    logFormatter = '%(message)s'
    logging.basicConfig(filename='miniProject.log', format=logFormatter, level=logging.DEBUG)
    main()
