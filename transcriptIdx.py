from Bio import SeqIO
from Bio import Entrez
from Bio import SearchIO
from Bio import GenBank
from Bio.Seq import Seq
import os
os.chdir("C:\\Users\\katie\\Downloads\\")



def getTranscIndex():
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
    print(count)

getTranscIndex()
