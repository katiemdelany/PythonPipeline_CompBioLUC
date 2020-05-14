# Computational Biology 483 Loyola University Chicago
This repository contains Python code from an LUC Computational Biology project. Project focuses on Human herpesvirus 5.

### Software Requirements
* Linux/Unix
* Python3
* Biopython
* Kallisto
* Bowtie2
* SPAdes

### Main Python Script
#### miniMain.py
##### To run this script, please clone this github repository with:
```
git clone https://github.com/katiemdelany/compBio_miniProject.git
```
The repository contains a test dataset which are a subset of 30,000 lines from the input fastq files from InptFiles() function. The SRR IDs for the test data are: SRR5660030 SRR5660033 SRR5660044 SRR5660045 

##### To run the script:
Argument format: SRR run ID numbers
```
<SRR Run ID> <SRR Run ID> <SRR Run ID> <SRR Run ID>
```
##### Example:
Use this command in the command line to run the script. 
```
python3 miniMain.py SRR5660030 SRR5660033 SRR5660044 SRR5660045
```

### Files Included in Repo:
#### miniMain.py: 
miniMain.py runs the whole pipeline and contains helper function to complete each step. It outputs files to multiple folders and writes into a .log file to record key results.
##### Key output of miniMain.py:
* miniproject.log - contains key results of pipeline
* EF999921.fasta - HCMV genome
* EF999921_CDS.fasta - coding sequences of HCMV genome
* LargeContig.txt - all contigs over 1000 bp
* Assemble.fasta - assembly file of LargeContig.txt

#### transcriptIdx.py:
This was used for development of the helper function used to retrieve and generate appropriate input for Kallisto.
#### sleuthInput.py:
This takes the output from Kallisto and creates a tab-delimited input file for sleuth. 
#### samtofastq.sh
A bash shell script that converts sam files from Bowtie2 to fastq files to be used in SPAdes

### Author
Kathleen M. Delany

#### Assignment Details
1. We would like to compare HCMV transcriptomes 2- and 6-days post-infection (dpi). First, retrieve the following transcriptomes from two patient donors from SRA and convert to paired-end fastq files. You can use wget (by constructing the path based on the SRR numbers for each of these samples). Donor 1 (2dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896360 Donor 1 (6dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896363 Donor 3 (2dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896374 Donor 3 (6dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896375
2. We will quantify TPM in each sample using kallisto, but first we need to build a transcriptome index for HCMV (NCBI accession EF999921). Use Biopython to retrieve and generate the appropriate input and then build the index with kallisto (https://pachterlab.github.io/kallisto/). You will need to extract the CDS features from the GenBank format. Write the following to your log file (replace # with the number of coding sequences in the HCMV genome): 2 The HCMV genome (EF99921) has # CDS.
3. Quantify the TPM of each CDS in each transcriptome using kallisto and use these results as input to find differentially expressed genes between the two timepoints (2pi and 6dpi) using the R package sleuth (https://pachterlab.github.io/sleuth/about). Write the following details for each significant transcript (FDR < 0.05) to your log file, include a header row, and tab-delimit each item: target_id test_stat pval qval
4. It has been proposed that HCMV disease and pathogenesis may be related to the genetic diversity of the virus (Renzette et al. https://www.ncbi.nlm.nih.gov/pubmed/25154343/). Which publicly available strains are most similar to these patient samples? To compare to other strains, we will assemble these transcriptome reads. We don’t expect assembly to produce the entire genome, but enough to be useful in BLAST. Virus sequencing experiments often include host DNAs. It is difficult to isolate the RNA of just the virus (as it only transcribes during infection of the host cell). Before assembly, let’s make sure our reads map to HCMV. Using Bowtie2, create an index for HCMV (NCBI accession EF999921). Next, save the reads that map to the HCMV index for use in assembly. Write to your log file the number of reads in each transcriptome before and after the Bowtie2 mapping. For instance, if I was looking at the Donor 1 (2dpi) sample, I would write to the log (numbers here are arbitrary): Donor 1 (2dpi) had 230000 read pairs before Bowtie2 filtering and 100000 read pairs after.
5. Using the Bowtie2 output reads, assemble all four transcriptomes together to produce 1 assembly via SPAdes. Write the SPAdes command used to the log file.
6. Write Python code to calculate the number of contigs with a length > 1000 and write the # out to the log file as follows: There are # contigs > 1000 bp in the assembly. For steps 7-9, you will only consider those contigs > 1000 bp in length.
7. Write Python code to calculate the length of the assembly (the total number of bp in all of the contigs > 1000 bp in length) and write this # out to the log file as follows: There are # bp in the assembly.
8. Write Python code to concatenate all of the contigs > 1000 bp in length into 1 fasta sequence. Separate each contig by a stretch of 50 N’s.
9. Using this concatenated fasta file, blast via NCBIWWW.qblast to query the nr database limited to members of the Herpesviridae family to identify the top 10 hits. For the top 10 hits, write the following to your log file: sequence title, alignment length, number of HSPs, and for the top HSP: HSP identities, HSP gaps, HSP bits, and HSP expect scores. Include the following header row and tab-delimit each item: seq_title align_len number_HSPs topHSP_ident topHSP_gaps topHSP_bits topHSP_expect

