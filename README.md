# compBio_miniProject
This repository contains Python code from the Comp Bio Mini Project. Project focuses on Human herpesvirus 5.

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
#### transcriptIdx.py:
This was used for development of the helper function used to retrieve and generate appropriate input for Kallisto.
#### sleuthInput.py:
This takes the output from Kallisto and creates a tab-delimited input file for sleuth. 
#### samtofastq.sh
A bash shell script that converts sam files from Bowtie2 to fastq files to be used in SPAdes



