#!/bin/bash

#this takes each sam file and converts it to fastq
grep -v ^@ ${1}.sam | awk '{print "@"$1"\n"$10"\n+\n"$11}' > $1_bow.fastq


