library(sleuth)
library(data.table)
setwd("/data/kdelany/compBio_miniProject")

##TABLE
help(read.table)
stab <- read.table("cov.txt", header=TRUE, stringsAsFactors = FALSE, sep='\t')
View(stab)

#sleuth object
so <- sleuth_prep(stab)

#differential expression analysis
#fit a model comparing conditions
so <- sleuth_fit(so, ~condition, 'full')
so
#reduced model
so <- sleuth_fit(so, ~1, 'reduced')
so
#likelihood ratio test for differential expression
so <- sleuth_lrt(so, 'reduced','full')
so

#extract results
library(dplyr)
sleuth_table <- sleuth_results(so, 'reduced:full','lrt',show_all=FALSE)

#filtering significant results
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05) %>% dplyr::arrange(pval)
sig_sleuth <- sleuth_significant %>% select(target_id, test_stat, pval, qval)
sig_sleuth
     
#write target id, test stat, pval and qval for significant transcript
#include header, tab-delimit
write.table(sig_sleuth, file="topten.txt",quote= FALSE,row.names= FALSE)
     