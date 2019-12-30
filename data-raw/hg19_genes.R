##########################################
# Goal: Generate the GRanges object for the hg19 varbins
# containing gene symbols
#
# # Output: hg19_genes object used in sysdata.rda
#
##########################################
# Author: Darlan Conterno Minussi
##########################################

# using package Homo.sapiens version 1.3.1 which returns the hg19 grange
# (can be verified in the seqinfo field)
library(Homo.sapiens)
txdb <- Homo.sapiens

hg19_genes <- GenomicFeatures::genes(txdb, columns = "SYMBOL")
