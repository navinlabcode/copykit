##########################################
# Goal: Generate the GRanges object for the hg19 varbins
# containing gene symbols
#
# # Output: hg19_genes object used in sysdata.rda
#
##########################################
# Author: Darlan Conterno Minussi
##########################################

library(AnnotationHub)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(Organism.dplyr)

src <- src_organism("TxDb.Hsapiens.UCSC.hg19.knownGene")
hg19_genes <- genes(src, columns = "symbol")

usethis::use_data(hg19_genes)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Tue Feb  2 10:04:12 2021
# hg38 genes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Tue Feb  2 10:04:21 2021

library(AnnotationHub)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(Organism.dplyr)

src <- src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene")
hg38_genes <- genes(src, columns = "symbol")

usethis::use_data(hg38_genes)
