##########################################
# Goal: Generate the GRanges object for the different varbins
# containing gene symbols
#
# Output: {genome_name}_genes object used in sysdata.rda
#
##########################################
# Author: Darlan Conterno Minussi & Anna K. Casasent
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
# mm10 genes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

library(AnnotationHub)
# BiocManager::install("Mus.musculus")
library(Mus.musculus)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
# BiocManager::install("Organism.dplyr")
library(Organism.dplyr)

src <- src_organism("TxDb.Mmusculus.UCSC.mm10.knownGene")
mm10_genes <- genes(src, columns = "symbol")
saveRDS(mm10_genes, "data/mm10_genes.rda")
usethis::use_data(mm10_genes)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
## mm39 genes
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
## This requires R 4.3 
## Plans to add when updating R
#library(AnnotationHub)
## BiocManager::install("Mus.musculus")
#library(Mus.musculus)
## BiocManager::install("TxDb.Mmusculus.UCSC.mm39.knownGene")
#library(TxDb.Mmusculus.UCSC.mm39.knownGene)
#library(Organism.dplyr)
#
#src <- src_organism("TxDb.Mmusculus.UCSC.mm39.knownGene")
#mm39_genes <- genes(src, columns = "symbol")
#saveRDS(mm39_genes, "data/mm39_genes.rda")
#usethis::use_data(mm39_genes)

