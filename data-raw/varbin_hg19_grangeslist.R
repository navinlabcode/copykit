##########################################
# Goal: Generate the GRangesList object for the hg19 varbins at 100k and
# 200k resolution.
#
# Inputs:
# https://github.com/navinlabcode/CNV_pipeline/tree/master/lib/
# - varbin.gc.content.100k.bowtie.k50.hg19.bin_not_removed.txt
# - varbin.gc.content.200k.bowtie.k50.hg19.bin_not_removed.txt
#
# Output: 'sysdata.rda' which additionally contains a GRangesList object
# named 'varbin_hg19_grangeslist'.
#
# 1. Convert the raw coordinate from a mixture of 0-based and 1-based to 1-based.
# For example, start=0, end=977835, bin_length=977386.
# Equivalently, pure 0-based: [0, 977836); pure 1-based: [1, 977836]
# 2. Provide the correct genomic locations and attache meta information (e.g.,
# GC content) for all varbins.
# 3. Chr-style is chr1, chr2, ..., chrX, chrY. (No mitochondrion chromosome)
#
##########################################
# Author: Yun Yan (yun.yan@uth.tmc.edu)
##########################################
library(devtools)
library(usethis)
library(GenomicRanges)
library(GenomeInfoDb)

# VarBins annotations are from the official pipeline's lib: UCSC style (chr1, chrX, chrY)
fpath_varbin <- list(
  `res_100k`='/volumes/seq/code/PIPELINES/CNA_pipeline_v1.4/lib/varbin.gc.content.100k.bowtie.k50.hg19.bin_not_removed.txt',
  `res_200k`='/volumes/seq/code/PIPELINES/CNA_pipeline_v1.4/lib/varbin.gc.content.200k.bowtie.k50.hg19.bin_not_removed.txt'
)

make_granges_from_varbin_file <- function(x){

  df <- read.delim(x, header = T, stringsAsFactors = F)
  colnames(df) <- gsub(pattern = '\\.', replacement = '_', x=colnames(df))
  # Input:
  # bin.chrom bin.start bin.end bin.length
  #       chr1         0  977835     977836
  #       chr1    977836 1200862     223027
  #       chr1   1200863 1455237     254375
  #       chr1   1455238 1758056     302819
  # ...
  #       chr2         0  237130     237131
  #       chr2    237131  454244     217114
  #       chr2    454245  660285     206041

  ## Convert to pure 0-based coordination
  df$bin_end <- df$bin_end + 1

  ## Convert to a legal GRange object (automatically switched to 1-based)
  gr <- makeGRangesFromDataFrame(
    df=df,
    keep.extra.columns = T,
    seqnames.field = 'bin_chrom',
    start.field = 'bin_start',
    end.field = 'bin_end',
    starts.in.df.are.0based = TRUE)
  return(gr)

  # Output:
  # seqnames          ranges strand | bin_length
  # [1]     chr1        1-977836      * |     977836
  # [2]     chr1  977837-1200863      * |     223027
  # [3]     chr1 1200864-1455238      * |     254375
  # ...
  # [1]     chr2      1-237131      * |     237131
  # [2]     chr2 237132-454245      * |     217114
  # [3]     chr2 454246-660286      * |     206041
}

list_gr <- lapply(fpath_varbin, make_granges_from_varbin_file)
varbin_hg19_grangeslist <- GRangesList(unlist(list_gr))
cat(length(varbin_hg19_grangeslist), 'resolutions are available:\n')
cat(names(varbin_hg19_grangeslist), '\n')
cat('How many varbins are available per resolution?:\n')
print(sapply(varbin_hg19_grangeslist, length))

hg19_seqinfo <- GenomeInfoDb::Seqinfo(genome = 'hg19')
GenomeInfoDb::seqinfo(varbin_hg19_grangeslist) <- hg19_seqinfo

# Export to package internal data
# devtools::load_all('./')
# usethis::use_data(
#   major_palette, minor_palette, hg19_genes, ## previously available
#   varbin_hg19_grangeslist, ## Added in this script
#   internal = TRUE, overwrite = TRUE)
