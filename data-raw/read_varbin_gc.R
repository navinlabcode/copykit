##########################################
# Goal: Generate the GRanges object for the varbins (200k) of hg19.
#
# Input: inst/extdata/varbin.gc.content.200k.bowtie.k50.hg19.bin_not_removed.txt
# Output: inst/extdata/hg19_200k_varbins_full.granges.rds
#
## Why this script and the internal data?
## 1. The coordination here is a mixture of 0-based and 1-based.
## For example, start=0, end=977835, bin_length=977386.
## Equivalently, pure 0-based: [0, 977836); pure 1-based: [1, 977836]
## 2. This provides the ground-truth of the genomic locaitons of all varbins.
#
##########################################
# Author: Yun Yan (yun.yan@uth.tmc.edu)
##########################################
library(devtools)
library(usethis)
library(GenomicRanges)

fpath_varbin_gc <- system.file(
  "extdata", "varbin.gc.content.200k.bowtie.k50.hg19.bin_not_removed.txt", package = "copykit")

df <- read.delim(fpath_varbin_gc, header = T, stringsAsFactors = F)
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
gr_hg19_200k_varbin_full <- makeGRangesFromDataFrame(
  df=df,
  keep.extra.columns = T,
  seqnames.field = 'bin_chrom',
  start.field = 'bin_start',
  end.field = 'bin_end',
  starts.in.df.are.0based = TRUE)

saveRDS(object = gr_hg19_200k_varbin_full,
        file = './inst/extdata/hg19_200k_varbins_full.granges.rds')

# Output:
# seqnames          ranges strand | bin_length
# [1]     chr1        1-977836      * |     977836
# [2]     chr1  977837-1200863      * |     223027
# [3]     chr1 1200864-1455238      * |     254375
# ...
# [1]     chr2      1-237131      * |     237131
# [2]     chr2 237132-454245      * |     217114
# [3]     chr2 454246-660286      * |     206041
