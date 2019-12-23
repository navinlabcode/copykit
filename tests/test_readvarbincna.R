library(copykit) ## load official version
library(devtools)
library(GenomicRanges)
library(SingleCellExperiment)
dir_demo <- "/volumes/seq/projects/CNA_projects/DT_CNA/snap_frozen/Breast/TNBC/TN20/TN20_output/final_result/"

ck <- readVarbinCNA(dir_demo, remove_Y = T)

gr <- rowRanges(ck)
table(seqnames(gr))
gr_chr1 <- gr[seqnames(gr) == '1', ]

bin_width <- width(gr_chr1)
cat('Max bin size is', max(bin_width), 'which is really strange as the VarBin defines a 200k width.')
# See the strange varbin:
print(gr_chr1[width(gr_chr1) == max(bin_width), ])
# GRanges object with 1 range and 1 metadata column:
#   seqnames              ranges strand |    abspos
# <Rle>           <IRanges>  <Rle> | <numeric>
#   515        1 120291655-145520501      * | 120291655

# The above varbin is a concatenation of 5 varbins.
# chrom start end length
# chr1	120055255	120291654	236400
# chr1	120291655	120508459	216805
# chr1	120508460	121326396	817937
# chr1	121326397	144071760	22745364
# chr1	144071761	145073280	1001520
# chr1	145073281	145520501	447221
# chr1	145520502	145823581	303080


