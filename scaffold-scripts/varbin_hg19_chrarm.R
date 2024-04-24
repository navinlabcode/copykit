##########################################
# Goal: Add chr arm information to hg19 scaffold.
#
# Inputs:
# hg19_rg.rda
#
# Output: use_data('hg19_rg.rda')
#

library(AnnotationHub)
library(GenomicRanges)
library(dplyr)
library(stringr)
library(devtools)

load_all()

hub <- AnnotationHub()

hub_hg19 <- subset(
    hub,
    (hub$species == "Homo sapiens") & (hub$genome == "hg19")
)

hub_df <- data.frame(
    ahid = hub_hg19$ah_id,
    title = hub_hg19$title
)

g_cytoband <- hub_hg19[["AH5012"]]

g_hg19 <- makeGRangesFromDataFrame(hg19_rg, keep.extra.columns = T)

olaps <- findOverlaps(g_hg19, g_cytoband)

df_olaps <- as.data.frame(olaps)

df_olaps_dist <- df_olaps %>% distinct(queryHits, .keep_all = T)

g_cyto_df <- as.data.frame(g_cytoband[df_olaps_dist$subjectHits])

g_cyto_df <- g_cyto_df %>%
    dplyr::mutate(arm = stringr::str_extract(name, "[pq]"))

g_hg19$arm <- g_cyto_df$arm

hg19_rg <- g_hg19
use_data(hg19_rg, overwrite = TRUE)
