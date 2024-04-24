##########################################
# Goal: Generate the GRanges object for the different varbins
# at different resoltuions creates each and add to an object
# 
# Output: {genome_name}_grangeslist object used in sysdata.rda
#
##########################################
# Author: Anna K. Casasent
##########################################
library(IRanges)
library(GenomicRanges)
library(AnnotationHub)
library(BiocManager)
library(annotatr)
#library(bioframe)
library(GenomicFeatures)
library(dplyr)
#library(CNEr)
#library(RMariaDB)
#library(TxDb.Mm.mm10.ensembl.reg)

#the_size_list <- c(50000, 100000, 175000, 200000, 250000, 500000, 1000000, 2000000)
#names(the_size_list) <- c("50kb", "100kb", "200kb","175kb", "250kb", "500kb", "1Mb","2Mb")

#test
#the_size_list <- c(50000, 2000000)
#names(the_size_list) <- c("50kb","2Mb")


the_size_list <-c(50000, 100000, 175000, 200000, 250000, 500000, 1000000, 2000000) #c(50000, 100000, 175000, 250000, 500000, 1000000)
names(the_size_list) <- c("50kb", "100kb", "175kb","200kb", "250kb", "500kb", "1Mb","2Mb") #c("50kb", "100kb", "175kb", "250kb", "500kb", "1Mb")
the_aligner <- "bowtie" # using bowtie for now since that is what I tend to run 
the_seqlen <- 48 # using 48 bp length as this is the standard length used by the core 
the_genome <- "mm10"
the_species <- "Mus musculus"
the_sp_abb <- paste0(substring(strsplit(the_species," ")[[1]], 1,1), collapse = "")
the_sp_semiabb <- paste0(substring(the_species,1,1),strsplit(the_species," ")[[1]][2], collapse = "")

# may need to be edits for different genomes
# however this should work for more than just mouse
hub <- AnnotationHub()
my_genome_files <- query(hub, the_species, the_genome)
#EnsemblTxDb <- makeTxDbFromEnsembl(organism = the_species, release = 99)
#the_database_name <- paste0("TxDb.",the_sp_abb,".",the_genome,".ensembl.reg")
#the_database_name <- paste0("org",the_sp_abb,".",the_genome,".eg.db")
the_database_name <- paste0("TxDb.",the_sp_semiabb,".UCSC.",the_genome,".knownGene")

if (!require(the_database_name, character.only = TRUE)) {
  BiocManager::install(the_database_name, force = TRUE)
  library(the_database_name, character.only = TRUE)
}

the_BS_name <- paste0("BSgenome.",the_sp_semiabb,".UCSC.",the_genome)
if (!require(the_BS_name, character.only = TRUE)) {
  BiocManager::install(the_BS_name, force = TRUE)
  library(the_BS_name, character.only = TRUE)
}
  
txdb <- get(the_database_name) #TxDb.Mmusculus.UCSC.mm10.knownGene
gene_ranges <- genes(txdb)
promoter_ranges <- promoters(gene_ranges, upstream = 2000, downstream = 200)
my_BSgenome <- get(the_BS_name)


### Function to get chrom arm
get_chrom_arm <- function(gr, the_hub, the_species, the_genome) {
my_hub <- subset(
  the_hub,
  (the_hub$species == the_species) & (the_hub$genome == the_genome)
)

my_hub_df <- data.frame(
  ahid = my_hub$ah_id,
  title = my_hub$title
)
my_adid_band <- my_hub_df$ahid[grep("Band",my_hub_df$title)]
print(my_adid_band)
if(length(my_adid_band) == 1)
{
  g_cytoband <- my_hub[[my_adid_band]]
  print(g_cytoband)
}
if(length(my_adid_band) != 1)
{
  print("unable to add cytobands")
  return(NULL)
  break
}

olaps <- findOverlaps(gr, g_cytoband)

df_olaps <- as.data.frame(olaps)

df_olaps_dist <- df_olaps %>% distinct(queryHits, .keep_all = T)

g_cyto_df <- as.data.frame(g_cytoband[df_olaps_dist$subjectHits])

g_cyto_df <- g_cyto_df %>%
  dplyr::mutate(arm = stringr::str_extract(name, "[pq]"))

return(g_cyto_df$arm)
}

### Function for cpg islands 

get_chrom_cpi <- function(the_species, the_genome) {
  library(annotatr)
  library(GenomicRanges)
  
  my_annotations = build_annotations(genome = the_genome, annotations = paste0(the_genome,"_cpgs"))
  
  # Access the CpG island annotations
  cpg_islands = my_annotations[my_annotations$type == paste0(the_genome,"_cpg_islands")] 
  return(cpg_islands)
}

### Function to get chrom arm
get_chrom_arm <- function(gr, the_hub, the_species, the_genome) {
my_hub <- subset(
  the_hub,
  (the_hub$species == the_species) & (the_hub$genome == the_genome)
)

my_hub_df <- data.frame(
  ahid = my_hub$ah_id,
  title = my_hub$title
)
my_adid_band <- my_hub_df$ahid[grep("Band",my_hub_df$title)]
print(my_adid_band)
if(length(my_adid_band) == 1)
{
  g_cytoband <- my_hub[[my_adid_band]]
  print(g_cytoband)
}
if(length(my_adid_band) != 1)
{
  print("unable to add cytobands")
  return(NULL)
  break
}

olaps <- findOverlaps(gr, g_cytoband)

df_olaps <- as.data.frame(olaps)

df_olaps_dist <- df_olaps %>% distinct(queryHits, .keep_all = T)

g_cyto_df <- as.data.frame(g_cytoband[df_olaps_dist$subjectHits])

g_cyto_df <- g_cyto_df %>%
  dplyr::mutate(arm = stringr::str_extract(name, "[pq]"))

return(g_cyto_df$arm)
}
#cgi_ranges <- findCpGIslands(promoter_ranges, my_BSgenome)

#### Functions 

GCContentByRange <- function(ranges, genome) {
  ## Check inputs
  #if (!is(ranges, "GRanges")) {
  #  stop("'ranges' must be a GRanges object")
  #}
  #if (!is(genome, "BSgenome")) {
  #  stop("'genome' must be a BSgenome object")
  #}
  
  ## Extract the sequences
  seqs <- getSeq(genome, ranges)
  
  ## Calculate GC content
  # not a percentage
  gc_content <-  (vcountPattern("G", seqs) + vcountPattern("C", seqs)) / 
                                (width(seqs))
  return(gc_content)
}


#### Start process
# expect a folder in data-raw 
# folder should contain 
# lengths - file that lists all chromosomes and their lengths
# centromeres - file that lists the centromeres start/stop
# all the scafold based on a ginko output for each chromosome, using length for sequence read, aligner and variable

the_centromeres <- read.delim(paste0("data-raw/",the_genome,"/centromeres"), header = FALSE)
colnames(the_centromeres) <- c("CHR","START","END")
the_lengths <- read.delim(paste0("data-raw/",the_genome,"/lengths"), sep=" ", header = FALSE)
colnames(the_lengths) <- c("CHR","END")
cpi_ranges <- get_chrom_cpi(the_species, the_genome)

for(my_size in the_size_list)
{
  my_size_char <- format(the_size_list[which(the_size_list==my_size)], scientific = FALSE)
  
  # contained variable bins to within this setting
  my_min_size <- floor(as.numeric(my_size_char)/2) # removing bins half the size of other bins
  my_max_size <- ceiling(as.numeric(my_size_char)*3)# removing bins three time the size of other bins
  
  #my_breakpoints <- read.delim(paste0("data-raw/",the_genome,"/variable_",my_size, "_", the_seqlen,"_", the_aligner))
  my_breakpoints <- read.delim(paste0("data-raw/",the_genome,"/variable_",format(my_size, scientific = FALSE), "_", the_seqlen,"_", the_aligner))
  colnames(my_breakpoints) <- c("CHR", "END")
  my_varbin_file <- data.frame("bin.chrom"= my_breakpoints$CHR, "bin.start"=0, "bin.end"= my_breakpoints$END, "bin.length"=0)
  
  # input 
  # # chrom END
  for(index in 1:nrow(my_varbin_file))
    {
    my_chr <- my_varbin_file$bin.chr[index]
    my_centromeres_df <- the_centromeres[which(the_centromeres$CHR == my_chr),]
    my_centromeres <- IRanges(start = my_centromeres_df$START, end = my_centromeres_df$END)
    my_length <- the_lengths[which(the_lengths$CHR == my_chr),]
    
    # not in the first or the last index 
    if(index > 1)
    {
      before_same = my_varbin_file$bin.chrom[index] == my_varbin_file$bin.chrom[index-1]
      # after_same = my_varbin_file$bin.chrom[index] == my_varbin_file$bin.chrom[index+1]
      # check chrom
      if(before_same)
      {
        # if chrom before index is the same both before and after
        # Note we are not checking for centemeres here... We might need to
        #print('changed the start bin')
        my_varbin_file$bin.start[index] = my_varbin_file$bin.end[index-1]
        my_varbin_file$bin.length[index] = my_varbin_file$bin.end[index] - my_varbin_file$bin.start[index] 
      }
    }
    # check that bin.start and bin.end does overlap the centermere 
    my_bin <- IRanges(start = my_varbin_file$bin.start[index], end = my_varbin_file$bin.end[index])
    if(overlapsAny(my_centromeres, my_bin))
    {
      #print(paste("My bin chr: ", my_varbin_file$bin.chrom[index]))
      #print(paste("My bin start: ", my_varbin_file$bin.start[index]))
      #print(paste("My bin end: ", my_varbin_file$bin.end[index]))
      
      #change my bin start to end of centermere 
      if(my_varbin_file$bin.start[index] <= my_centromeres_df$END & my_varbin_file$bin.end[index] > my_centromeres_df$END)
      {
        #print("Overlap: change start of bin")
        my_varbin_file$bin.start[index] <- my_centromeres_df$END
        my_varbin_file$bin.length[index] = my_varbin_file$bin.end[index] - my_varbin_file$bin.start[index] 
      }
      if(my_varbin_file$bin.start[index] < my_centromeres_df$START & my_varbin_file$bin.end[index] <= my_centromeres_df$END)
      {
        #print("Overlap: change end of bin")
        my_varbin_file$bin.end[index] <- my_centromeres_df$START
        my_varbin_file$bin.length[index] = my_varbin_file$bin.end[index] - my_varbin_file$bin.start[index] 
      }
    }
    
  }
  # remove bins that are less than half the size of the normal bin
  rm_bin_index <- which(my_varbin_file$bin.length <my_min_size | my_varbin_file$bin.length >my_max_size)
  if(length(rm_bin_index) >= 1)
    {
    my_varbin_file <- my_varbin_file[-rm_bin_index,]
    }
  
  # create a genome ranges scaffold 
  my_gr <- makeGRangesFromDataFrame(
    df = my_varbin_file,
    keep.extra.columns = T,
    seqnames.field = "bin.chrom",
    start.field = "bin.start",
    end.field = "bin.end",
    starts.in.df.are.0based = TRUE
  )
  
  # Add to each res
  # gene_count
  # for each in 
  #abspos
  mcols(my_gr)$gene_count <- countOverlaps(my_gr, gene_ranges, ignore.strand=FALSE)
  mcols(my_gr)$cgi_count <- countOverlaps(my_gr, cpi_ranges, ignore.strand=FALSE)
  
  # did not get the distanceToTelomere
  # mcols(my_gr)$dist_telomere <- distanceToTelomere(my_gr, my_BSgenome)
  
  mcols(my_gr)$gc_content <- GCContentByRange(my_gr, my_BSgenome)
  # make an ablouse position by stacking them end to end 
  mcols(my_gr)$abspos <- end(my_gr) + sapply(as.vector(chrom(my_gr)),function(my_chr)
      {
      index_to <- which(the_lengths$CHR== my_chr)
      to_add <- 0
      if(index_to - 1 > 0)
      {
        to_add <- sum(the_lengths$END[1:(index_to - 1)])
      }
      to_add
      })
    #countOverlaps(my_gr, gene_ranges, ignore.strand=FALSE)
  
  #internal function
  mcols(my_gr)$arm <- get_chrom_arm(my_gr,  the_hub = hub, the_species = the_species, the_genome = the_genome)
  
  print(paste0("created object called: res_",names(my_size_char)))
  # issue to fix not getting the name for my size
  assign(x=paste0("res_",names(my_size_char)),value=my_gr)

  }




# make the full object 
list_gr <- lapply(names(the_size_list), function(my_res)
{
  get(paste0("res_",my_res))
})

names(list_gr) <- paste(the_genome,names(the_size_list), sep="_")
varbin_genome_grangeslist <- GRangesList(unlist(list_gr))

cat(length(varbin_genome_grangeslist), "resolutions are available:\n")
cat(names(varbin_genome_grangeslist), "\n")
cat("How many varbins are available per resolution?:\n")
print(sapply(varbin_genome_grangeslist, length))

genome_seqinfo <- GenomeInfoDb::Seqinfo(genome = the_genome)
GenomeInfoDb::seqinfo(varbin_genome_grangeslist) <- genome_seqinfo



assign(paste0(the_genome, "_grangeslist"), value=varbin_genome_grangeslist )

my_grangelist <- get(paste0(the_genome, "_grangeslist"))
save(my_grangelist, file=paste0("data/",the_genome,"_grangeslist.rda"))

# in order to create different genome varbins of about these sizes we need to 
# get the GC content and from the genome of interest