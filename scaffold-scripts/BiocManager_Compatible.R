if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.19")
BiocManager::install(c("AnnotationHub", "GenomicRanges", "GenomicFeatures", "annotatr"), 
                     version = "3.19")


#pkgs <- rownames(installed.packages())
#BiocManager::install(pkgs, type = "source", checkBuilt = TRUE, force=TRUE, version="3.19") #"3.15"

#pkg <- c("bindSC", "CellTrek", "copykat", "copykit", "DoubletFinder", 
#         "maptools", "MASS", "Matrix", "presto", "rgeos", 
#         "sceasy", "SeuratData", "spatstat.core", "Cairo", "DNAcopy", "Rsubread", 
#         "qqconf", "textshaping", "igraph", "gert", "influenceR", "leiden", 
#         "RMariaDB", "bluster", "DiagrammeR", "nloptr", "randomForestSRC", 
#         "Seurat", "ggtree")

install.packages("remotes")