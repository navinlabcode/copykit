# copykit 0.1.1

* Reduced quality of heatmap raster that could quickly use all magick cache
* Fixed hg38 scaffold issue for lower resolutions 500kb, 1Mb and 2.8Mb in which the quality control of low quality bins was too strict and causing problems especially on chromosome X. 

# copykit 0.1.0

* CopyKit goes 0.1.0.

# copykit 0.0.0.9036

* Added a `NEWS.md` file to track changes to the package.

# copykit 0.0.0.9037

* Changed values of resolution argument in `runVarbin()` to more accurately 
reflect the variable genomic scaffolds. 
Resolutions are: '55kb', '110kb', '195kb', '220kb', '280kb', '500kb', '1Mb', '2.8Kb'

# copykit 0.0.0.9038

* Adding PCA and PCA related functions
 - `runPca()`
 - `plotPca()`
 - `plotScree()`
 
* Clustering functions can use n dimensions with either UMAP or PCA
 - added argument `ncomponents` to `findSuggestedK()` and `findClusters` 

# copykit 0.0.0.9039

* Fixed a critical bug in runSegmentation to set the correct log base call
during merge levels.

# copykit 0.0.0.9040
* Adding raster arguments to plotHeatmap().
