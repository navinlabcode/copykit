# copykit 0.1.2

**New features**
* New cell smoothing method with the function `knnSmooth()`. Uses k-nearest neighbors to smooth cells profiles and re-segments the datasets obtaining cleaner copy number profiles, with reduced overdispersion and improving downstream analysis. (Thanks to [Runmin Wei]([)https://github.com/WandeRum) for the helpful discussion.)
* scquantum method is available for `calcInteger()` and is now a CopyKit import ([scquantum](https://github.com/navinlabcode/scquantum) is a single cell ploidy estimation tool developed by [Alexander Davis](https://github.com/alex-l-m))
* `calcInteger()` now accepts option methdo = 'metadata'. To use this option the user can add custom values of ploidy to every cell in the colData column 'ploidy' and run `calcInteger(ck, method = 'metadata')` to obtain the integer matrix on the CopyKit object.
* runVst allows selection of the assay for the transformation

**Changes**
* `plotHeatmap()` order_cells argument now defaults to NULL. NULL option respects the order of the CopyKit object. order_cells argument can be set to 'consensus_tree' and 'hclust'.

* Method 'scquantum' from `calcInteger()` adds 3 elements to the colData.
1. ploidy: contanining the inferred ploidy call for each cell
2. confidence_ratio: ratio from scquantum inferred ploidy to scquantum theoretical ploidy
3. ploidy_score: Score derived from the confidence ratio. Values closer to 0 indicate a better fit of the ploidy call

* Significance thresholds for CBS alpha segmentation and Merge levels were reduced to increase sensitivity to focal amplifications.

**Removed**
* option 'phylogeny' from function argument `plotHeatmap()' 'order_cells' has been removed.

**Bug Fixes**
* Fixed error in plotGeneCopy not returning plots with geom violin and barplot. (Thanks to @Romeo1-1)
* Fixed error in plotGeneCopy with duplicated sample names on a merged object. Now it warns the user of merged sample names. (Thanks to @Romeo1-1)
* Allowing control of parameter merge_levels_alpha on `runVarbin()` and `runSegmentation()` to control the significance level of merge levels when merging not significant segments.


# copykit 0.1.1

* Reduced quality of heatmap raster that could quickly use all magick cache
* Fixed hg38 scaffold issue for lower resolutions 500kb, 1Mb and 2.8Mb in which the quality control of low quality bins was too strict and causing problems especially on chromosome X. (Thanks to @Romeo1-1)

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
