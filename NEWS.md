# copykit (development version)

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
