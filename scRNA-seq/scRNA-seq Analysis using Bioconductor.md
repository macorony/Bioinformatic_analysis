# scRNA-seq Analysis using Bioconductor
[Credit resource] (https://www.singlecellcourse.org/scrna-seq-analysis-with-bioconductor.html)
## About bioconductor
1. Bioconductor releases its packages as a cohort on a half-yearly cycle ensuring that different packages work smoothly.  

## Packages
1. SingleCellExperiment (R package)
2. Seurat (R package)
3. Scanpy (python package)
## Formats of read counts 
1. tabular formats
```R
mat <- as.matrix(read.delim("dir/file"))

```
