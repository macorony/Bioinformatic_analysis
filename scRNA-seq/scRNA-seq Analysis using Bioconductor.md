# scRNA-seq Analysis using Bioconductor
[Credit resource] (https://www.singlecellcourse.org/scrna-seq-analysis-with-bioconductor.html)
## About bioconductor
1. Bioconductor releases its packages as a cohort on a half-yearly cycle ensuring that different packages work smoothly.  

## Packages
1. SingleCellExperiment (R package)
2. Seurat (R package)
3. Scanpy (python package)

## Formats of read counts 
1. tabular format
```R
mat <- as.matrix(read.delim("dir/file"))
```
or
```R
library(scuttle)
sparse.mat <- readSparseCounts("dir/file")
```
2. Cellranger format
```R
library(DropletUtils)
sce <- read10xCounts("dir/file")
```
3. HDF5-based formats
4. Binary-compressed format (Salmon-Alevin output)
Salmon-alevin produces a per-cell level gene-count matrixt in a binary-compressed format with the row and column indexes in a separate file. 
```R
library(Seurat)
library(tximport)
txi <- tximport("dir/file", type="alevin")
pbmc <- CreateSeuratObject(counts= txi$counts, min.cells=3, min.features=200, project="10X_PBMC")
```

## The SingleCellExperiment class



