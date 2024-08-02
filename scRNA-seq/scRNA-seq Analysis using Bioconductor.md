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
## Quality control
1. Low-quality can be resulted from cell damage during dissociation or failure in library preparation. These can be identified with low total counts, few expressed genes and high mitochondrial or spike-in proportions. 
    1. The damaged cells form their own distinct clusters because of increased mitochondrial proportions or enrichment for nuclear RNAs after cell damage. 
    2. They distort the characterization of population heterogeneity during variance estimation or principal component analysis. The first few principal component will capture differences in quality rather than biology. Similarly, genes with the largest variances will be driven by differences between low- and high-quality cells. 
    3. They contain genes that appear to be strongly "upregulated" due to aggressive scaling to normalize for small library sizes. Increased scaling in low-quality libraries transforms small counts for these transcripts in large normalized expression values, resulting in upregulation compared to other cells.
2. QC metrics
Several common QC metrics are used to identify low-quality cells based on their expression profiles.
    1. The libary size. Cells with small library sizes are of low quality as the RNA has been lost at some point during preparation. 
    2. The number of expressed genes. Cells with very few expressed genes are likely to be of poor quality as the diverse transcript population has not been captured. 
    3. The proportion of reads mapped to spike-in transcripts is calculated relative to the total count across all features for each cell. As the same amount of spike-in RNA should have been added to each cell, any enrichment in spike-in counts is symptomatic of loss of endogenous RNA. Thus, high proportions are indicative of poor-quality cells. 
    4. In the absence of spike-in transcripts, the proportion of reads mapped to genes in the mitochondrial genome can be used. High proportions are indicative of poor-quality cells.
3. Using QC metrics to identify low-quality cells
    1. Using fixed threshold. We may consider cells to be low quality if the library size is below 100000 reads; expresses fewer than 5000 genes; have spike-in proportion above 10% or mitochondrial proportions above 10%. Thresholds for read count-based data are simply not applicable for UMI-based data, and vice versa. Differences in mitochondrial activity or total RNA content require constant adjustment of mitochondrial and spike-in thresholds.
    2. Using adaptive threshold. 
        1. Identifying outliers
        To obtain an adaptive threshold, we assume that most of the dataset consists of high-quality cells. We then identify cell that are outlier for the various QC metrics based on the median absolute deviation (MAD) from the median value of each metric across all cells. Specifically, a value is considered an outlier if it is more than 3 MADs from the median in the problematic direction. 
        2. Assumptions of outlier detection
        Outlier detection assumes that most cells are of acceptable quality. Another assumption is that the QC metrics are independent of the biological state of each cell.
        3. Considering experimental factors
4. Checking diagnostic plots
5. Removing low-quality cells
6. Marking low-quality cells 

























