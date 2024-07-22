if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SingleCellExperiment")
BiocManager::install(c("scuttle", "scran", "scater", "uwot", "rtracklayer"))
BiocManager::install("DropletUtils")
BiocManager::install("scater")
BiocManager::install("scRNAseq")
BiocManager::install("batchelor")
BiocManager::install("bluster")
install.packages("uwot")
library(uwot)
library(SingleCellExperiment)
library(scater)
library(scuttle)


# Reading from tabular formats
mat <- as.matrix(read.delim("data/GSE85241_cellsystems_dataset_4donors_updated.csv.gz"))
# scuttle package
library(scuttle)
sparse.mat <- readSparseCounts("data/GSE85241_cellsystems_dataset_4donors_updated.csv.gz")
dim(sparse.mat)
object.size(sparse.mat) < object.size(mat)
# From Excel files
library(R.utils)
excel_file2 = "data/GSE61533_HTSEQ_count_results.xls.gz"
gunzip(excel_file2, destname = "data/GSE61533_HTSEQ_count_results.xls", remove=FALSE, overwrite=TRUE)
library(readxl)
all.counts <- read_excel("data/GSE61533_HTSEQ_count_results.xls")
gene.names <- all.counts$ID
all.counts <- as.matrix(all.counts[,-1])
rownames(all.counts) <- gene.names
head(all.counts)
dim(all.counts)
# From CellRanger output
library(DropletUtils)
untar("data/pbmc4k_raw_gene_bc_matrices.tar.gz", exdir = "data/tenx-pbmc4k")
sce <- read10xCounts("data/tenx-pbmc4k/raw_gene_bc_matrices/GRCh38/")
sce$Barcode

# Construct SingleCellExperiment object
mat <- read.delim("data/E-MTAB-5522/counts_Calero_20160113.tsv", header = T, row.names = 1, check.names = F)
spike.mat <- mat[grepl("^ERCC-", rownames(mat)), ]
gene.length <- mat[, 1]
mat <- as.matrix(mat[,-1])
sce <- SingleCellExperiment(assay = list(counts = mat))
# Another Example
ncells <- 100
u <- matrix(rpois(20000, 5), ncol=ncells)
v <- log2(u + 1)
pca <- matrix(runif(ncells*5), ncells)
tsne <- matrix(rnorm(ncells*2), ncells)
sce <- SingleCellExperiment(assays = list(counts=u, logcounts=v), 
                            reducedDims = SimpleList(PCA=pca, tSNE=tsne))
assayNames(sce)

# Metadata
coldata <- read.delim("data/E-MTAB-5522/E-MTAB-5522.sdrf.txt", check.names = F)
coldata <- coldata[coldata[, "Derived Array Data File"] == "counts_Calero_20160113.tsv", ]
coldata <- DataFrame(genotype=coldata[, "Characteristics[genotype]"], 
                     phenotype=coldata[, "Characteristics[phenotype]"], 
                     spike_in=coldata[, "Factor Value[spike-in addition]"], 
                     row.names = coldata[, "Source Name"])
sce <- SingleCellExperiment(assays = list(counts=mat), colData=coldata)
colData(sce)
head(sce$phenotype)
head(sce$Factor.Value.phenotype.)

sce <- SingleCellExperiment(list(counts=mat))
colData(sce) <- coldata
sce
sce <- SingleCellExperiment(list(counts=mat))
sce$phenotype <- coldata$phenotype
colData(sce)
# rows of the coldata refer to the same cells as the columns of the count matrix
identical(rownames(coldata), colnames(mat))

sce <- scuttle::addPerCellQC(sce)
colData(sce)

rowData(sce)$Length <- gene.length
rowData(sce)
rowRanges(sce)

# subsetting and combining
first.10 <- sce[, 1:10]
ncol(counts(first.10))
colData(first.10)
coding.only <- sce[rowData(sce)$gene_biotype == "protein_coding", ]
counts(coding.only)
rowData(coding.only)

sce2 <- cbind(sce, sce)
ncol(counts(sce2))
sce2 <- rbind(sce, sce)
nrow(counts(sce2))

sce <- scater::logNormCounts(sce)
sce <- scater::runPCA(sce)
dim(reducedDim(sce, "PCA"))
sce <- scater::runTSNE(sce, perplexity = 0.1)
reducedDim(sce, "TSNE")

u <- uwot::umap(t(logcounts(sce)), n_neighbors = 2)
reducedDim(sce, "UMAP_UWOT") <- U
reducedDims(sce)

# workflow example
library(scRNAseq)
# Quality control (using mitochondrial genes)
library(scater)
sce <- MacoskoRetinaData()
is.mito <- grepl("^MT-", rownames(sce))
qcstats <- perCellQCMetrics(sce, subsets = list(Mito=is.mito))
filtered <- quickPerCellQC(qcstats, percent_subsets = "subsets_Mito_percent")
sce <- sce[, !filtered$discard]
# Normalization
sce <- logNormCounts(sce)
# Feature selection, blocking on the individual of origin
library(scran)
dec <- modelGeneVar(sce)
hvg <- getTopHVGs(dec, prop=0.1)
# Batch correction
library(scater)
set.seed(1234)
sce <- runPCA(sce, ncomponents=25, subset_row=hvg)
library(bluster)
colLabels(sce) <- clusterCells(sce, use.dimred = "PCA", BLUSPARAM = NNGraphParam(cluster.fun = "louvain"))
# Visualization
sce <- runUMAP(sce, dimred = "PCA")
plotUMAP(sce, colour_by = "label")

