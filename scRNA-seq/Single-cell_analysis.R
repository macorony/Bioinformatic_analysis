if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SingleCellExperiment")
BiocManager::install(c("scuttle", "scran", "scater", "uwot", "rtracklayer"))
BiocManager::install("DropletUtils")
BiocManager::install("scater")
BiocManager::install("scRNAseq")
BiocManager::install("batchelor")
BiocManager::install("bluster")
BiocManager::install("ensembldb")
install.packages("uwot")
library(uwot)
library(SingleCellExperiment)
library(scater)
library(scuttle)
library(ensembldb)


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

# human pancreas single cell analysis
# Data loading
library(scRNAseq)
sce.seger <- SegerstolpePancreasData()
library(AnnotationHub)
edb <- AnnotationHub()[["AH73881"]]
symbols <- rowData(sce.seger)$symbol
ens.id <- mapIds(edb, keys=symbols, keytype="SYMBOL", column="GENEID")
ens.id <- ifelse(is.na(ens.id), symbols, ens.id)

keep <- !duplicated(ens.id)
sce.seger <- sce.seger[keep,]
rownames(sce.seger) <- ens.id[keep]

emtab.meta <- colData(sce.seger)[, c("cell type", "disease", "individual", "single cell well quality")]
colnames(emtab.meta) <- c("CellType", "Disease", "Donor", "Quality")
colData(sce.seger) <- emtab.meta

sce.seger$CellType <- gsub(" cell", "", sce.seger$CellType)
sce.seger$CellType <- paste0(
  toupper(substr(sce.seger$CellType, 1, 1)),
  substring(sce.seger$CellType, 2))

# Quality control
unfiltered <- sce.seger
low.qual <- sce.seger$Quality == "low quality cell"

library(scater)
stats <- perCellQCMetrics(sce.seger)
qc <- quickPerCellQC(stats, percent_subsets="altexps_ERCC_percent", batch=sce.seger$Donor, 
                     subset=!sce.seger$Donor %in% c("HP1504901", "HP1509101"))
sce.seger <- sce.seger[, !(qc$discard | low.qual)]

colData(unfiltered) <- cbind(colData(unfiltered), stats)
unfiltered$discard <- qc$discard

gridExtra::grid.arrange(
  plotColData(unfiltered, x="Donor", y="sum", colour_by = "discard") +
    scale_y_log10() + ggtitle("Total count") +
    theme(axis.text.x = element_text(angle = 90)),
  plotColData(unfiltered, x="Donor", y="detected", colour_by = "discard") +
    scale_y_log10() + ggtitle("Detected features") +
    theme(axis.text.x = element_text(angle = 90)),
  plotColData(unfiltered, x="Donor", y="altexps_ERCC_percent", colour_by = "discard") +
    ggtitle("ERCC percent") +
    theme(axis.text.x = element_text(angle = 90)), 
  ncol = 2
)

# Normalization
library(scran)
clusters <- quickCluster(sce.seger)
sce.seger <- computeSumFactors(sce.seger, clusters=clusters)
sce.seger <- logNormCounts(sce.seger)

summary(sizeFactors(sce.seger))
plot(librarySizeFactors(sce.seger), sizeFactors(sce.seger), pch=16, 
     xlab = "Library size factors", ylab = "Deconvolution factors", log="xy")

# Variance modeling
for.hvg <- sce.seger[, librarySizeFactors(altExp(sce.seger)) > 0 & sce.seger$Donor!="AZ"]
dec.seger <- modelGeneVarWithSpikes(for.hvg, "ERCC", block=for.hvg$Donor)
chosen.hvgs <- getTopHVGs(dec.seger, n=2000)

par(mfrow=c(3,3))
blocked.stats <- dec.seger$per.block
for (i in colnames(blocked.stats)) {
  current <- blocked.stats[[i]]
  plot(current$mean, current$total, main=i, pch=16, cex=0.5, 
       xlab="Mean of log-expression", ylab="Variance of log-expression")
  curfit <- metadata(current)
  points(curfit$mean, curfit$var, col="red", pch=16)
  curve(curfit$trend(x), col="dodgerblue", add=T, lwd=2)
}

# Dimensionality reduction
library(BiocSingular)
set.seed(100)
sce.seger <- runPCA(sce.seger, subset_row=chosen.hvgs, ncomponents=25)
sce.seger <- runTSNE(sce.seger, dimred="PCA")

# Clustering
library(bluster)
clust.out <- clusterRows(reducedDim(sce.seger, "PCA"), NNGraphParam(), full = T)
snn.gr <- clust.out$objects$graph
colLabels(sce.seger) <- clust.out$clusters
tab <- table(Cluster=colLabels(sce.seger), Donor=sce.seger$Donor)
library(pheatmap)
pheatmap(log10(tab+10), color=viridis::viridis(100))

gridExtra::grid.arrange(
  plotTSNE(sce.seger, colour_by = "label"), 
  plotTSNE(sce.seger, colour_by = "Donor"), 
  ncol=2
)

library(batchelor)
set.seed(1000)
corrected <- fastMNN(sce.seger, batch=sce.seger$Donor, subset.row=chosen.hvgs)
corrected <- runTSNE(corrected, dimred="corrected")
colLabels(corrected) <- clusterRows(reducedDim(corrected, "corrected"), NNGraphParam())
tab <- table(Cluster=colLabels(corrected), Donor=corrected$batch)
tab

gridExtra::grid.arrange(
  plotTSNE(corrected, colour_by = "label"), 
  plotTSNE(corrected, colour_by = "batch"),
  ncol = 2
)

# Lung single cell analysis (Smart-seq2)
# Data loading
library(scRNAseq)
sce.416b <- LunSpikeInData(which="416b")
sce.416b$block <- factor(sce.416b$block)
library(AnnotationHub)
ah <- AnnotationHub()
query(ah, c("Mus musculus","Ensembl 98"))

# ens.mm.v97 <- AnnotationHub()[["AH73905"]]
ens.mm.v98 <- AnnotationHub()[["AH75036"]]

rowData(sce.416b)$ENSEMBL <- rownames(sce.416b) 
rowData(sce.416b)$SYMBOL <- mapIds(ens.mm.v98, keys=rownames(sce.416b), 
                                   keytype = "GENEID", column = "SYMBOL")
rowData(sce.416b)$SEQNAME <- mapIds(ens.mm.v98, keys = rownames(sce.416b), 
                                    keytype = "GENEID", column = "SEQNAME")
library(scater)
length(rowData(sce.416b)$ENSEMBL)
length(rowData(sce.416b)$SYMBOL)
rownames(sce.416b) <- uniquifyFeatureNames(rowData(sce.416b)$ENSEMBL, rowData(sce.416b)$SYMBOL)
# Quality control
unfiltered <- sce.416b

mito <- which(rowData(sce.416b)$SEQNAME == "MT")
stats <- perCellQCMetrics(sce.416b, subsets = list(Mt=mito))
qc <- quickPerCellQC(stats, percent_subsets=c("subsets_Mt_percent", 
                                              "altexps_ERCC_percent"), batch=sce.416b$block)
sce.416b <- sce.416b[, !qc$discard]

colData(unfiltered) <- cbind(colData(unfiltered), stats)
unfiltered$block <- factor(unfiltered$block)
unfiltered$discard <- qc$discard

gridExtra::grid.arrange(
  plotColData(unfiltered, x="block", y="sum", 
              colour_by="discard") + scale_y_log10() + ggtitle("Total count"), 
  plotColData(unfiltered, x="block", y="detected", 
              colour_by="discard") + scale_y_log10() + ggtitle("Detected features"), 
  plotColData(unfiltered, x="block", y="subsets_Mt_percent", 
              colour_by="discard") + ggtitle("Mito percent"),
  plotColData(unfiltered, x="block", y="altexps_ERCC_percent", 
              colour_by="discard") + ggtitle("ERCC percent"),
  nrow=2,
  ncol=2
)