# Single-cell analysis
Single-cell analysis characterize individual cells, allowing their clustering and characterization in an unsupervised manner. The analysis mainly refers to genomic and transcriptomic analysis. The technique can be applied to exploring cell types within a tissue, identify unknown cell types, gene expression change during differentiation processes and so forth. 

# Technical challenges
1. Large volume of data
2. Low depth of sequencing per cell
3. Technical variability accross cells/samples
4. Biological variability accross cells/samples

# Analysis workflow
1. [Best practics for single-cell analysis across modalities](https://www.nature.com/articles/s41576-023-00586-w)
2. Raw count analysis
    1. Filtering low-quality cells and noise correction
    2. Normalization and variance stabilization
    3. Removing confounding sources of variation
    4. Selecting informative features and reducing dimensionality
3. Clustering analysis
    1. Clustering sing cells into groups
    2. Mapping cell clusters to cell identities
    3. From discrete states to continuous processes
4. Mechanism study
    1. Differential gene expression analysis
    2. Gene set enrichment analysis
    3. Compositional analysis
    4. Inferring perturbation effects

![Schematic of a typical scRNA-seq analysis workflow](Images/workflow.png)


# Sample preparation protocols
[credit resource] (https://www.singlecellcourse.org/introduction-to-single-cell-rna-seq.html#overview-of-single-cell-rna-seq)
there are currently a wide diversity of protocols for preparing scRNA-seq data, which can be categorized mainly based on the two most important aspects of cell capture or isolation and transcription quantification.
1. Cell capture
    1. microtitre-plate-based
        1. Cell image and individual cell information before library preparation.
        2. Low-throughput and considerable amounts of work.
    2. microfluidic-array-based
        1. Higher throughput than microtitre-plate-based.
        2. Not appropriate dealing with rare cell types.
    3. microfluidic-droplet-based
        1. The highest throughput. 
2. Transcript quantification
There are two types of transcript quantification: full-length and tag-based. The choice of quantification method has important implications for what types of analyses the data can be used for. 
    1. Full-length
        1. Restricted to plate-based protocols. 
        2. Able to detect splice variants.
    2. Tag-based
        1. One of the ends (3' or 5') of the transcript is sequenced. 
        2. Unique molecular identifiers (UMIs) can improve the quantification accuracy.
        3. Most protocols are tag-based. 

# Read alignment and quantification in droplet-based scRNA-seq data
Most scRNA-seq technologies generate read sequences containing three key pieces of information
1. cDNA fragment
2. Cell barcode
3. Unique molecular identifiers
4. A classical scRNA-seq workflow contains four main steps:
    1. Mapping the cDNA fragments to a reference
    2. Assigning reads to genes
    3. Assigning reads to cells
    4. Counting the number of UMI

# Data precessing and downstream analysis
1. To compute quality control metrics to remove low-quality cells.
2. To normalize count matrix across different cells to eliminate cell-specific biases. 
3. To select a subset of features of interests for downstream analysis.
4. To apply dimensionality reduction analysis. 
5. To cluster cells into groups according to similarities in their normalized expression profiles.  







