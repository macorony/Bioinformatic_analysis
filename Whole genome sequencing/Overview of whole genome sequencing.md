# Overview of whole genome sequencing analysis 
(https://www.basepairtech.com/knowledge-center/whole-genome-sequencing-analysis-an-overview/)
The WGS analysis is to identify variants by mapping the short fragment reads to a reference genome. 
## Analysis Pipeline
### Quality Control
#### Tools
1. fastqc
2. multiqc
#### Metrics
1. Yield: The total number of reads per run
2. Read analysis: Read length, GC content, adapter content and duplication
3. Error rate: The percentage of bases incorrectly called. As read length increases, error rate also increase. 
4. Clusters passing filter: An indiction of signal purity. 
5. Phas/Prephas

[The link](https://frontlinegenomics.com/how-to-ngs-quality-control/)

#### The QC process
1. Remove Phix sequences (what is PhiX sequence: https://www.illumina.com/products/by-type/sequencing-kits/cluster-gen-sequencing-reagents/phix-control-v3.html)

2. Adapter trimming
3. Quality trimming of reads
4. Quality assessment

### De novo assembly
De novo genome assembly is a strategy for genome assembly, representing the genome assembly of a novel genome from scratch without the aid of reference genomic data.  
### Alignment
When the reference genome is known, the alignment of short reads to the reference genome requires a genome indexing step aiming to reduce and improve computational efficiency during the mapping process. Next, the reads are mapped to the reference sequence using BWA.  
### Post-processing Analysis
After obtaining the SAM or BAM files, reads that are uniquely mapped to the reference genome must be sorted and filtered, which minimizes errors during variant calling. Samtools, sambamba and Picard are software tools widely used to manipulated BAM and SAM files. 
### Variant Calling and Filtering
This step aims to identify the polymorphic regions in the DNA of a sample. The algorithms are based on the likelihood that a given variant (SNV or INDEL) exist in a position of BAM file. The identified variants are stored in a VCF file format. The next step is to filter the variants to only retains those that meeting the minimum quality criteria required, such as base quality and depth. Finally, an additional annotation step can be used to integrate information and imporve variant filtration and prioritization.