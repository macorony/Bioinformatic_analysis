# Overview of whole genome sequencing analysis
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
[the reference link](https://frontlinegenomics.com/how-to-ngs-quality-control/)
### Alignment
### Post-processing Analysis
### Variant Calling and Filtering