# Raw data to count matrix 
## Library preparation
### Amplification method
#### Pooled PCR amplification
1. Drop-seq
2. InDrop
3. CEL-seq
4. MARS-seq
5. SCRB-seq
6. CEL-seq2
#### Individual cell amplification
1. SMART-seq
2. SMART-seq2
3. STRT-seq
## Generation of count matrix using Alevin
[credit resource] (https://combine-lab.github.io/alevin-tutorial/2018/setting-up-resources/)
Example data
```bash
wget https://cg.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_fastqs.tar

tar -xvf pbmc4k_fastqs.tar
```

Reference transcriptome
```bash
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.pc_transcripts.fa.gz
```

Index transcriptome
```bash
salmon index -i index -k 31 --gencode -p 4 -t gencode.v31.pc_transcripts.fa.gz
```

Transcript to gene mapping
```bash
bioawk -c gff '$feature=="transcript" {print $group}' <(gunzip -c gencode.v31.primary_assembly.annotation.gtf.gz) | awk -F ' ' '{print substr($4,2,length($4)-3) "\t" substr($2,2,length($2)-3)}' - > txp2gene.tsv
```

Running Alevin
```bash
salmon alevin -lISR -1 fastqs/pbmc4k_S1_L001_R1_001.fastq.gz -2 fastqs/pbmc4k_S1_L001_R2_001.fastq.gz --chromium -i index -p 8 -o alevin_output1 --tgMap txp2gene.tsv
```

```bash
salmon alevin -lISR -1 fastqs/pbmc4k_S1_L002_R1_001.fastq.gz -2 fastqs/pbmc4k_S1_L002_R2_001.fastq.gz --chromium -i index -p 8 -o alevin_output2 --tgMap txp2gene.tsv
```
## About Salmon
[credit resource] (https://salmon.readthedocs.io/en/latest/salmon.html)
Salmon is a tool for wicked-fast transcript quantification from RNA-seq data. Salmon supports dual modes in quantification: The mapping-based mode and the alignment-based mode. 
The mapping-based mode of salmon runs in two phases: indexing and quantification. 
The alignment-based mode does not require indexing. The quantification only needs a FASTA file of transcripts and a SAM/BAM files containing the alignments. 








