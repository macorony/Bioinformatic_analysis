# Raw data to count matrix 
_**Resource credit:**_
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
Example data
```bash
wget https://cg.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_fastqs.tar

tar -xvf pbmc4k_fastqs.tar
```

Reference transcriptome
```bash
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.pc_transcripts.fa.gz

salmon index -i index -k 31 --gencode -p 4 -t gencode.v31.pc_transcripts.fa.gz

```



