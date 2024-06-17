# Allocate resource
salloc --time=05:00:00 -mem-per-cpu=8G --ntasks=4
# Creat directory
mkdir ./WGS_tutorial/quality_control
# Run FastQC
fastqc -o quality_control mini.fastq
# Run MultiQC
multiqc quality_control


# De novo assembly
spades.py -o assembly/spades/ --careful -1 data/evol1_R1.fastq.gz data/evol1_R2.fastq.gz

# BWA index
bwa index assembly/scaffolds.fasta

# Mapping 
bwa mem assembly/scaffolds.fasta data/evol1_R1.fastq.gz data/evol1_R2.fastq.gz > mapping/evol1.sam






# fixmate and compress to bam
samtools sort -n -O sam mappings/evol1.sam | samtools fixmate -m -O bam - mappings/evol1.fixmate.bam
rm mappings/evol1.sam
# sort
samtools sort -O bam -o mappings/evol1.sorted.bam mappings/evol1.fixmate.bam
rm mappings/evol1.fixmate.bam
# mark duplicates
samtools markdup -r -S mappings/evol1.sorted.bam mappings/evol1.sorted.dedup.bam
rm mappings/evol1.sorted.bam
# extract q20 mappers
samtools view -h -b -q 20 mappings/evol1.sorted.dedup.bam > mappings/evol1.sorted.dedup.q20.bam
# extract unmapped
samtools view -b -f 4 mappings/evol1.sorted.dedup.bam > mappings/evol1.sorted.unmapped.bam
rm mappings/evol1.sorted.dedup.bam
# covert to fastq
samtools fastq -1 mappings/evol1.sorted.unmapped.R1.fastq.gz -2 mappings/evol1.sorted.unmapped.R2.fastq.gz mappings/evol1.sorted.unmapped.bam
# delete not needed files
rm mappings/evol1.sorted.unmapped.bam

#
# Evol 2
#

samtools sort -n -O sam mappings/evol2.sam | samtools fixmate -m -O bam - mappings/evol2.fixmate.bam
rm mappings/evol2.sam
samtools sort -O bam -o mappings/evol2.sorted.bam mappings/evol2.fixmate.bam
rm mappings/evol2.fixmate.bam
samtools markdup -r -S mappings/evol2.sorted.bam mappings/evol2.sorted.dedup.bam
rm mappings/evol2.sorted.bam
samtools view -h -b -q 20 mappings/evol2.sorted.dedup.bam > mappings/evol2.sorted.dedup.q20.bam
rm mappings/evol2.sorted.dedup.bam