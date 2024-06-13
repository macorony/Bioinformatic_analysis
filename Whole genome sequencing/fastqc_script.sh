# Creat directory
mkdir ./WGS_tutorial/quality_control
# Run FastQC
fastqc -o quality_control mini.fastq
# Run MultiQC
multiqc quality_control
