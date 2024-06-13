# Allocate resource
salloc --time=05:00:00 -mem-per-cpu=8G --ntasks=4
# Creat directory
mkdir ./WGS_tutorial/quality_control
# Run FastQC
fastqc -o quality_control mini.fastq
# Run MultiQC
multiqc quality_control
