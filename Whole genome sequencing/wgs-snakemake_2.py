rule all:
    input:
        expand("variants/{sample}.freebayes.filtered.vcf.gz", sample=["evol1", "evol2"])

rule fastp:
    input:
        r1="data/{sample}_R1.fastq.gz",
        r2="data/{sample}_R2.fastq.gz"
    output:
        r1_trim="trimmed/{sample}_R1.fastq.gz",
        r2_trim="trimmed/{sample}_R2.fastq.gz",
        html="trimmed/{sample}.fastp.html",
        json="trimmed/{sample}.fastp.json"
    params:
        threads=2
    shell:
        """
        fastp --detect_adapter_for_pe --overrepresentation_analysis --correction \
              --cut_right --thread {params.threads} \
              -i {input.r1} -I {input.r2} \
              -o {output.r1_trim} -O {output.r2_trim} \
              --html {output.html} --json {output.json}
        """

rule fastqc:
    input:
        "trimmed/{sample}_R1.fastq.gz",
        "trimmed/{sample}_R2.fastq.gz"
    output:
        "trimmed-fastqc/{sample}_R1_fastqc.html",
        "trimmed-fastqc/{sample}_R2_fastqc.html"
    shell:
        "fastqc -o trimmed-fastqc {input}"

rule mapping:
    input:
        r1="trimmed/{sample}_R1.fastq.gz",
        r2="trimmed/{sample}_R2.fastq.gz",
        ref="assembly/scaffolds.fasta"
    output:
        "mappings/{sample}.sam"
    shell:
        "bwa mem {input.ref} {input.r1} {input.r2} > {output}"
