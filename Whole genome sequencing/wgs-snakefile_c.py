# Snakefile
configfile: "config/config.yaml"

# Load sample information
import pandas as pd
samples = pd.read_table("config/samples.tsv").set_index("sample", drop=False)

# Final output files
rule all:
    input:
        expand("results/variants/{sample}.freebayes.filtered.vcf.gz", sample=samples.index),
        "results/qc/multiqc_report.html"

# Include rules
include: "rules/qc.smk"
include: "rules/assembly.smk"
include: "rules/mapping.smk"
include: "rules/variants.smk"

# Quality control rule
rule fastp:
    input:
        r1 = "data/{sample}_R1.fastq.gz",
        r2 = "data/{sample}_R2.fastq.gz"
    output:
        r1 = "results/trimmed/{sample}_R1.fastq.gz",
        r2 = "results/trimmed/{sample}_R2.fastq.gz",
        html = "results/trimmed/{sample}.fastp.html",
        json = "results/trimmed/{sample}.fastp.json"
    threads: config["fastp"]["threads"]
    log: "logs/fastp/{sample}.log"
    shell:
        """
        fastp --detect_adapter_for_pe \
            --overrepresentation_analysis \
            --correction --cut_right \
            --thread {threads} \
            --html {output.html} \
            --json {output.json} \
            -i {input.r1} -I {input.r2} \
            -o {output.r1} -O {output.r2} \
            2> {log}
        """

# Add more rules for each step...
