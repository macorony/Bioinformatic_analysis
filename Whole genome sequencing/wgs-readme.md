# Whole Genome Sequencing Analysis Pipeline

A Snakemake-based pipeline for analyzing whole genome sequencing data, including quality control, de novo assembly, mapping, and variant calling.

## Features

- Quality control and read trimming with FastP
- De novo assembly using SPAdes
- Read mapping with BWA-MEM
- Variant calling using FreeBayes
- Automated workflow management with Snakemake
- Configurable parameters through YAML files
- Comprehensive quality reporting with MultiQC

## Prerequisites

- Conda or Mamba package manager
- Git
- 8GB+ RAM (recommended: 16GB+)
- 50GB+ free disk space

## Installation

1. Clone the repository:
```bash
git clone https://github.com/yourusername/wgs-pipeline.git
cd wgs-pipeline
```

2. Create the conda environment:
```bash
conda env create -f envs/environment.yaml
conda activate wgs-pipeline
```

## Usage

1. Configure your samples:
   - Add your sample information to `config/samples.tsv`
   - Adjust parameters in `config/config.yaml`

2. Prepare your data:
   - Place your raw FASTQ files in the `data/` directory
   - Files should be named as: `{sample}_R1.fastq.gz` and `{sample}_R2.fastq.gz`

3. Run the pipeline:
```bash
# Dry run to check workflow
snakemake -n

# Run the pipeline using 4 cores
snakemake --use-conda -j 4

# Run with cluster submission (e.g., SLURM)
snakemake --cluster "sbatch --mem=8G" -j 100
```

## Pipeline Steps

1. **Quality Control**
   - Adapter trimming
   - Quality trimming
   - Overrepresentation analysis
   - Quality reports generation

2. **De Novo Assembly**
   - SPAdes assembly with --careful mode
   - Scaffold generation

3. **Mapping**
   - BWA index creation
   - Read mapping to assembly
   - SAM to BAM conversion
   - Duplicate marking
   - Quality filtering

4. **Variant Calling**
   - FreeBayes variant calling
   - Variant filtering
   - VCF compression and indexing

## Output Structure

```
results/
├── trimmed/           # Trimmed FASTQ files
├── fastqc/            # FastQC reports
├── assembly/          # De novo assembly
├── mappings/          # BWA alignments
└── variants/          # Called variants
```

## Configuration

Edit `config/config.yaml` to modify:
- Resource allocation
- Tool parameters
- Quality thresholds
- Output paths

## Contributing

1. Fork the repository
2. Create a feature branch
3. Commit your changes
4. Push to the branch
5. Create a Pull Request

## Troubleshooting

Common issues and solutions:

1. **Insufficient memory**:
   - Adjust memory in config.yaml
   - Use SLURM or other cluster system

2. **Missing files**:
   - Check input file naming
   - Verify sample sheet format

3. **Failed jobs**:
   - Check log files in `logs/`
   - Verify input data quality

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contact

Your Name - your.email@example.com
Project Link: https://github.com/yourusername/wgs-pipeline

## Acknowledgments

- FastP
- SPAdes
- BWA
- FreeBayes
- Snakemake
