# LMD_rnaseq Usage Guide

## Overview

The LMD_rnaseq pipeline is designed for RNA sequencing analysis of low input samples, particularly those obtained through Laser Microdissection (LMD). This pipeline performs quality control, trimming, pseudo-alignment with Kallisto, and generates comprehensive reports using MultiQC.

## Pipeline Features

- **Quality Control**: FastQC analysis of raw and trimmed reads
- **Adapter Trimming**: TrimGalore for adapter removal and quality trimming
- **Pseudo-alignment**: Kallisto for fast and accurate quantification
- **DESeq2 Normalization**: Automated gene expression normalization and filtering
- **Advanced QC Analysis**: Sample correlation, PCA, and distance analysis
- **Comprehensive Reporting**: MultiQC for integrated quality control reports
- **Low Input Optimization**: Specifically designed for LMD samples with limited RNA

## Prerequisites

### Software Requirements
- Nextflow (≥21.10.3)
- Docker, Singularity, or Conda for dependency management
- Git for cloning the repository

### System Requirements
- Minimum 8GB RAM
- 50GB free disk space (depending on dataset size)
- Multi-core processor recommended

## Installation

1. **Clone the repository:**
```bash
git clone https://github.com/pdichiaro/LMD_rnaseq.git
cd LMD_rnaseq
```

2. **Install Nextflow (if not already installed):**
```bash
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

## Input Preparation

### Sample.txt File Format

The pipeline requires a tab-separated sample sheet (`Sample.txt`) located in the `assets/` directory. This file contains metadata about your samples and the location of FASTQ files.

#### Required Columns:

| Column | Description | Example |
|--------|-------------|---------|
| `SID` | Sample ID (unique identifier) | LCM_1, LCM_2 |
| `SHORT_NAME` | Short sample name for output files | 20942_1D3_A1 |
| `REPLICATE` | Replicate identifier | R1, R2, R3 |
| `RID` | Run ID from sequencing | S61303 |
| `strandedness` | Library strandedness | reverse, forward, unstranded |
| `FASTQ_FOLDER` | Full path to directory containing FASTQ files | /path/to/fastq/files |

#### Example Sample.txt:
```
SID	SHORT_NAME	REPLICATE	RID	strandedness	FASTQ_FOLDER
LCM_1	20942_1D3_A1	R1	S61303	reverse	/data/fastq/Sample_S61303_LCM_1
LCM_2	20942_1D3_B1	R2	S61304	reverse	/data/fastq/Sample_S61304_LCM_2
LCM_3	20942_1D3_C1	R3	S61305	reverse	/data/fastq/Sample_S61305_LCM_3
```

### Creating Your Sample.txt File

1. **Prepare your data structure:**
   - Organize FASTQ files in separate directories for each sample
   - Ensure FASTQ files follow naming convention: `*_R1_*.fastq.gz` and `*_R2_*.fastq.gz` for paired-end data

2. **Create the Sample.txt file:**
```bash
# Navigate to the assets directory
cd assets/

# Create or edit Sample.txt
nano Sample.txt
```

3. **Fill in the required information:**
   - Use absolute paths for `FASTQ_FOLDER`
   - Ensure strandedness matches your library preparation protocol
   - Use consistent naming conventions

### FASTQ File Organization

Your FASTQ files should be organized as follows:
```
/path/to/data/
├── Sample_S61303_LCM_1/
│   ├── sample1_R1_001.fastq.gz
│   └── sample1_R2_001.fastq.gz
├── Sample_S61304_LCM_2/
│   ├── sample2_R1_001.fastq.gz
│   └── sample2_R2_001.fastq.gz
└── Sample_S61305_LCM_3/
    ├── sample3_R1_001.fastq.gz
    └── sample3_R2_001.fastq.gz
```

## Reference Genome Setup

### Option 1: Using iGenomes (Recommended)
```bash
nextflow run . --genome GRCh38 --input assets/Sample.txt --outdir results
```

### Option 2: Custom Reference Files
```bash
nextflow run . \
  --input assets/Sample.txt \
  --fasta /path/to/genome.fa \
  --gtf /path/to/annotation.gtf \
  --reference /path/to/reference.txt \
  --outdir results
```

### Option 3: Pre-built Kallisto Index
```bash
nextflow run . \
  --input assets/Sample.txt \
  --kallisto_index /path/to/kallisto.idx \
  --gtf /path/to/annotation.gtf \
  --outdir results
```

## Running the Pipeline

### Basic Command
```bash
nextflow run . \
  --input assets/Sample.txt \
  --genome GRCh38 \
  --outdir results \
  -profile docker
```

### Advanced Configuration
```bash
nextflow run . \
  --input assets/Sample.txt \
  --genome GRCh38 \
  --reference /path/to/reference.txt \
  --outdir results \
  --aligner kallisto \
  --trimmer trimgalore \
  --extra_trimgalore_args '--clip_R2 3 --quality 20 --stringency 3 --length 20' \
  --kallisto_quant_fraglen 200 \
  --kallisto_quant_fraglen_sd 200 \
  --max_cpus 16 \
  --max_memory '32.GB' \
  -profile docker \
  -resume
```

### Profile Options

Choose the appropriate profile for your system:

- **Docker**: `-profile docker` (recommended for most users)
- **Singularity**: `-profile singularity` (for HPC environments)
- **Conda**: `-profile conda` (if containers are not available)

### Key Parameters

#### Input/Output
- `--input`: Path to Sample.txt file (default: `assets/Sample.txt`)
- `--outdir`: Output directory for results (required)

#### Reference Genome
- `--genome`: iGenomes reference (e.g., GRCh38, GRCm38)
- `--fasta`: Custom genome FASTA file
- `--gtf`: Custom GTF annotation file
- `--kallisto_index`: Pre-built Kallisto index
- `--reference`: Reference annotation file for gene expression matrix creation

#### Trimming Options
- `--trimmer`: Trimming tool (default: `trimgalore`)
- `--extra_trimgalore_args`: Additional TrimGalore arguments
- `--skip_trimming`: Skip adapter trimming step

#### Alignment Options
- `--aligner`: Alignment method (default: `kallisto`)
- `--kallisto_quant_fraglen`: Fragment length for single-end reads (default: 200)
- `--kallisto_quant_fraglen_sd`: Fragment length standard deviation (default: 200)

#### Quality Control
- `--skip_fastqc`: Skip FastQC analysis
- `--skip_multiqc`: Skip MultiQC report generation

#### DESeq2 Normalization (Advanced)
- `--min_reads`: Minimum read count threshold for gene filtering (default: 10)

#### Resource Limits
- `--max_cpus`: Maximum number of CPUs (default: 10)
- `--max_memory`: Maximum memory usage (default: 16.GB)
- `--max_time`: Maximum runtime (default: 240.h)

## Output Structure

The pipeline generates the following output structure:

```
results/
├── pipeline_info/          # Pipeline execution reports
├── fastqc/                 # FastQC reports (raw reads)
├── trimgalore/            # TrimGalore outputs
├── fastqc_trim/           # FastQC reports (trimmed reads)
├── kallisto/              # Kallisto quantification results
├── 6_Norm_folder/         # Normalization plots and analysis
│   └── Read_Distribution/ # Before/after normalization distribution plots
├── 7_Counts_folder/       # Raw and normalized gene expression matrices
├── 8_Quality_folder/      # Advanced QC analysis (PCA, correlations, distances)
├── multiqc/               # MultiQC integrated report
└── multiqc_data/          # MultiQC data files
```

### Merged Gene Expression Matrix

The pipeline automatically creates a merged gene expression matrix (`EX_reads_RAW.txt`) in the `7_Counts_folder/` directory. This matrix contains:

- **Rows**: Genes/transcripts from the reference annotation
- **Columns**: Sample names (based on SHORT_NAME_REPLICATE from Sample.txt)
- **Values**: Raw expression counts from Kallisto quantification
- **Format**: Tab-separated text file suitable for downstream analysis (DESeq2, edgeR, etc.)

This merged matrix is created using the `create_reference_db.R` script, which combines individual Kallisto quantification results into a single expression matrix for all samples.

### DESeq2 Normalization and Quality Control

After creating the raw gene expression matrix, the pipeline automatically performs:

1. **Gene Filtering**: 
   - Filters genes by type (protein-coding genes by default)
   - Removes low-expressed genes (< 10 reads by default)
   - Keeps genes expressed in at least half of the samples

2. **DESeq2 Normalization**:
   - Size factor estimation and normalization
   - Variance stabilizing transformation (VST) for visualization

3. **Quality Control Analysis**:
   - **Read Distribution Plots**: Before and after normalization
   - **Sample Distance Analysis**: Euclidean and Poisson distance heatmaps
   - **Correlation Analysis**: Sample-to-sample correlation heatmaps
   - **Principal Component Analysis**: 2D and 3D PCA plots
   - **Variance Analysis**: PC variance contribution plots

4. **Output Files**:
   - `7_Counts_folder/EX_reads_NORM_filt.txt`: DESeq2 normalized expression matrix
   - `7_Counts_folder/rLog_reads.txt`: Variance-stabilized data for visualization
   - `7_Counts_folder/Normalisation_Parameters.txt`: Size factors used for normalization
   - Various QC plots in PDF and HTML formats in `8_Quality_folder/`

### Key Output Files

- **MultiQC Report**: `results/multiqc/multiqc_report.html` - Comprehensive quality control report
- **Gene Counts**: `results/kallisto/*/abundance.tsv` - Individual sample gene expression quantification
- **Raw Gene Expression Matrix**: `results/7_Counts_folder/EX_reads_RAW.txt` - Combined raw gene expression matrix
- **Normalized Gene Expression Matrix**: `results/7_Counts_folder/EX_reads_NORM_filt.txt` - DESeq2 normalized and filtered matrix
- **rlog Transformed Matrix**: `results/7_Counts_folder/rLog_reads.txt` - Variance-stabilized expression data for visualization
- **QC Analysis**: `results/8_Quality_folder/` - PCA plots, correlation heatmaps, and distance analysis
- **Normalization Plots**: `results/6_Norm_folder/Read_Distribution/` - Before/after normalization distribution plots
- **Pipeline Report**: `results/pipeline_info/` - Execution statistics and resource usage

## Example Workflow

### 1. Prepare Your Environment
```bash
# Clone the repository
git clone https://github.com/pdichiaro/LMD_rnaseq.git
cd LMD_rnaseq

# Verify Nextflow installation
nextflow -version
```

### 2. Prepare Sample Sheet
```bash
# Edit the sample sheet
nano assets/Sample.txt

# Example content:
SID	SHORT_NAME	REPLICATE	RID	strandedness	FASTQ_FOLDER
LMD_1	Sample_A	R1	S001	reverse	/data/fastq/Sample_A
LMD_2	Sample_B	R2	S002	reverse	/data/fastq/Sample_B
```

### 3. Run the Pipeline
```bash
# Basic run with human genome
nextflow run . \
  --input assets/Sample.txt \
  --genome GRCh38 \
  --outdir results \
  -profile docker

# Monitor progress
tail -f .nextflow.log
```

### 4. Review Results
```bash
# Open MultiQC report
firefox results/multiqc/multiqc_report.html

# Check pipeline execution report
firefox results/pipeline_info/execution_report.html
```

## Troubleshooting

### Common Issues

#### 1. Sample.txt Format Errors
**Error**: `No such file or directory: Sample.txt`
**Solution**: Ensure Sample.txt is in the assets/ directory with correct tab-separated format

#### 2. FASTQ File Not Found
**Error**: `No such file or directory: /path/to/fastq`
**Solution**: Verify FASTQ_FOLDER paths are absolute and accessible

#### 3. Memory Issues
**Error**: `Process exceeded available memory`
**Solution**: Increase memory limits or reduce parallel processes:
```bash
--max_memory '32.GB' --max_cpus 8
```

#### 4. Reference Genome Issues
**Error**: `Genome 'XXX' not found`
**Solution**: Use supported genomes or provide custom reference files:
```bash
--fasta genome.fa --gtf annotation.gtf
```

### Performance Optimization

#### For Large Datasets
```bash
nextflow run . \
  --input assets/Sample.txt \
  --genome GRCh38 \
  --outdir results \
  --max_cpus 32 \
  --max_memory '64.GB' \
  -profile docker \
  -resume
```

#### For Small/Test Datasets
```bash
nextflow run . \
  --input assets/Sample.txt \
  --genome GRCh38 \
  --outdir results \
  --max_cpus 4 \
  --max_memory '8.GB' \
  -profile docker
```

## Advanced Usage

### Custom Configuration

Create a custom configuration file (`custom.config`):
```nextflow
params {
    // Custom trimming parameters
    extra_trimgalore_args = '--clip_R2 5 --quality 25 --stringency 5 --length 25'
    
    // Kallisto-specific parameters
    kallisto_quant_fraglen = 180
    kallisto_quant_fraglen_sd = 20
    
    // Resource allocation
    max_cpus = 24
    max_memory = '48.GB'
}

process {
    // Increase memory for specific processes
    withName: KALLISTO_QUANT {
        memory = '16.GB'
    }
}
```

Run with custom configuration:
```bash
nextflow run . -c custom.config --input assets/Sample.txt --genome GRCh38 --outdir results
```

### Resume Failed Runs

If a pipeline run fails or is interrupted, you can resume from the last successful step:
```bash
nextflow run . \
  --input assets/Sample.txt \
  --genome GRCh38 \
  --outdir results \
  -profile docker \
  -resume
```

## Support and Citation

### Getting Help
- Check the [GitHub Issues](https://github.com/pdichiaro/LMD_rnaseq/issues) for common problems
- Review the Nextflow documentation for general pipeline questions
- Contact the pipeline author for specific LMD_rnaseq questions

### Citation
If you use this pipeline in your research, please cite:
- The LMD_rnaseq pipeline: [GitHub repository](https://github.com/pdichiaro/LMD_rnaseq)
- Nextflow: Di Tommaso et al. (2017) Nature Biotechnology
- Kallisto: Bray et al. (2016) Nature Biotechnology
- MultiQC: Ewels et al. (2016) Bioinformatics

---

*This documentation was generated for LMD_rnaseq pipeline version 0.0.1*