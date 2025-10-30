# LMDseq Usage Guide

## Overview

The LMD_RNAseq pipeline is a specialized bioinformatics workflow designed for RNA sequencing analysis of **low input samples obtained through Laser Microdissection (LMD)**. This pipeline addresses the unique analytical challenges of processing RNA-seq data from microdissected tissue samples, which typically yield very small amounts of RNA and require optimized computational approaches.

## Pipeline Features

- **LMD-Optimized**: Specifically designed for low input RNA samples from laser microdissection
- **Kallisto Quantification**: Fast and accurate pseudo-alignment optimized for LMD samples
- **Quality Control**: Comprehensive FastQC analysis of raw and trimmed reads
- **Adapter Trimming**: TrimGalore for adapter removal and quality trimming
- **Gene Expression Matrix**: Automated creation of merged expression matrices across samples
- **DESeq2 Normalization**: Size factor estimation and variance stabilizing transformation
- **Advanced QC Analysis**: PCA, correlation analysis, and sample distance metrics
- **Comprehensive Reporting**: MultiQC integration for unified quality control reports
- **Reproducible Workflow**: Containerized environment ensuring consistent results


## Input Preparation

### Sample.txt File Format

The pipeline requires a **tab-separated** sample sheet (`Sample.txt`) located in the `assets/` directory. This file contains metadata about your LMD samples and the location of FASTQ files.

#### Required Columns:

| Column | Description | Example | Notes |
|--------|-------------|---------|-------|
| `SID` | Sample ID (unique identifier) | LCM_1, LCM_2 | Must be unique across all samples |
| `SHORT_NAME` | Short sample name for output files | 20942_1D3_A1 | Used in final sample naming |
| `REPLICATE` | Replicate identifier | R1, R2, R3 | Biological or technical replicate |
| `RID` | Run ID from sequencing | S61303 | Sequencing run identifier |
| `strandedness` | Library strandedness | reverse, forward, unstranded | Must match library prep protocol |
| `FASTQ_FOLDER` | Full path to directory containing FASTQ files | /path/to/fastq/files | Must be absolute path |

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
   - Use consistent naming across all samples

2. **Create the Sample.txt file:**
```bash
# Navigate to the assets directory
cd assets/

# Create or edit Sample.txt
nano Sample.txt
```

3. **Fill in the required information:**
   - Use **absolute paths** for `FASTQ_FOLDER`
   - Ensure `strandedness` matches your library preparation protocol
   - Use consistent naming conventions for `SHORT_NAME`
   - Verify all FASTQ directories exist and are accessible

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

### Option 1: Using iGenomes
```bash
nextflow run pdichiaro/LMDseq \
  --input assets/Sample.txt \
  --genome GRCh38 \
  --outdir results \
  -profile docker
```

### Option 2: Custom Reference Files
```bash
nextflow run pdichiaro/LMDseq \
  --input assets/Sample.txt \
  --fasta /path/to/genome.fa \
  --gtf /path/to/annotation.gtf \
  --reference /path/to/reference.txt \
  --outdir results \
  -profile docker
```

### Option 3: Pre-built Kallisto Index
```bash
nextflow run pdichiaro/LMDseq \
  --input assets/Sample.txt \
  --kallisto_index /path/to/kallisto.idx \
  --gtf /path/to/annotation.gtf \
  --outdir results \
  -profile docker
```

### Test Run
```bash
nextflow run pdichiaro/LMDseq \
  -profile test_lmd \
  --outdir test_results
```

## Parameters

### Essential Parameters

#### Input/Output
- `--input`: Path to Sample.txt file (required)
- `--outdir`: Output directory for results (required)

#### Reference Genome
- `--genome`: iGenomes reference (e.g., GRCh38, GRCm38, GRCz10)
- `--fasta`: Custom genome FASTA file
- `--gtf`: Custom GTF annotation file
- `--kallisto_index`: Pre-built Kallisto index
- `--reference`: Reference annotation file for enhanced gene expression matrix

### LMD-Specific Parameters

#### Kallisto Quantification
- `--aligner`: Alignment method (default: `kallisto`)
- `--extra_kallisto_quant_args`: Additional Kallisto arguments (default: `-b 50 --single-overhang --genomebam`)
- `--extra_kallisto_index_args`: Additional Kallisto index arguments
- `--kallisto_quant_fraglen`: Fragment length for single-end reads (default: 200)
- `--kallisto_quant_fraglen_sd`: Fragment length standard deviation (default: 200)
- `--kallisto_kmer_size`: K-mer size for Kallisto index (default: 31)

#### DESeq2 Analysis
- `--min_reads`: Minimum read count threshold for gene filtering (default: 10)

#### Read Processing
- `--trimmer`: Trimming tool (default: `trimgalore`)
- `--extra_trimgalore_args`: Additional TrimGalore arguments (default: `--clip_R2 3 --quality 20 --stringency 3 --length 20`)

### Process Control

#### Quality Control
- `--skip_fastqc`: Skip FastQC analysis
- `--skip_multiqc`: Skip MultiQC report generation

#### Process Skipping
- `--skip_trimming`: Skip adapter trimming step
- `--skip_alignment`: Skip alignment processes
- `--skip_gtf_filter`: Skip GTF filtering (default: true)


## Output Structure

The pipeline generates a structured output directory with numbered folders for easy navigation:

```
results/
├── 1_FASTQ/                    # Merged FASTQ files
├── 2_QC/                       # FastQC quality control reports
├── 3_TRIM/                     # Adapter trimming outputs and reports
├── 4_BAM/                      # Kallisto quantification results and pseudo-BAM files
├── 5_Bw/                       # BigWig coverage files for genome browser visualization
├── 6_Norm_folder/              # Normalization plots and read distribution analysis
│   └── Read_Distribution/      # Before/after normalization distribution plots
├── 7_Counts_folder/            # Gene expression matrices (raw and normalized)
├── 8_Quality_folder/           # Advanced QC analysis (PCA, correlations, distances)
├── multiqc/                    # Integrated quality control reports
└── pipeline_info/              # Pipeline execution information and resource usage
```

### Key Output Files

#### Gene Expression Data
- **Raw Expression Matrix**: `7_Counts_folder/EX_reads_RAW.txt`
- **Normalized Expression Matrix**: `7_Counts_folder/EX_reads_NORM_filt.txt`
- **VST-transformed Data**: `7_Counts_folder/rLog_reads.txt`
- **Normalization Parameters**: `7_Counts_folder/Normalisation_Parameters.txt`

#### Quality Control Reports
- **Comprehensive QC Report**: `multiqc/multiqc_report.html`
- **PCA Analysis**: `8_Quality_folder/PCA_rlogTransformedID.pdf`
- **Sample Correlations**: `8_Quality_folder/Heatmap_sampleTosample_correlation.pdf`
- **Distance Analysis**: `8_Quality_folder/Heatmap_sampleTosample_distances_vstTransformed.pdf`
- **Read Distribution**: `6_Norm_folder/Read_Distribution/`

#### Individual Sample Results
- **Kallisto Quantification**: `4_BAM/*/abundance.tsv`
- **Quality Control**: `2_QC/*_fastqc.html`
- **Trimming Reports**: `3_TRIM/*_trimming_report.txt`
- **Coverage Files**: `5_Bw/*.bigwig`


## Support and Citation

### Getting Help
- **Issues**: [GitHub Issues](https://github.com/pdichiaro/LMD_rnaseq/issues)
- **Documentation**: Complete documentation in `docs/` directory
- **Contact**: Pierluigi Di Chiaro (AUSL-IRCCS Reggio Emilia)


### Citation
If you use LMD_RNAseq in your research, please cite:

> **Di Chiaro P, et al.** "Mapping functional to morphological variation reveals the basis of regional extracellular matrix subversion and nerve invasion in pancreatic cancer." *Cancer Cell* 2024.
> **Di Chiaro P, et al.** "A framework to mine laser microdissection-based omics data and uncover regulators of pancreatic cancer heterogeneity." *Gigascience* 2025.
