# Integration Summary: create_reference_db.R

## Overview
The `create_reference_db.R` script from the `bin/` folder has been successfully integrated into the LMD_rnaseq pipeline to generate a merged gene expression matrix after the Kallisto subworkflow.

## Integration Details

### 1. R Script Fixes
- **Fixed undefined `samples` variable**: Added `samples <- basename(BAM_DIR)` to extract sample names from directory paths
- **Script location**: `bin/create_reference_db.R`
- **Dependencies**: Uses optparse, tximport, and readr R packages

### 2. TXIMPORT Module
- **Location**: `modules/local/tximport/main.nf`
- **Purpose**: Wraps the R script as a Nextflow process
- **Inputs**: 
  - `all_folders`: Collection of Kallisto result directories
  - `gtf`: GTF annotation file
  - `reference`: Reference transcript-to-gene mapping file
- **Output**: `EX_reads_RAW.txt` - merged gene expression matrix
- **Environment**: Uses conda environment with R and required packages

### 3. Kallisto Subworkflow Integration
- **Location**: `subworkflows/local/kallisto/main.nf`
- **Integration point**: After KALLISTO_QUANT process
- **Channel handling**: Uses `.collect()` to gather all Kallisto results before passing to TXIMPORT
- **Output**: Emits the merged matrix as `matrix` channel

### 4. Main Workflow Integration
- **Location**: `workflows/lmd_rnaseq.nf`
- **Usage**: Calls KALLISTO subworkflow and captures matrix output
- **Variable**: `ch_raw_matrix = KALLISTO.out.matrix`
- **Emit**: Added `raw_matrix` to the workflow emit section for downstream use

## Workflow Flow

```
FASTQ Files
    ↓
FASTQ_FASTQC_TRIMGALORE
    ↓
KALLISTO Subworkflow:
    ├── KALLISTO_INDEX (if needed)
    ├── KALLISTO_QUANT (per sample)
    └── TXIMPORT (merge all samples)
        └── create_reference_db.R
            └── EX_reads_RAW.txt
```

## Output Files

### Generated Matrix File
- **Filename**: `EX_reads_RAW.txt`
- **Format**: Tab-separated values
- **Content**: Gene expression counts matrix
- **Rows**: Genes (from GTF annotation)
- **Columns**: Samples (from input Sample.txt)

### Location in Results
The merged gene expression matrix will be available in:
```
results/
└── tximport/
    └── EX_reads_RAW.txt
```

## Usage

The integration is automatic when running the pipeline. The merged gene expression matrix is generated after all individual sample quantifications are complete.

### Example Command
```bash
nextflow run . \
  --input assets/Sample.txt \
  --genome GRCh38 \
  --outdir results \
  -profile docker
```

### Required Files
1. **Sample.txt**: Sample metadata file in assets/
2. **GTF file**: Gene annotation (via --gtf or --genome)
3. **Reference file**: Transcript-to-gene mapping (generated automatically)

## Technical Notes

### R Script Parameters
- `-g, --gtf`: GTF annotation file path
- `-r, --reference`: Reference transcript-to-gene mapping file
- `-i, --input`: Space-separated list of Kallisto result directories
- `-o, --output`: Output matrix filename (default: EX_reads_RAW.txt)

### Dependencies
- R packages: tximport, readr, optparse
- Conda environment: `modules/local/tximport/environment.yml`
- Docker container: `pdichiaro/lmdseq:latest`

### Channel Handling
- Individual Kallisto results are collected using `.collect()`
- All sample directories are passed as a single input to TXIMPORT
- The R script processes all samples simultaneously to create the merged matrix

## Verification

The integration has been verified by:
1. ✅ Syntax checking the R script
2. ✅ Validating Nextflow workflow syntax
3. ✅ Confirming proper channel flow
4. ✅ Testing pipeline help functionality

## Future Enhancements

Potential improvements for the integration:
1. Add parameter to customize output matrix filename
2. Include additional matrix formats (e.g., TPM, FPKM)
3. Add quality control metrics for the merged matrix
4. Generate summary statistics for the expression data

---

*Integration completed on: $(date)*
*Pipeline version: 0.0.1*