# Folder Structure Changes Summary

## Overview
Removed unnecessary nested folder structures from the DESeq2 normalization and QC outputs to simplify the pipeline output organization.

## Changes Made

### 1. Removed `9_Stat_folder`
- **Before**: Pipeline created `9_Stat_folder/` directory
- **After**: No `9_Stat_folder` is created
- **Impact**: Eliminates unused statistical output directory

### 2. Removed `ALL` subfolder in `7_Counts_folder`
- **Before**: Files saved to `7_Counts_folder/ALL/`
  - `7_Counts_folder/ALL/EX_reads_RAW_filt.txt`
  - `7_Counts_folder/ALL/EX_reads_NORM_filt.txt`
  - `7_Counts_folder/ALL/rLog_reads.txt`
  - `7_Counts_folder/ALL/Normalisation_Parameters.txt`
- **After**: Files saved directly to `7_Counts_folder/`
  - `7_Counts_folder/EX_reads_RAW_filt.txt`
  - `7_Counts_folder/EX_reads_NORM_filt.txt`
  - `7_Counts_folder/rLog_reads.txt`
  - `7_Counts_folder/Normalisation_Parameters.txt`

### 3. Removed `ALL_CONDITIONS/ALL` subfolders in `8_Quality_folder`
- **Before**: QC files saved to `8_Quality_folder/ALL_CONDITIONS/ALL/`
  - `8_Quality_folder/ALL_CONDITIONS/ALL/Heatmap_sampleTosample_distances_vstTransformed.pdf`
  - `8_Quality_folder/ALL_CONDITIONS/ALL/Heatmap_Poisson_Distance.pdf`
  - `8_Quality_folder/ALL_CONDITIONS/ALL/Heatmap_sampleTosample_correlation.pdf`
  - `8_Quality_folder/ALL_CONDITIONS/ALL/PCA_rlogTransformedID.pdf`
  - `8_Quality_folder/ALL_CONDITIONS/ALL/Percent_var_PCA.pdf`
  - `8_Quality_folder/ALL_CONDITIONS/ALL/All_pca.pdf`
  - `8_Quality_folder/ALL_CONDITIONS/ALL/3D_PCA_*.html`
- **After**: QC files saved directly to `8_Quality_folder/`
  - `8_Quality_folder/Heatmap_sampleTosample_distances_vstTransformed.pdf`
  - `8_Quality_folder/Heatmap_Poisson_Distance.pdf`
  - `8_Quality_folder/Heatmap_sampleTosample_correlation.pdf`
  - `8_Quality_folder/PCA_rlogTransformedID.pdf`
  - `8_Quality_folder/Percent_var_PCA.pdf`
  - `8_Quality_folder/All_pca.pdf`
  - `8_Quality_folder/3D_PCA_*.html`

### 4. Simplified `6_Norm_folder` structure
- **Before**: Files saved to `6_Norm_folder/ALL/Read_Distribution/`
- **After**: Files saved to `6_Norm_folder/Read_Distribution/`

## Files Modified

### 1. `bin/deseq2_norm_qc.R`
- Removed `Stat_folder` creation and references
- Updated all file paths to remove `ALL` and `ALL_CONDITIONS` subfolders
- Simplified directory structure creation

### 2. `modules/local/deseq2_norm_qc/main.nf`
- Removed `9_Stat_folder/**` from outputs
- Updated file copy commands to reflect new paths
- Updated stub section with simplified folder structure

### 3. `conf/modules.config`
- Removed `9_Stat_folder/**` from publishing patterns
- Updated publishing configuration to match new structure

### 4. `docs/usage.md`
- Updated output structure documentation
- Corrected file paths in examples

## Benefits

1. **Simplified Structure**: Eliminates unnecessary nested folders
2. **Easier Navigation**: Files are more accessible with flatter directory structure
3. **Cleaner Output**: Reduces clutter in results directory
4. **Maintained Functionality**: All analysis outputs are preserved, just better organized

## New Output Structure

```
results/
├── 6_Norm_folder/
│   └── Read_Distribution/
│       ├── Read_Distribution_Raw.pdf
│       └── Read_Distribution_Norm_Filt.pdf
├── 7_Counts_folder/
│   ├── EX_reads_RAW_integer.txt
│   ├── EX_reads_RAW_filt.txt
│   ├── EX_reads_NORM_filt.txt
│   ├── rLog_reads.txt
│   └── Normalisation_Parameters.txt
├── 8_Quality_folder/
│   ├── Heatmap_sampleTosample_distances_vstTransformed.pdf
│   ├── Heatmap_Poisson_Distance.pdf
│   ├── Heatmap_sampleTosample_correlation.pdf
│   ├── PCA_rlogTransformedID.pdf
│   ├── Percent_var_PCA.pdf
│   ├── All_pca.pdf
│   ├── 3D_PCA_Group.html
│   ├── 3D_PCA_Samples.html
│   └── 3D_PCA_ID.html
└── deseq2_results/
    ├── EX_reads_NORM_filt.txt
    └── rLog_reads.txt
```

---

*Changes implemented on: $(date)*
*Pipeline version: LMD_rnaseq*