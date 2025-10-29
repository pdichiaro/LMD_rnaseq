**LMDseq** is a bioinformatics pipeline designed for RNA sequencing analysis of **low input samples obtained through Laser Microdissection (LMD)**. This pipeline addresses the unique analytical challenges of processing RNA-seq data from microdissected tissue samples, which typically yield very small amounts of RNA and require specialized computational approaches to ensure robust and reproducible results.

The pipeline implements the comprehensive methodology described in:

> **Di Chiaro P, et al.** "Mapping functional to morphological variation reveals the basis of regional extracellular matrix subversion and nerve invasion in pancreatic cancer." *Cancer Cell* 2024.

This work demonstrates the application of LMD-RNAseq for spatially resolved transcriptomic analysis, providing insights into tumor heterogeneity and the tumor microenvironment.


### Analytical Workflow

```
Raw FASTQ Files
    ↓
Quality Control (FastQC)
    ↓
Adapter Trimming (TrimGalore)
    ↓
Pseudo-alignment (Kallisto)
    ↓
Gene Expression Matrix Creation
    ↓
DESeq2 Normalization & QC
    ↓
Comprehensive Reporting (MultiQC)
```



