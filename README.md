# HSC_dormancy_Ikba_via_retinoic_acid

This repository includes scripts required for the bulk RNAseq and CUT&Tag data analysis included in Thambyrajah et al. (2024). All scripts include comments so they are self-explanatory.

The repository is organized in the following subfolders:

## Main figures folder

## RNAseq data analysis folder

Scripts required to reproduce the complete RNAseq data analysis, specifically:

- Data preprocessing: to obtain a raw expression matrix from FASTQ files. Scripts from 1 to 4.
- Downstream analysis: to conduct differential expression analysis and functional analysis (GSEA). Scripts 5 and 6.

To conduct data preprocessing, original FASTQ files are required. Please check GEO accession no. GSE188523. All required scripts were executed using Singularity images (v3.8.3) per each required tool.

It is possible to directly conduct the downstream analysis if the corresponding raw expression matrix from GEO is downloaded. NOTE: Gene annotation files are required. For this analysis, these were retrieved from Ensembl (release 102, mm10).

## CUT&Tag data analysis folder
