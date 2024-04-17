# HSC_dormancy_Ikba_via_retinoic_acid

This repository includes scripts required for the bulk RNAseq and CUT&Tag data analysis included in Thambyrajah et al. (2024). All scripts include comments so they are self-explanatory.

The repository is organized in the following subfolders:

## RNAseq data analysis folder

Scripts required to reproduce the complete RNAseq data analysis, specifically:

- Data preprocessing: to obtain a raw expression matrix from FASTQ files. Scripts from 1 to 4.
- Downstream analysis: to conduct differential expression analysis and functional analysis (GSEA). Scripts 5 and 6.

To conduct data preprocessing, original FASTQ files are required. Please check GEO accession no. GSE188523 [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE188523]. All required scripts were executed using Singularity images (v3.8.3) per required tool.

It is possible to directly conduct the downstream analysis if the corresponding raw expression matrix from GEO is downloaded. 
NOTE: Gene annotation files are required. For this analysis, these were retrieved from Ensembl (release 102, mm10).

REMARK: Script to conduct PreRanked GSEA on selected signatures also includes the code to generated the enrichment plots included in the publication.

## CUT&Tag data analysis folder

CUT&Tag dataset includes two assays: one of IκBα and other for H3K27me3. Data was analyzed in the following manner:

- QC, trimming and alignment steps were the same for both assays. Scripts 1 to 3.
- BigWig files generation differently configured per assay. Scripts 4a and 4b.
- Rest of scripts referred only to H3K27me3 assay where peak calling and TF discovery with STREME were conducted. Scripts 5 to 10.
  
Please check GEO accession no. GSE188524 [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE188524] to access raw data (FASTQ) or processed files (i.e. BigWig or called peaks). All required scripts were executed using Singularity images (v3.8.3) per required tool.

NOTE: Gene annotation files are required for annotating peaks to genes. For this analysis, these were retrieved from Ensembl (release 102, mm10).

## Other (figures related)
