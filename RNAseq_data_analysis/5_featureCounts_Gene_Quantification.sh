#!/bin/bash

#SBATCH --job-name=fCounts 
#SBATCH --partition=long
#SBATCH --nodes=1  
#SBATCH --output=logs/fCounts_quantification.out
#SBATCH --error=logs/fCounts_quantification.err

#=========================
# User defined parameters: relevant paths
#=========================
# SPECIFY Root directory in the cluster (usually /projects/cancer)
ROOTDIR="/projects/cancer"

# SPECIFY your project working directory. 
# NOTE: It is assumed that Raw data is within 'raw data' folder
#       and Trimmed data is within 'trimmed_data' folder

WKD=$ROOTDIR'/HSC_IKB'

# Specify project name. It will be used for gene quantification file
PROJECT='RNAseq_HSC_IKB'

# SPECIFY your GTF annotation reference genome (full path)
# REMARK: It should be the same as the one used for STAR index
GTF=$ROOTDIR'/db_files/Genomes/Ensembl/mouse/mm10/release-102/Mus_musculus.GRCm38.102.gtf'

#=========================
# General configuration
#=========================
START=$(date +%s)
# Enable Singularity image to look into the general path (equivalent to -B)
export SINGULARITY_BIND=$ROOTDIR 
# Path to images folder in cluster
IMAGES_PATH=$ROOTDIR"/images"
# Path to databases folder in cluster
DB_PATH=$ROOTDIR"/db_files"

################################################################################
##       Gene reads quantification with featureCounts
################################################################################

# Link to featureCounts tutorial
# http://subread.sourceforge.net/featureCounts.html

###########################
## 1. Other relevant paths
###########################

# Folder where BAM files are located
DATA=$WKD'/STAR_align/BAM'

# Folder where featureCounts output should be stored
OUT=$WKD'/quantification'

#################################################
## 2. Singularity image and Tool Parametrization
#################################################

# Specify image/s name to be used (tool-related)
FCOUNTS='featureCounts_subRead_v2.0.1.sif'  #This image inludes featureCounts from subRead v2.0.1

# Specify any particular tool parameters
# Number of threads                         
T='8' 

# Stranded library
# 0 is for unstranded, 1 for stranded and 2 for reversely stranded
# Typically = 2, if succesfully assigned reads is extremely low, switch to 1
# THIS IS SINGLE-END READS => TYPICALLY UNSTRANDED
ST='0'                                          

# Feature type (-t)
ft='exon'   #exon is by default
#ft='gene'   #gene feature to account for intronic reads

# Attribute type to group features based on GTF annotations (gene id by default)
at='gene_id'

# Other interesting parameters:
# -p only applicable for paired-end: fragments to be counted instead of reads.
# -C for NOT counting chimeric reads (although already discarded in STAR, by default)
# -B only count read pairs that have both ends aligned
# -M if multi-mappers have to be also counted

################################################################################
## 3. Gene quantification with featureCounts
################################################################################

# Listing all BAM files from which you wanna quantify genes
BAM=$(for f in $DATA/*.bam; do printf '%s ' $f; done)

echo ''
echo 'Begin Read Quantification with featureCounts....................... '`date`
echo ''

singularity exec $IMAGES_PATH/$FCOUNTS featureCounts  -T $T -s $ST -t $ft -g $at -a $GTF -C -o $OUT/$PROJECT'_gene_quantification.txt' $BAM

################################################################################
## 4. End Gene Quantification
################################################################################

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "Processing Time: $DIFF seconds"