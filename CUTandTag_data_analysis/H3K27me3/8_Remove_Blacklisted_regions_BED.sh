#!/bin/bash

#SBATCH --job-name=Black_Regions
#SBATCH --partition=long
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --nodes=1  
#SBATCH --output=logs/Black_Regions.out
#SBATCH --error=logs/Black_Regions.err
#SBATCH --array=1-6%2

#=========================
# User defined parameters: relevant paths
#=========================
# SPECIFY Root directory in the cluster (usually /projects/cancer)
ROOTDIR="/projects/cancer"

WKD=$ROOTDIR'/IkBa_HSC_AGM_FL_ABigas/Cut_and_Tag_Ikba'

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
##       Remove BlackListed regions: ENCODE Black list
## https://github.com/Boyle-Lab/Blacklist
## These includes: High Signal Regions and Low Mappability region
################################################################################

###########################
## 1. Other relevant paths
###########################

# Folder where the BED data to be filtered is stored 
DATA=$WKD'/SEACR_peak_calling_woDups'

# Folder where to store BED files without blacklisted regions
OUT=$WKD'/SEACR_peak_calling_woDups'

# BED file including blacklisted regions
BL=$DB_PATH'/Blacklisted_Regions/mm10-blacklist.v2_woCHR.bed'

#################################################
## 2. Singularity image and Tool Parametrization
#################################################

# Specify image/s name to be used (tool-related)
SAMTOOLS='samtools_v1.15.sif'   
BEDTOOLS='bedtools_v2.30.0.sif'

# Specify any particular tool parameters
T=4  # Number of threads for samtools

################################################################################
## 3. Remove black listed regions: create command
################################################################################

# Command for samtools execution -> easier to read for below command
SAMTOOLS_exec="singularity exec $IMAGES_PATH/$SAMTOOLS samtools"

for FILENAME in $DATA/*.top0.5.stringent.bed
do
        NAME=${FILENAME%.bed}
        SAMPLE=$(basename $NAME)

        echo "singularity exec $IMAGES_PATH/$BEDTOOLS bedtools intersect -a $FILENAME -b $BL -v  > $OUT/$SAMPLE.filtered_blacklisted.bed"

done > $WKD'/scripts/cmd/Remove_BL_bis.cmd'

################################################################################
## 4. Execute removing black-listed regions
################################################################################

echo "-----------------------------------------"
echo "Removing black-listed regions"
echo "-----------------------------------------"

DATE=$(date +%m-%d-%Y--%T)
echo "  Creating BAM index in array mode: $DATE"
echo " "

SEEDFILE=$WKD'/scripts/cmd/Remove_BL_bis.cmd'
SEED=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SEEDFILE)
eval $SEED

DATE=$(date +%m-%d-%Y--%T)
echo "Blacklisted regions removed: $DATE"

################################################################################
## 5. End
################################################################################

END=$(date +%s)
DIFF=$(( $END - $START ))
echo 'BL removed' 
echo "Processing Time: $DIFF seconds"
