#!/bin/bash

#SBATCH --job-name=BAM_Markdup 
#SBATCH --partition=long
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --nodes=1  
#SBATCH --output=logs/BAM_Markdup.out
#SBATCH --error=logs/BAM_Markdup.err
#SBATCH --array=1-6%3

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
##       Mark Duplicates in BAM files using Samtools
################################################################################

# ChIP-seq data will contain duplicates that should not be considered for peak-calling
# In this step we are going to mark them so we can also check the %duplication rate

# Link to Samtools mark duplicates
# http://www.htslib.org/doc/samtools-markdup.html

###########################
## 1. Other relevant paths
###########################

# Folder where the BAM data to be marked-duplicate (i.e 'Bowtie_align/BAM')
DATA=$WKD'/Bowtie_align/BAM'

# Folder where BAM files MARKDUP will be finally stored
OUTBAM=$WKD'/Bowtie_align/BAM_Markdup'

# Folder where Dup stats files will be finally stored
OUTSTATS=$WKD'/Bowtie_align/Other_results'

#################################################
## 2. Singularity image and Tool Parametrization
#################################################

# Specify image/s name to be used (tool-related)
SAMTOOLS='samtools_v1.15.sif' # This image includes SAMTOOLS 1.15  

# Specify any particular tool parameters
T=4  # Number of threads for samtools

################################################################################
## 3. Mark Duplicates: create command and execute
## NOTE: Duplicates will be directly removed to be used for SEACR analysis
################################################################################

# Command for samtools execution -> easier to read for below command
SAMTOOLS_exec="singularity exec $IMAGES_PATH/$SAMTOOLS samtools"

for FILENAME in $DATA/*.sorted.unique.bam 
do
         NAME=${FILENAME%.sorted.unique.bam}
         SAMPLE=$(basename $NAME)

         # Construct the samtools mark duplicates execution (full for paired-end data)
        # Operations sequence: sort by name rather than position, add ms and MC tags for markdup to use later, sort, mark duplicates (-r flag to remove)

        # PAIRED-END: remove duplicates directly as recommended by SEACR
        echo "$SAMTOOLS_exec collate -Ou $FILENAME | $SAMTOOLS_exec fixmate -m - - | $SAMTOOLS_exec sort - | $SAMTOOLS_exec markdup -r -s - $OUTBAM/$SAMPLE.NoDups.bam 2> $OUTSTATS/$SAMPLE.NoDups.stats.txt "

done > $WKD'/scripts/cmd/BAM_markdup_samples.cmd'

echo "-----------------------------------------"
echo "Starting Marking Duplicates in BAM Files"
echo "-----------------------------------------"

DATE=$(date +%m-%d-%Y--%T)
echo "  Samples Marking Duplicates in array mode: $DATE"
echo " "

SEEDFILE=$WKD'/scripts/cmd/BAM_markdup_samples.cmd'
SEED=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SEEDFILE)
eval $SEED

DATE=$(date +%m-%d-%Y--%T)
echo "  All samples marked-duplicate: $DATE"


################################################################################
## 5. Sort by name (required for SEACR) and create index for BAM removed-duplicates
################################################################################


echo "-----------------------------------------"
echo "BAM Files: Sort by name"
echo "-----------------------------------------"


# for FILENAME in $OUTBAM/*.NoDups.bam
# do
#     NAME=${FILENAME%.NoDups.bam}
#     SAMPLE=$(basename $NAME)
#     
#     echo "$SAMTOOLS_exec sort -n $FILENAME -o $OUTBAM/$SAMPLE.NoDups.sorted.bam"
#     
# done > $WKD'/scripts/cmd/BAM_sorted.cmd'
# 
# SEEDFILE=$WKD'/scripts/cmd/BAM_sorted.cmd'
# SEED=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SEEDFILE)
# eval $SEED


################################################################################
## 5. End
################################################################################

END=$(date +%s)
DIFF=$(( $END - $START ))
echo 'Samtools Removed Duplicates completed' 
echo "Processing Time: $DIFF seconds"
