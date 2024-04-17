#!/bin/bash

#SBATCH --job-name=BigWig 
#SBATCH --partition=long
#SBATCH --cpus-per-task=4 
#SBATCH --mem=12G
#SBATCH --nodes=1  
#SBATCH --output=logs/BigWig_Paper_H3K27me3.out
#SBATCH --error=logs/BigWig_Paper_H3K27me3.err
#SBATCH --array=1-6%3

#=========================
# User defined parameters: relevant paths and analysis type
#=========================

# SPECIFY Root directory in the cluster (usually /projects/cancer)
ROOTDIR="/projects/cancer"

# SPECIFY your project working directory. 
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
##       Create BigWig files with bamCoverage from deepTools
################################################################################

# Take alignment of reads or fragments (BAM format) and generate a coverage track
# in bigWig format.
# BigWig will be individually created (not compared i.e. to an input control)

# Link to deepTools > bamCoverage 
# https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html

###########################
## 1. Other relevant paths
###########################

# Folder where input BAM files are available
DATA=$WKD'/Bowtie_align/BAM'

# Folder where BigWig files will be stored: 
OUTBW=$WKD'/BigWig_H3K27me3'

#################################################
## 2. Singularity image and Tool Parametrization
#################################################

# Specify image/s name to be used (tool-related)
DEEPTOOLS='deepTools_v3.5.1.simg'  #This image inludes deepTools v3.5.1

# Specify any particular tool parameters
# Effective genome size
# NOTE: Effective genome size is directly incorporated in the samplesheet
# since we have to apply a diff value depending on the sample
GSIZE=2308125349 # mm10 50bp read length 
GSIZE=2407883243 # mm10 75bp read length 

# Type of normalization to be used
NORM=RPGC  # Other options: CPM, RPKM, BPM, RPGC or None

# Number of processors
T=4

# BinSize (by default is 50)
BS=5

# Smoothing (should be larger than BinSize)
SMOOTH=60

# Reads extension (bp)
# Not need to specify extend value since they are paired-end
#EXTEND=200

# Specify the SampleSheet
SAMPLESHEET=$WKD'/BigWig_PAPER_H3K27me3/Samples_Input_PAPER_BW.txt'

################################################################################
## 3. Command file preparation: to execute batch array
################################################################################

# Duplicates kept!

while IFS=";" read -r sample gsize;
do
   # Sample name
   NAME=${sample%.sorted.unique.bam}

   # Conf: small bin size, smoothing, no centering (--centerReads), extension (-e)
  echo "singularity exec $IMAGES_PATH/$DEEPTOOLS bamCoverage -b $DATA/$sample --effectiveGenomeSize $gsize --binSize $BS --smoothLength $SMOOTH --extendReads  --normalizeUsing $NORM -of bigwig -p $T -o $OUTBW/$NAME.$NORM.final.bw"

done < $SAMPLESHEET > $WKD'/scripts/cmd/BigWig_H3K27me3.cmd'


################################################################################
## 4. Create BigWig files
################################################################################

# echo "-----------------------------------------"
# echo "Starting BigWig creation"
# echo "-----------------------------------------"
# 
# DATE=$(date +%m-%d-%Y--%T)
# echo "  BigWig creation in array mode: $DATE"
# echo " "
# 
SEEDFILE=$WKD'/scripts/cmd/BigWig_H3K27me3.cmd'
SEED=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SEEDFILE)
eval $SEED
# 
# DATE=$(date +%m-%d-%Y--%T)
# echo "  All bigwig files created: $DATE"


################################################################################
## 5. End
################################################################################

END=$(date +%s)
DIFF=$(( $END - $START ))
echo 'BigWig creation completed' 
echo "Processing Time: $DIFF seconds"


