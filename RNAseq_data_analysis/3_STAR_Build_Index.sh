#!/bin/bash

#SBATCH --job-name=index_STAR 
#SBATCH --partition=fast
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --nodes=1  
#SBATCH --output=logs/index_STAR.out
#SBATCH --error=logs/index_STAR.err

#=========================
# User defined parameters: relevant paths
#=========================
# SPECIFY Root directory in the cluster (usually /projects/cancer)
ROOTDIR="/projects/cancer"

# SPECIFY your project working directory. STAR index will be stored there
WKD=$ROOTDIR'/HSC_IKB'

# SPECIFY your FASTA sequence reference genome (full path)
FASTA=$ROOTDIR'/db_files/Genomes/Ensembl/mouse/mm10/release-102/Mus_musculus.GRCm38.dna.primary_assembly.fa'

# SPECIFY your GTF annotation reference genome (full path)
GTF=$ROOTDIR'/db_files/Genomes/Ensembl/mouse/mm10/release-102/Mus_musculus.GRCm38.102.gtf'

# ReadLength-1 for sjdbOverhand  (ReadLength corresponds to Maximum read length available)
LENGTH='49' 

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
##       STAR build index: required for STAR alignment
################################################################################

# STAR index SHALL be re-computed in following cases:
#	 .- different reference genome
#	 .- different STAR version
#	 .- your reads (maximum length) considerably differ (i.e. 50bp vs 150bp). In 
#	    case you are interested in splicing junctions, you SHOULD always re-build
# 		the index with the proper length

# Link to STAR manual
# https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

###########################
## 1. Other relevant paths
###########################

# Folder where to store the STAR index
INDEX=$WKD'/STAR_align/Index'

#################################################
## 2. Singularity image and Tool Parametrization
#################################################

# Specify image/s name to be used (tool-related)
STAR='star_2.7.8a.sif'  #This image inludes START 2.7.8

# Specify any particular tool parameters
# Number of threads                         
T='8' 


################################################################################
## 3. Build Genome Index
################################################################################

DATE=$(date +%m-%d-%Y--%T)
echo "Starting building Genome Index: $DATE"
echo ''

# NOTE: Include this parameter if using gff3 instead of gtf annotations
#--sjdbGTFtagExonParentTranscript Parent

singularity exec $IMAGES_PATH/$STAR STAR --runMode genomeGenerate --genomeDir $INDEX --genomeFastaFiles $FASTA --sjdbOverhang $LENGTH --sjdbGTFfile $GTF --runThreadN $T 

################################################################################
## 4. End
################################################################################

END=$(date +%s)
DIFF=$(( $END - $START ))
echo 'Genome Index Built' 
echo "Processing Time: $DIFF seconds"
