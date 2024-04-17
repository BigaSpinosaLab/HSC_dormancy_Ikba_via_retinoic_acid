#!/bin/bash

#SBATCH --job-name=Bg_SEACR
#SBATCH --partition=long
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --nodes=1  
#SBATCH --output=logs/Bedgraph_x_SEACR_woDups_Ikba.out
#SBATCH --error=logs/Bedgraph_x_SEACR_woDups_Ikba.err

# No array execution is considered in this case

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
##       Bedgraph from BAM files for SEACR analysis
################################################################################

# Bedgraph files should reflect density across read pairs rather than individual reads. 
# If starting from BAM files, we recommend converting to paired end BED files using 
# bedtools bamtobed with the -bedpe flag, then selecting the 5' and 3' coordinates of 
# the read pair to generate a new BED3 file, and finally converting that file to a bedgraph 
# using bedtools genomecov.

# Link to SEACR indications
# https://github.com/FredHutch/SEACR

###########################
## 1. Other relevant paths
###########################

# Folder where the input BAM (with no duplicates) to be converted to Bedgraph
DATA=$WKD'/Bowtie_align/BAM_Markdup'

# Folder where Bedgraphs to be stored
#OUT=$WKD'/Bedgraph_wDups_sortedbyName'
OUT=$WKD'/Bedgraph_woDups_Ikba'

#################################################
## 2. Singularity image and Tool Parametrization
#################################################

# Specify image/s name to be used (tool-related)
BEDTOOLS='bedtools_v2.30.0.sif' # This image includes BEDTOOLS 2.30  

# Specify any particular tool parameters
# Location of the genome file used during aligment
CHROM_SIZES=$DB_PATH'/Genomes/Ensembl/mouse/mm10/release-102/mm10_Ensembl_r102_chrom.sizes'

# NOTE: If genome file is missing, it can be created as:
#samtools faidx mygenome.fa
#cut -f 1,2 mygenome.fa.fai > chrom.sizes

################################################################################
## 3. Bedgraphs creation according to SEACR indications: create command and execute
################################################################################

# NOTE: Assuming paired-end data

# Command for samtools execution -> easier to read for below command
BEDTOOLS_exec="singularity exec $IMAGES_PATH/$BEDTOOLS bedtools"

for FILENAME in $DATA/*.NoDups.sorted.bam
do
         NAME=${FILENAME%.NoDups.sorted.bam}
         SAMPLE=$(basename $NAME)

        # Starting from BAM files, convert to paired end BED files using bedtools bamtobed with the -bedpe flag
        $BEDTOOLS_exec bamtobed -bedpe -i $FILENAME > $OUT/$SAMPLE.bed
        
        # Select the 5' and 3' coordinates of the read pair to generate a new BED3 file
        awk '$1==$4 && $6-$2 < 1000 {print $0}' $OUT/$SAMPLE.bed > $OUT/$SAMPLE.clean.bed
        cut -f 1,2,6 $OUT/$SAMPLE.clean.bed | sort -k1,1 -k2,2n -k3,3n > $OUT/$SAMPLE.fragments.bed
        
        # Convert previous file to a bedgraph using bedtools genomecov
        $BEDTOOLS_exec genomecov -bg -i $OUT/$SAMPLE.fragments.bed -g $CHROM_SIZES > $OUT/$SAMPLE.fragments.bedgraph

done

################################################################################
## 4. End
################################################################################

END=$(date +%s)
DIFF=$(( $END - $START ))
echo 'Bedgraph creation completed' 
echo "Processing Time: $DIFF seconds"
