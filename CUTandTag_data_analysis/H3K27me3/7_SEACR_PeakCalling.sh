#!/bin/bash

#SBATCH --job-name=SEACR
#SBATCH --partition=long
#SBATCH --cpus-per-task=2 
#SBATCH --mem=12G
#SBATCH --nodes=1  
#SBATCH --output=logs/SEACR.woDups.out
#SBATCH --error=logs/SEACR.woDups.err
#SBATCH --array=1-8%2

#=========================
# User defined parameters: relevant paths and analysis type
#=========================

# SPECIFY Root directory in the cluster (usually /projects/cancer)
ROOTDIR="/projects/cancer"

# SPECIFY your project working directory. 
WKD=$ROOTDIR'/IkBa_HSC_AGM_FL_ABigas/Cut_and_Tag_Ikba'

# SPECIFY the file name where the sample;input is included. Include a Return in 
# the last row file!
SAMPLESHEET=$WKD"/SEACR_peak_calling_woDups/Samples_Input_Cut_and_Tag.txt"

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
##       SEACR peak calling
################################################################################

# SEACR is intended to call peaks and enriched regions from sparse CUT&RUN or 
# chromatin profiling data in which background is dominated by "zeroes" (i.e. 
# regions with no read coverage)

# Link to SEACR peak caller manual
# https://github.com/FredHutch/SEACR

###########################
## 1. Other relevant paths
###########################

# Folder where input Bedgraph files are available
DATA=$WKD'/Bedgraph_woDups'

# Folder where SEACR output results will be stored: 
OUTPEAKS=$WKD'/SEACR_peak_calling_woDups'

#################################################
## 2. Singularity image and Tool Parametrization
#################################################

# Specify image/s name to be used (tool-related)
SEACR='seacr_v1.3.sif'  #This image inludes SEACR 1.3

# Specify the (full path) script name for SEACR (from github)
SEACR_script=$WKD'/scripts/SEACR_1.3.sh'

# Specify any particular tool parameters
MODE="stringent" # Stringent (recommended) or relaxed

################################################################################
## 3. Command file preparation: to execute batch array
################################################################################

while IFS=";" read -r sample input; 
do
  # Sample name
  NAME=${sample%.fragments.bedgraph}
  
  # Peak calling with SEACR:
    # Without IgG control => applying a threshold
  echo "singularity exec $IMAGES_PATH/$SEACR bash $SEACR_script $DATA/$sample 0.005 norm $MODE $OUTPEAKS/$NAME.woDups.top0.5"
  
done < $SAMPLESHEET > $WKD'/scripts/cmd/SEACR_peak_calling_samples_woDups.cmd'


################################################################################
## 4. SEACR Peak calling
################################################################################

echo "-----------------------------------------"
echo "Starting SEACR Peak Calling"
echo "-----------------------------------------"

DATE=$(date +%m-%d-%Y--%T)
echo "  Samples peak calling in array mode: $DATE"
echo " "

SEEDFILE=$WKD'/scripts/cmd/SEACR_peak_calling_samples_woDups.cmd'
SEED=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SEEDFILE)
eval $SEED

DATE=$(date +%m-%d-%Y--%T)
echo "  All samples peak-called: $DATE"


################################################################################
## 5. End
################################################################################

END=$(date +%s)
DIFF=$(( $END - $START ))
echo 'SEACR peak calling completed' 
echo "Processing Time: $DIFF seconds"


