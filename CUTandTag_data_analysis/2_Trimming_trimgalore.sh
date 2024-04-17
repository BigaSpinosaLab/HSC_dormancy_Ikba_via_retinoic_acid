#!/bin/bash

#SBATCH --job-name=trimming 
#SBATCH --partition=long
#SBATCH --nodes=1  
#SBATCH --output=logs/trim.out
#SBATCH --error=logs/trim.err
#SBATCH --array=1-10%2


#=========================
# User defined parameters: relevant paths
#=========================
# SPECIFY Root directory in the cluster (usually /projects/cancer)
ROOTDIR="/projects/cancer"

# SPECIFY your project working directory. 
# NOTE: It is assumed that Raw data is within 'raw data' folder
#       and Trimmed data is within 'trimmed_data' folder

WKD=$ROOTDIR'/IkBa_HSC_AGM_FL_ABigas/Cut_and_Tag_Ikba_LTHSC_Roshana'

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
##       Reads trimming with trim-galore
################################################################################

# Link to trim galore manual
# https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md

###########################
## 1. Other relevant paths
###########################

# Folder where raw data is stored
RAWDATA=$WKD'/raw_data'

# Folder where trimmed data will be stored
TRIMDATA=$WKD'/trimmed_data'

#################################################
## 2. Singularity image and Tool Parametrization
#################################################

# Specify image/s name to be used (tool-related)
TRIMGALORE='trimgalore_v0.6.6.sif'  #This image inludes TRIMGALORE 0.6.6

# Specify any particular tool parameters
# Quality criteria: by default is 20. For ChIPseq it is not mandatory => less restringent criteria                           
MIN_QUAL='15' 

# Stringency parameter: minimum overlap with adapter sequence to be trimmed (3'). Default: 1b
STRINGENCY='3'

################################################################################
## 3. Command file preparation: to execute batch array
################################################################################

cd $RAWDATA

for FILENAME in *_R1.fastq.gz
do
    BASENAME=${FILENAME%_R1.fastq.gz}
    
    # Construct the full execution command
    echo singularity exec $IMAGES_PATH/$TRIMGALORE trim_galore -o $TRIMDATA \
                                                  --stringency $STRINGENCY \
                                                  -q $MIN_QUAL \
                                                  --paired $RAWDATA/$BASENAME'_R1.fastq.gz' $RAWDATA/$BASENAME'_R2.fastq.gz'
                                                  
done > $WKD'/scripts/cmd/trimgalore_samples_Ikba.cmd'

################################################################################
## 4. TrimGalore execution for all available samples
################################################################################

DATE=$(date +%m-%d-%Y--%T)
echo "Starting Trimming in array mode: $DATE"
echo ''

SEEDFILE=$WKD'/scripts/cmd/trimgalore_samples_Ikba.cmd'
SEED=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SEEDFILE)
eval $SEED

################################################################################
## 4. End
################################################################################

END=$(date +%s)
DIFF=$(( $END - $START ))
echo 'Trimmind completed' 
echo "Processing Time: $DIFF seconds"
