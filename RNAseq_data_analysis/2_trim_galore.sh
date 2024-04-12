#!/bin/bash
#SBATCH -p fast
#SBATCH -c 2
#SBATCH -N 1 
#SBATCH --mem-per-cpu 5000
#SBATCH -t 05:00:00 
#SBATCH -o trim.out
#SBATCH -e trim.err

ROOT="/projects/cancer"
RAWDATA="/projects/cancer/HSC_IKB/raw_data"
IMAGES="/projects/cancer/images"

cd ${ROOT}/HSC_IKB

mkdir Trim_data


cd ${RAWDATA}

singularity exec -B ${ROOT}:${ROOT} ${IMAGES}/trimgalore.simg trim_galore -q 30 -o ${ROOT}/HSC_IKB/Trim_data INPUT
