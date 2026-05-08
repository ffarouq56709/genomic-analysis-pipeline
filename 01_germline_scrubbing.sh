#!/bin/bash
#SBATCH --job-name=germline_scrub
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G               
#SBATCH --time=01:00:00         
#SBATCH --array=0-3
#SBATCH --output=scrub_%A_%a.out
#SBATCH --error=scrub_%A_%a.err

source ~/.bashrc
conda activate germline_sv_env  

PATIENTS=("PAPKHU" "PASJJH" "PASNSF" "PATNBV")
PATIENT=${PATIENTS[$SLURM_ARRAY_TASK_ID]}
MANTA="${PATIENT}_Germline_Filter.vcf"

echo "Scrubbing germline variants for $PATIENT..."

# 1. Scrub cuteSV
bedtools intersect -v -a ${PATIENT}_cutesv.vcf -b $MANTA -f 0.5 -r > ${PATIENT}_cutesv_somatic_only.vcf

# 2. Scrub Sniffles
bedtools intersect -v -a ${PATIENT}_SVs.vcf -b $MANTA -f 0.5 -r > ${PATIENT}_sniffles_somatic_only.vcf

echo "Germline scrubbing complete for cuteSV and Sniffles!"
