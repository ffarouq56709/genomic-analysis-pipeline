#!/bin/bash
#SBATCH --job-name=run_sniffles
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=0-3

source ~/.bashrc
conda activate sniffles_env

PATIENTS=("PAPKHU" "PASJJH" "PASNSF" "PATNBV")
PATIENT=${PATIENTS[$SLURM_ARRAY_TASK_ID]}
TUMOR_BAM="${PATIENT}.CCS.reads.T2T_Aligned_extracthifi.MAPQ10.bam"

echo "Running Sniffles for $PATIENT..."

sniffles -i $TUMOR_BAM \
         -v ${PATIENT}_SVs.vcf \
         -t 8 \
         --minsupport 3
