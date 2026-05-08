#!/bin/bash
#SBATCH --job-name=bwa_bam_qc
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4       
#SBATCH --mem=8G                
#SBATCH --time=04:00:00         
#SBATCH --array=0-3
#SBATCH --output=bam_qc_%A_%a.out
#SBATCH --error=bam_qc_%A_%a.err

source ~/.bashrc
conda activate germline_sv_env  

PATIENTS=("PAPKHU" "PASJJH" "PASNSF" "PATNBV")
PATIENT=${PATIENTS[$SLURM_ARRAY_TASK_ID]}
BAM="${PATIENT}_Germline_T2T.bam"
OUT="${PATIENT}_flagstat.txt"

echo "Running samtools QA on $BAM..."
samtools flagstat -@ 4 $BAM > $OUT
echo "Finished QA for $PATIENT."
