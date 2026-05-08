#!/bin/bash
#SBATCH --job-name=run_savana
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --array=1-4

source ~/.bashrc
conda activate savana-test  

# ==========================================
# 1. SETUP
# ==========================================
REF_FASTA="/home/ffarouq/stingy_delly/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna"

SAMPLES=(PAPKHU PASJJH PASNSF PATNBV)
CURRENT_SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID - 1]}

BAM_FILE="${CURRENT_SAMPLE}.CCS.reads.T2T_Aligned_extracthifi.MAPQ10.bam"
OUT_DIR="${CURRENT_SAMPLE}_savana_output"

echo "================================================="
echo "Running SAVANA (Tumour-Only Mode) for: $CURRENT_SAMPLE"
echo "================================================="

# ==========================================
# 2. RUN SAVANA
# ==========================================
savana to \
    --tumour "$BAM_FILE" \
    --ref "$REF_FASTA" \
    --outdir "$OUT_DIR"

echo "DONE with $CURRENT_SAMPLE!"
