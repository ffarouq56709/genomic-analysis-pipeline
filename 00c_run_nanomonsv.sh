#!/bin/bash
#SBATCH --job-name=run_nanomonsv
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --array=1-4

source ~/.bashrc
conda activate nanomonsv  # Make sure nanomonsv is installed in this env

# ==========================================
# 1. SETUP
# ==========================================
REF_FASTA="/home/ffarouq/stingy_delly/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna"

SAMPLES=(PAPKHU PASJJH PASNSF PATNBV)
CURRENT_SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID - 1]}

BAM_FILE="${CURRENT_SAMPLE}.CCS.reads.T2T_Aligned_extracthifi.MAPQ10.bam"
NORMAL_BAM="${CURRENT_SAMPLE}_Germline_T2T.bam"

PARSE_PREFIX="${CURRENT_SAMPLE}_nanomonsv_parsed"
FINAL_PREFIX="${CURRENT_SAMPLE}_nanomonsv_final"

echo "================================================="
echo "Running nanomonSV for: $CURRENT_SAMPLE"
echo "================================================="

# ==========================================
# 2. STEP 1: PARSE THE TUMOR BAM
# ==========================================
echo "Parsing BAM file for SV signatures..."
nanomonsv parse "$BAM_FILE" "$PARSE_PREFIX"

# ==========================================
# 3. STEP 2: GET SOMATIC SVS
# ==========================================
echo "Extracting somatic SVs against the normal baseline..."
nanomonsv get "$PARSE_PREFIX" "$NORMAL_BAM" "$REF_FASTA" \
    --out_prefix "$PARSE_PREFIX" \
    --processes 8

echo "Finished processing $CURRENT_SAMPLE."
