#!/bin/bash
#SBATCH --job-name=ultimate_merge
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G               
#SBATCH --time=02:00:00         
#SBATCH --array=0-3
#SBATCH --output=ultimate_%A_%a.out
#SBATCH --error=ultimate_%A_%a.err

source ~/.bashrc
conda activate germline_sv_env  

PATIENTS=("PAPKHU" "PASJJH" "PASNSF" "PATNBV")
PATIENT=${PATIENTS[$SLURM_ARRAY_TASK_ID]}

# --- TIER 1: Unify SAVANA (1-caller rule) ---
SAVANA_SOMATIC="${PATIENT}_savana_output/${PATIENT}.CCS.reads.T2T_Aligned_extracthifi.MAPQ10.classified.somatic.vcf"
SAVANA_CLASS="${PATIENT}_savana_output/${PATIENT}.CCS.reads.T2T_Aligned_extracthifi.MAPQ10.classified.vcf"
SAVANA_BND="${PATIENT}_savana_output/${PATIENT}.CCS.reads.T2T_Aligned_extracthifi.MAPQ10.sv_breakpoints.vcf"

SAVANA_LIST="${PATIENT}_savana_only_list.txt"
echo -e "$SAVANA_SOMATIC\n$SAVANA_CLASS\n$SAVANA_BND" > $SAVANA_LIST

SAVANA_RAW_MASTER="${PATIENT}_SAVANA_RAW_MASTER.vcf"
SURVIVOR merge $SAVANA_LIST 1000 1 1 1 0 50 $SAVANA_RAW_MASTER

# --- TIER 2: Scrub unified SAVANA ---
MANTA_FILTER="${PATIENT}_Germline_Filter.vcf"
SAVANA_CLEAN="${PATIENT}_savana_somatic_only.vcf"
bedtools intersect -v -a $SAVANA_RAW_MASTER -b $MANTA_FILTER -f 0.5 -r > $SAVANA_CLEAN

# --- TIER 3: Final 4-Caller Consensus (2-caller rule) ---
CUTESV="${PATIENT}_cutesv_somatic_only.vcf"
SNIFFLES="${PATIENT}_sniffles_somatic_only.vcf"
NANOMONSV="${PATIENT}_nanomonsv_parsed.nanomonsv.result.vcf"

FINAL_LIST="${PATIENT}_final_4caller_list.txt"
echo -e "$CUTESV\n$SNIFFLES\n$NANOMONSV\n$SAVANA_CLEAN" > $FINAL_LIST

FINAL_VCF="${PATIENT}_ULTIMATE_consensus_50bp.vcf"
SURVIVOR merge $FINAL_LIST 1000 2 1 1 0 50 $FINAL_VCF

echo "Unbiased, ultra-clean consensus generated: $FINAL_VCF"
