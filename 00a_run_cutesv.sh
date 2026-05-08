#!/bin/bash
#SBATCH --job-name=run_cutesv
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=0-3

source ~/.bashrc
conda activate germline_sv_env

PATIENTS=("PAPKHU" "PASJJH" "PASNSF" "PATNBV")
PATIENT=${PATIENTS[$SLURM_ARRAY_TASK_ID]}
REF="/home/ffarouq/stingy_delly/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna"
TUMOR_BAM="${PATIENT}.CCS.reads.T2T_Aligned_extracthifi.MAPQ10.bam"

echo "Running cuteSV for $PATIENT..."

mkdir -p ./${PATIENT}_cutesv_tmp

cuteSV $TUMOR_BAM $REF ${PATIENT}_cutesv.vcf ./${PATIENT}_cutesv_tmp \
    -t 8 -s 3 --genotype \
    --max_cluster_bias_INS 1000 \
    --diff_ratio_merging_INS 0.9 \
    --max_cluster_bias_DEL 1000 \
    --diff_ratio_merging_DEL 0.5

# Clean up temp folder
rm -r ./${PATIENT}_cutesv_tmp
