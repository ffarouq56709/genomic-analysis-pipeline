#!/bin/bash
#SBATCH --job-name=iAMP21_Master_Pipeline
#SBATCH --output=iAMP21_Master_Pipeline_%j.out
#SBATCH --error=iAMP21_Master_Pipeline_%j.err
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

echo "==================================================="
echo "  Starting iAMP21 Master Validation Pipeline "
echo "==================================================="

# Load your specific conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate sv_val_env

# 1. Setup the workspace
GNOMAD_DB="gnomad_common_SVs.vcf.gz"
CHAIN="chm13v2-hg38.over.chain.gz"
REF="hg38.fa"

# Clean up any leftover temp files from aborted runs before starting
rm -rf temp_*
mkdir -p Final_Master_Validation

# 2. Master Loop for the 4 Gold-Standard Patients
for patient in PAPKHU PASJJH PASNSF PATNBV; do
    echo "==================================================="
    echo "       STARTING MASTER PIPELINE FOR $patient       "
    echo "==================================================="

    declare -A callers=(
        ["Sniffles"]="${patient}_SVs.vcf"
        ["cuteSV"]="${patient}_cutesv.vcf"
        ["nanomonSV"]="${patient}_nanomonsv_parsed.nanomonsv.result.vcf"
        ["SAVANA"]="${patient}_SAVANA_RAW_MASTER.vcf"
    )

    > ${patient}_gnomAD_caller_list.txt

    # -------------------------------------------------------------
    # PHASE 1: Scrub each raw caller through the Virtual Normal
    # -------------------------------------------------------------
    for caller_name in "${!callers[@]}"; do
        raw_vcf="${callers[$caller_name]}"
        
        echo "--> Processing $caller_name ($raw_vcf)..."

        # Sanitize Raw Caller
        sed -e 's/CIPOS=[^;]*//g' -e 's/CIEND=[^;]*//g' -e 's/;;*/;/g' -e 's/;\t/\t/g' -e 's/\t;/\t/g' $raw_vcf > temp_${caller_name}_clean.vcf

        # Rename and Translate
        bcftools annotate --rename-chrs rename_chrs.txt temp_${caller_name}_clean.vcf -o temp_${caller_name}_renamed.vcf
        CrossMap vcf $CHAIN temp_${caller_name}_renamed.vcf $REF temp_${caller_name}_hg38.vcf
        
        # Sort and compress
        bcftools sort temp_${caller_name}_hg38.vcf -Oz -o temp_${caller_name}_hg38.sorted.vcf.gz
        tabix -p vcf temp_${caller_name}_hg38.sorted.vcf.gz

        # Clear Truvari output directory and run Truvari
        rm -rf temp_gnomad_${caller_name}_eval
        truvari bench -b $GNOMAD_DB -c temp_${caller_name}_hg38.sorted.vcf.gz -o temp_gnomad_${caller_name}_eval --passonly

        # Extract pure somatic variants and unzip them
        FINAL_CALLER_VCF="Final_Master_Validation/${patient}_${caller_name}_gnomAD_somatic.vcf"
        gunzip -c temp_gnomad_${caller_name}_eval/fp.vcf.gz > $FINAL_CALLER_VCF
        
        # Add the scrubbed file to the SURVIVOR list
        echo "$FINAL_CALLER_VCF" >> ${patient}_gnomAD_caller_list.txt
        
        # Clean up
        rm -rf temp_${caller_name}* temp_gnomad_${caller_name}_eval
    done

    # -------------------------------------------------------------
    # PHASE 2: SURVIVOR Merge (gnomAD Consensus)
    # -------------------------------------------------------------
    echo "--> Running SURVIVOR Merge for $patient..."
    # Generate the raw merge
    SURVIVOR merge ${patient}_gnomAD_caller_list.txt 1000 2 1 1 0 50 temp_gnomAD_4caller_consensus_raw.vcf

    # THE FIX: Sanitize SURVIVOR's output before continuing
    sed -e 's/CIPOS=[^;]*//g' -e 's/CIEND=[^;]*//g' -e 's/;;*/;/g' -e 's/;\t/\t/g' -e 's/\t;/\t/g' temp_gnomAD_4caller_consensus_raw.vcf > Final_Master_Validation/${patient}_gnomAD_4caller_consensus.vcf

    # Sort and compress the sanitized consensus
    bcftools sort Final_Master_Validation/${patient}_gnomAD_4caller_consensus.vcf -Oz -o Final_Master_Validation/${patient}_gnomAD_4caller_consensus.sorted.vcf.gz
    tabix -p vcf Final_Master_Validation/${patient}_gnomAD_4caller_consensus.sorted.vcf.gz

    # -------------------------------------------------------------
    # PHASE 3: The Ultimate Bake-Off vs MANTA Consensus
    # -------------------------------------------------------------
    echo "--> Preparing historical MANTA consensus for comparison..."
    HISTORICAL_VCF="${patient}_4caller_consensus_50bp.vcf"
    
    # Sanitize the historical SURVIVOR consensus as well
    sed -e 's/CIPOS=[^;]*//g' -e 's/CIEND=[^;]*//g' -e 's/;;*/;/g' -e 's/;\t/\t/g' -e 's/\t;/\t/g' $HISTORICAL_VCF > temp_manta_clean.vcf
    
    # Rename and LiftOver
    bcftools annotate --rename-chrs rename_chrs.txt temp_manta_clean.vcf -o temp_manta_renamed.vcf
    CrossMap vcf $CHAIN temp_manta_renamed.vcf $REF temp_manta_hg38.vcf
    
    bcftools sort temp_manta_hg38.vcf -Oz -o Final_Master_Validation/${patient}_MANTA_consensus_hg38.sorted.vcf.gz
    tabix -p vcf Final_Master_Validation/${patient}_MANTA_consensus_hg38.sorted.vcf.gz
    rm -f temp_manta_*

    echo "--> Running Final Truvari Comparison..."
    rm -rf Final_Master_Validation/${patient}_ULTIMATE_COMPARISON
    
    truvari bench -b Final_Master_Validation/${patient}_MANTA_consensus_hg38.sorted.vcf.gz \
                  -c Final_Master_Validation/${patient}_gnomAD_4caller_consensus.sorted.vcf.gz \
                  -o Final_Master_Validation/${patient}_ULTIMATE_COMPARISON --passonly

    echo "Finished $patient!"
    echo "==================================================="
done
