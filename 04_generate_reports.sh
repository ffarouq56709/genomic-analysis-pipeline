#!/bin/bash

REPORT_SV="iAMP21_Comprehensive_Cohort_Report.txt"
REPORT_BAM="iAMP21_BWA_Alignment_Quality_Report.txt"

# 1. GENERATE SV REPORT
echo "=================================================================================" > $REPORT_SV
echo "                 iAMP21 Cohort: Comprehensive SV Pipeline Report                 " >> $REPORT_SV
echo "=================================================================================" >> $REPORT_SV
echo -e "\n--- 1. Germline Filtering (MANTA & bedtools) ---\n" >> $REPORT_SV

for p in PAPKHU PASJJH PASNSF PATNBV; do 
    manta=$(grep -v "^#" ${p}_Germline_Filter.vcf 2>/dev/null | wc -l || echo "0")
    cute_raw=$(grep -v "^#" ${p}_cutesv.vcf 2>/dev/null | wc -l || echo "0")
    cute_clean=$(grep -v "^#" ${p}_cutesv_somatic_only.vcf 2>/dev/null | wc -l || echo "0")
    cute_scrubbed=$((cute_raw - cute_clean))
    
    snif_raw=$(grep -v "^#" ${p}_SVs.vcf 2>/dev/null | wc -l || echo "0")
    snif_clean=$(grep -v "^#" ${p}_sniffles_somatic_only.vcf 2>/dev/null | wc -l || echo "0")
    snif_scrubbed=$((snif_raw - snif_clean))
    
    sava_raw=$(grep -v "^#" ${p}_SAVANA_RAW_MASTER.vcf 2>/dev/null | wc -l || echo "0")
    sava_clean=$(grep -v "^#" ${p}_savana_somatic_only.vcf 2>/dev/null | wc -l || echo "0")
    sava_scrubbed=$((sava_raw - sava_clean))

    echo "Patient: $p | MANTA Germline SVs Found: $manta" >> $REPORT_SV
    printf "  > %-10s : %-7s raw -> %-7s clean (%-7s scrubbed)\n" "cuteSV" "$cute_raw" "$cute_clean" "-$cute_scrubbed" >> $REPORT_SV
    printf "  > %-10s : %-7s raw -> %-7s clean (%-7s scrubbed)\n" "Sniffles" "$snif_raw" "$snif_clean" "-$snif_scrubbed" >> $REPORT_SV
    printf "  > %-10s : %-7s raw -> %-7s clean (%-7s scrubbed)\n" "SAVANA" "$sava_raw" "$sava_clean" "-$sava_scrubbed" >> $REPORT_SV
    echo "---------------------------------------------------------------------------------" >> $REPORT_SV
done

echo -e "\n--- 2. Final Cohort Consensus Counts (Post-Merge) ---\n" >> $REPORT_SV
printf "%-8s | %-20s | %-20s\n" "Patient" "Final Somatic SVs" "Status" >> $REPORT_SV
echo "---------------------------------------------------------------------------------" >> $REPORT_SV
for p in PAPKHU PASJJH PASNSF PATNBV; do 
    count=$(grep -v "^#" ${p}_ULTIMATE_consensus_50bp.vcf 2>/dev/null | wc -l || echo "Missing")
    printf "%-8s | %-20s | %-20s\n" "$p" "$count" "Validated" >> $REPORT_SV
done

# 2. GENERATE BAM REPORT
echo "===========================================================================" > $REPORT_BAM
echo "              iAMP21 Cohort: BWA-MEM Normal BAM Quality                    " >> $REPORT_BAM
echo "===========================================================================" >> $REPORT_BAM
printf "\n%-10s | %-15s | %-12s | %-15s\n" "Patient" "Total Reads" "Mapping %" "Properly Paired" >> $REPORT_BAM
echo "---------------------------------------------------------------------------" >> $REPORT_BAM

for p in PAPKHU PASJJH PASNSF PATNBV; do 
    FILE="${p}_flagstat.txt"
    if [ -f "$FILE" ]; then
        total=$(grep "in total" $FILE | head -n 1 | awk '{print $1}')
        map_pct=$(grep "mapped (" $FILE | head -n 1 | awk -F'[(%]' '{print $2}')
        pair_pct=$(grep "properly paired (" $FILE | head -n 1 | awk -F'[(%]' '{print $2}')
        printf "%-10s | %-15s | %-12s | %-15s\n" "$p" "$total" "${map_pct}%" "${pair_pct}%" >> $REPORT_BAM
    else
        printf "%-10s | %-15s | %-12s | %-15s\n" "$p" "Missing" "Missing" "Missing" >> $REPORT_BAM
    fi
done

echo "Both reports successfully generated."
