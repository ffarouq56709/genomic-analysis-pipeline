import json
import os

patients = ["PAPKHU", "PASJJH", "PASNSF", "PATNBV"]
results_dir = "Final_Master_Validation"

# Output file
output_file = "iAMP21_Cohort_Validation_Report.txt"

with open(output_file, 'w') as out:
    def write_and_print(text):
        print(text)
        out.write(text + "\n")

    write_and_print("="*105)
    write_and_print(f"{'iAMP21 Cohort: 4-Caller Consensus Validation (MANTA Germline vs. gnomAD Virtual Normal)':^105}")
    write_and_print("="*105)
    write_and_print(f"{'Patient':<10} | {'MANTA (True)':<15} | {'gnomAD (Test)':<15} | {'Exact Overlap':<15} | {'Precision':<10} | {'Recall':<10} | {'F1-Score':<10}")
    write_and_print("-" * 105)

    total_base = 0
    total_comp = 0
    total_tp = 0
    total_f1 = 0
    valid_patients = 0

    for patient in patients:
        json_path = os.path.join(results_dir, f"{patient}_ULTIMATE_COMPARISON", "summary.json")
        
        if os.path.exists(json_path):
            with open(json_path, 'r') as f:
                data = json.load(f)
                
                # Extract True Baseline (MANTA) and Test (gnomAD) counts
                base_cnt = data.get("base cnt", 0)
                comp_cnt = data.get("comp cnt", 0)
                tp = data.get("TP-base", 0)
                
                # Extract Statistics
                precision = data.get("precision", 0)
                recall = data.get("recall", 0)
                f1 = data.get("f1", 0)
                
                # Accumulate for averages
                total_base += base_cnt
                total_comp += comp_cnt
                total_tp += tp
                total_f1 += f1
                valid_patients += 1
                
                # Format to 2 decimal places for percentages
                prec_str = f"{precision*100:.1f}%" if precision else "N/A"
                rec_str = f"{recall*100:.1f}%" if recall else "N/A"
                f1_str = f"{f1:.3f}" if f1 else "N/A"
                
                write_and_print(f"{patient:<10} | {base_cnt:<15} | {comp_cnt:<15} | {tp:<15} | {prec_str:<10} | {rec_str:<10} | {f1_str:<10}")
        else:
            write_and_print(f"{patient:<10} | ERROR: summary.json not found in {json_path}")

    write_and_print("="*105)
    write_and_print("COHORT AVERAGES & TOTALS")
    write_and_print("-" * 105)
    
    if valid_patients > 0:
        avg_f1 = total_f1 / valid_patients
        overall_precision = (total_tp / total_comp) * 100 if total_comp > 0 else 0
        overall_recall = (total_tp / total_base) * 100 if total_base > 0 else 0
        
        write_and_print(f"Total True Somatic SVs Evaluated (MANTA Baseline): {total_base}")
        write_and_print(f"Total Somatic SVs Called by gnomAD Pipeline:       {total_comp}")
        write_and_print(f"Total SVs Verified by BOTH Pipelines (Overlap):    {total_tp}")
        write_and_print("")
        write_and_print(f"Overall Cohort Precision: {overall_precision:.1f}% (Percentage of gnomAD calls that were true)")
        write_and_print(f"Overall Cohort Recall:    {overall_recall:.1f}% (Percentage of MANTA calls successfully recovered)")
        write_and_print(f"Average Cohort F1-Score:  {avg_f1:.3f}")
    write_and_print("="*105)

print(f"\nReport successfully saved to: {output_file}")
