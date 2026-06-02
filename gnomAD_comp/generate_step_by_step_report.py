import os
import json

# Your Historical MANTA Pipeline Data (Post-Scrubbing, Pre-Merge)
manta_data = {
    "PAPKHU": {"cuteSV": 26979, "Sniffles": 17764, "nanomonSV": 20329, "SAVANA": 21660, "Consensus": 17508},
    "PASJJH": {"cuteSV": 20968, "Sniffles": 13549, "nanomonSV": 16519, "SAVANA": 17971, "Consensus": 13846},
    "PASNSF": {"cuteSV": 17129, "Sniffles": 10380, "nanomonSV": 13295, "SAVANA": 13930, "Consensus": 10900},
    "PATNBV": {"cuteSV": 14804, "Sniffles": 8723,  "nanomonSV": 11083, "SAVANA": 11848, "Consensus": 9293}
}

patients = ["PAPKHU", "PASJJH", "PASNSF", "PATNBV"]
callers = ["cuteSV", "Sniffles", "nanomonSV", "SAVANA"]
results_dir = "Final_Master_Validation"
output_file = "iAMP21_Step_By_Step_Validation.txt"

def count_vcf_variants(filepath):
    """Safely counts the number of variant lines in an unzipped VCF."""
    if not os.path.exists(filepath):
        return "N/A"
    try:
        with open(filepath, 'r') as f:
            return sum(1 for line in f if not line.startswith('#') and line.strip())
    except Exception:
        return "N/A"

with open(output_file, 'w') as out:
    def write_and_print(text):
        print(text)
        out.write(text + "\n")

    write_and_print("=========================================================================================")
    write_and_print("             iAMP21 Pipeline Validation: Step-by-Step Caller Tracking")
    write_and_print("=========================================================================================\n")
    
    for patient in patients:
        write_and_print(f"### PATIENT: {patient} ###")
        write_and_print(f"{'Caller':<15} | {'MANTA Scrubbed (T2T)':<22} | {'gnomAD Scrubbed (hg38)':<22} | {'Difference'}")
        write_and_print("-" * 85)
        
        # 1. Individual Caller Breakdown
        for caller in callers:
            manta_clean = manta_data[patient][caller]
            
            # Count the newly generated gnomAD file for this caller
            gnomad_file = os.path.join(results_dir, f"{patient}_{caller}_gnomAD_somatic.vcf")
            gnomad_clean = count_vcf_variants(gnomad_file)
            
            if isinstance(gnomad_clean, int):
                diff = gnomad_clean - manta_clean
                diff_str = f"{diff:+,}"
            else:
                diff_str = "N/A"
                
            write_and_print(f"{caller:<15} | {manta_clean:<22,} | {gnomad_clean:<22} | {diff_str}")
        
        # 2. Consensus Breakdown
        write_and_print("-" * 85)
        manta_consensus = manta_data[patient]["Consensus"]
        
        # Count the final gnomAD merged consensus
        gnomad_consensus_file = os.path.join(results_dir, f"{patient}_gnomAD_4caller_consensus.vcf")
        gnomad_consensus = count_vcf_variants(gnomad_consensus_file)
        
        if isinstance(gnomad_consensus, int):
            con_diff = gnomad_consensus - manta_consensus
            con_diff_str = f"{con_diff:+,}"
        else:
            con_diff_str = "N/A"
            
        write_and_print(f"{'FINAL MERGE':<15} | {manta_consensus:<22,} | {gnomad_consensus:<22} | {con_diff_str}")
        
        # 3. Final Truvari Bake-Off Results
        json_path = os.path.join(results_dir, f"{patient}_ULTIMATE_COMPARISON", "summary.json")
        if os.path.exists(json_path):
            with open(json_path, 'r') as f:
                data = json.load(f)
                tp = data.get("TP-base", 0)
                precision = data.get("precision", 0) * 100
                write_and_print(f"\n>> OVERLAP: {tp:,} variants survived BOTH pipelines perfectly (Precision: {precision:.1f}%)")
        
        write_and_print("=========================================================================================\n")

print(f"\nStep-by-step report generated and saved to {output_file}")
