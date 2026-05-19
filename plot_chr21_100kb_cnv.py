import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# The 4 patients we are analyzing
patients = ['PAPKHU', 'PASJJH', 'PASNSF', 'PATNBV']

# In T2T-CHM13, autosomes are NC_060925.1 (chr1) through NC_060946.1 (chr22)
# Chromosome 21 is NC_060945.1
autosomes = [f"NC_0609{i}.1" for i in range(25, 47)]
chr21_id = "NC_060945.1"

for patient in patients:
    file_path = f"{patient}_100kb_depth.regions.bed.gz"
    print(f"Processing {patient} Chromosome 21...")
    
    try:
        # Load the mosdepth output
        df = pd.read_csv(file_path, sep='\t', header=None, names=['chrom', 'start', 'end', 'depth'])
        
        # Filter for all autosomes to calculate the true global baseline
        df_auto = df[df['chrom'].isin(autosomes)].copy()
        global_median = df_auto['depth'].median()
        
        # Now, isolate ONLY Chromosome 21 for plotting
        df_chr21 = df_auto[df_auto['chrom'] == chr21_id].copy()
        
        # Calculate Log2 Ratio
        df_chr21['log2_ratio'] = np.log2((df_chr21['depth'] + 1) / global_median)
        
        # Convert starting base pairs to Megabases (Mb) for a cleaner X-axis
        df_chr21['start_mb'] = df_chr21['start'] / 1_000_000
        
        # Set up the plot
        plt.figure(figsize=(10, 5))
        plt.scatter(df_chr21['start_mb'], df_chr21['log2_ratio'], s=15, color='dodgerblue', alpha=0.8)
        
        # Add the biological reference lines
        plt.axhline(0, color='black', linewidth=1.5, linestyle='-', label='Normal (2 copies)')
        plt.axhline(-1.0, color='red', linewidth=1.5, linestyle='--', label='Deletion (1 copy)')
        plt.axhline(0.58, color='green', linewidth=1.5, linestyle='--', label='Duplication (3 copies)')
        
        # Formatting
        plt.ylim(-2.5, 2.5)
        plt.title(f"Chromosome 21 Copy Number Profile (100kb Bins) - {patient}", fontsize=14)
        plt.xlabel("Position on Chromosome 21 (Mb)", fontsize=12)
        plt.ylabel("Log2 Ratio (Bin/Median)", fontsize=12)
        plt.legend(loc='upper right')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        # Save the figure
        out_file = f"{patient}_chr21_100kb_CNV.png"
        plt.savefig(out_file, dpi=300)
        plt.close()
        print(f" -> Saved {out_file}\n")
        
    except Exception as e:
        print(f" -> Error processing {patient}: {e}\n")
