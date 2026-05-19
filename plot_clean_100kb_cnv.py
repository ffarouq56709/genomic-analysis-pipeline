import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# Map T2T RefSeq accessions to human-readable chromosome numbers
chrom_map = {
    "NC_060925.1": "1", "NC_060926.1": "2", "NC_060927.1": "3", "NC_060928.1": "4",
    "NC_060929.1": "5", "NC_060930.1": "6", "NC_060931.1": "7", "NC_060932.1": "8",
    "NC_060933.1": "9", "NC_060934.1": "10", "NC_060935.1": "11", "NC_060936.1": "12",
    "NC_060937.1": "13", "NC_060938.1": "14", "NC_060939.1": "15", "NC_060940.1": "16",
    "NC_060941.1": "17", "NC_060942.1": "18", "NC_060943.1": "19", "NC_060944.1": "20",
    "NC_060945.1": "21", "NC_060946.1": "22"
}

patients = ["PAPKHU", "PASJJH", "PASNSF", "PATNBV"]

for patient in patients:
    filename = f"{patient}_100kb_depth.regions.bed.gz"
    if not os.path.exists(filename):
        print(f"Skipping {patient} - could not find {filename}")
        continue
        
    print(f"Plotting {patient}...")
    
    # Load Data
    df = pd.read_csv(filename, sep='\t', names=['chrom', 'start', 'end', 'depth'])
    df = df[df['chrom'].isin(chrom_map.keys())]
    
    # Baseline & Log2 Math
    median_depth = df['depth'].replace(0, np.nan).median()
    df['log2_ratio'] = np.log2((df['depth'] + 1) / median_depth)
    df = df[(df['log2_ratio'] > -2.5) & (df['log2_ratio'] < 2.5)]
    
    # Plotting Setup
    plt.figure(figsize=(16, 6))
    colors = ['#2c3e50', '#7f8c8d']
    x_labels = []
    x_ticks = []
    current_pos = 0
    
    # Plot each chromosome
    for i, (chrom, chrom_name) in enumerate(chrom_map.items()):
        chrom_data = df[df['chrom'] == chrom].copy()
        if chrom_data.empty: continue
        
        chrom_data['plot_pos'] = range(current_pos, current_pos + len(chrom_data))
        
        # Draw the dots (s=3 makes them visible blocks for 100kb)
        plt.scatter(chrom_data['plot_pos'], chrom_data['log2_ratio'], 
                    color=colors[i % 2], s=3, alpha=0.7)
        
        x_ticks.append(current_pos + (len(chrom_data) // 2))
        x_labels.append(chrom_name)
        current_pos += len(chrom_data)

    # Reference Lines & Polish
    plt.axhline(y=0, color='black', linestyle='-', linewidth=1)
    plt.axhline(y=-1, color='#c0392b', linestyle='--', linewidth=1.5, alpha=0.8) # Deletion 
    plt.axhline(y=0.58, color='#2980b9', linestyle='--', linewidth=1.5, alpha=0.8) # Duplication
    
    plt.title(f"Whole Genome Copy Number Profile (100kb Bins) - {patient}", fontsize=16, fontweight='bold')
    plt.ylabel("Log2 Ratio (Bin/Median)", fontsize=14)
    plt.xlabel("Chromosome", fontsize=14)
    plt.xticks(x_ticks, x_labels)
    plt.ylim(-2.0, 2.0)
    
    plt.tight_layout()
    output_img = f"{patient}_100kb_Log2_CNV.png"
    plt.savefig(output_img, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  -> Saved {output_img}")

print("Done! Check your folder for the new PNG images.")
