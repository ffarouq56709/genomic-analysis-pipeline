import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# T2T Accession Map
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
        
    print(f"Extracting Chr 21 column for {patient}...")
    
    # 1. Load Whole Genome Data
    df = pd.read_csv(filename, sep='\t', header=None, names=['chrom', 'start', 'end', 'depth'])
    df = df[df['chrom'].isin(chrom_map.keys())].copy()
    
    # 2. Calculate the exact same global Log2 math
    median_depth = df['depth'].replace(0, np.nan).median()
    df['log2_ratio'] = np.log2((df['depth'] + 1) / median_depth)
    
    # 3. Isolate ONLY Chromosome 21
    chr21_df = df[df['chrom'] == 'NC_060945.1'].copy()
    chr21_df = chr21_df[(chr21_df['log2_ratio'] > -2.5) & (chr21_df['log2_ratio'] < 2.5)]
    
    if chr21_df.empty: 
        continue
        
    # Stack the points sequentially (just like the whole genome plot)
    chr21_df['plot_pos'] = range(len(chr21_df))
        
    # 4. Plotting Setup - Make it a tall, narrow column!
    plt.figure(figsize=(4, 6))
    
    # Use the red color from the whole genome plot to tie them together visually
    plt.scatter(chr21_df['plot_pos'], chr21_df['log2_ratio'], color='#e74c3c', s=15, alpha=0.8, edgecolor='none')
    
    # Reference Lines
    plt.axhline(y=0, color='black', linestyle='-', linewidth=1)
    plt.axhline(y=-1, color='#c0392b', linestyle='--', linewidth=1.5, alpha=0.5)
    plt.axhline(y=0.58, color='#2980b9', linestyle='--', linewidth=1.5, alpha=0.5)
    
    plt.title(f"{patient}", fontsize=18, fontweight='bold', pad=15)
    plt.ylabel("Log2 Ratio", fontsize=14, fontweight='bold')
    
    # Create the single "Chr 21" label perfectly centered at the bottom
    center_pos = len(chr21_df) // 2
    plt.xticks([center_pos], ["Chr 21"], fontsize=14, fontweight='bold')
    
    # Hide the little tick marks on the X-axis for a cleaner look
    plt.tick_params(axis='x', length=0) 
    
    # Set limits to keep it perfectly scaled with the whole-genome plot
    plt.ylim(-2.0, 2.0)
    plt.xlim(-10, len(chr21_df) + 10) # Tiny bit of padding on the sides
    
    # Clean up borders
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['bottom'].set_visible(False) # Removes the bottom line so it looks floating
    
    plt.tight_layout()
    output_img = f"{patient}_Chr21_Column.png"
    plt.savefig(output_img, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  -> Saved {output_img}")
