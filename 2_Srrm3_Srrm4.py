# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.4
# ---

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

os.chdir("/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_SRRM3_Tasic")

def analyze_srrm_expression(metadata, srrm_counts):
    # Reshape srrm_counts from wide to long format
    srrm_counts_long = srrm_counts.melt(
        id_vars=['gene_name'],
        var_name='sample_name',
        value_name='TPM'
    )
    
    # Merge expression data with metadata
    merged_data = pd.merge(srrm_counts_long, metadata, 
                          on='sample_name')
    
    # Create expression analysis plots
    plt.figure(figsize=(15, 6))
    
    # Plot expression by cell type
    sns.boxplot(data=merged_data, 
                x='cell_class', 
                y='TPM',
                hue='gene_name')
    plt.xticks(rotation=90)
    plt.title('Srrm3 and Srrm4 Expression by Cell Class')
    plt.tight_layout()
    plt.savefig('processed_data/srrm_expression.png')
    
    return merged_data

def analyze_splicing_patterns(merged_data):
    # Check if PSI column exists
    if 'PSI' not in merged_data.columns:
        print("Warning: PSI column not found in the data. Skipping splicing analysis.")
        return None
        
    # Group by gene and cell class to analyze splicing patterns
    splicing_patterns = merged_data.groupby(['gene_name', 'cell_class'])['PSI'].mean().reset_index()
    
    # Create splicing analysis plots
    plt.figure(figsize=(15, 6))
    sns.barplot(data=splicing_patterns,
                x='cell_class',
                y='PSI',
                hue='gene_name')
    plt.xticks(rotation=90)
    plt.title('Srrm3 and Srrm4 Splicing Patterns by Cell Class')
    plt.tight_layout()
    plt.savefig('processed_data/srrm_splicing.png')
    
    return splicing_patterns

# Load metadata
metadata = pd.read_csv("raw_data/GSE115746_complete_metadata_28706-cells.csv.gz")

# Load processed data
srrm_counts = pd.read_csv("processed_data/srrm_counts.csv")

# Run the analysis
merged_data = analyze_srrm_expression(metadata, srrm_counts)

# Only run splicing analysis if we have PSI data
if 'PSI' in merged_data.columns:
    splicing_patterns = analyze_splicing_patterns(merged_data)
else:
    print("Note: PSI data not available. Skipping splicing analysis.")

# Print summary statistics
print("\nSummary Statistics:")
print(merged_data.groupby('gene_name')['TPM'].describe())
