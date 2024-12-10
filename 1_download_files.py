# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.4
# ---

import os
import gzip
import requests
import pandas as pd
from pathlib import Path

os.chdir("/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_SRRM3_Tasic")

# Create directories for data organization
def create_directories():
    Path("raw_data").mkdir(exist_ok=True)
    Path("processed_data").mkdir(exist_ok=True)

# Download function with progress tracking
def download_file(url, filename):
    print(f"Downloading {filename}...")
    response = requests.get(url, stream=True)
    total_size = int(response.headers.get('content-length', 0))
    
    with open(filename, 'wb') as f:
        for chunk in response.iter_content(chunk_size=8192):
            if chunk:
                f.write(chunk)
    print(f"Downloaded {filename}")

# URLs for the GSE115746 dataset
base_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115746/suppl/"
files_to_download = {
    "metadata": "GSE115746_complete_metadata_28706-cells.csv.gz",
    "exon_counts": "GSE115746_cells_exon_counts.csv.gz",
    "intron_counts": "GSE115746_cells_intron_counts.csv.gz",
    "accession_table": "GSE115746_accession_table.csv.gz"
}

# Create directories
create_directories()

# Download files
for key, filename in files_to_download.items():
    url = base_url + filename
    output_path = os.path.join("raw_data", filename)
    if not os.path.exists(output_path):
        download_file(url, output_path)

# Process the data
# Load metadata
metadata = pd.read_csv("raw_data/GSE115746_complete_metadata_28706-cells.csv.gz")

# Load exon counts and set gene names correctly
exon_counts = pd.read_csv("raw_data/GSE115746_cells_exon_counts.csv.gz")

exon_counts.head(3)

exon_counts = exon_counts.rename(columns={'Unnamed: 0': 'gene_name'})

exon_counts.head(3)

# Filter for Srrm3 and Srrm4
srrm_genes = ['Srrm3', 'Srrm4']
srrm_counts = exon_counts[exon_counts['gene_name'].isin(srrm_genes)]

srrm_counts

# Save processed data
srrm_counts.to_csv("processed_data/srrm_counts.csv", index=False)
