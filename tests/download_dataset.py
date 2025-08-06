import os
import requests

# Download the dataset from Zenodo
def download_dataset(url, filename):
    response = requests.get(url)
    response.raise_for_status()
    with open(filename, 'wb') as f:
        f.write(response.content)
    print(f"Downloaded {filename} successfully.")
"""
# Define the dataset URL and filenames

"""
# Define the dataset URL and filenames
topology_url = "https://zenodo.org/records/16753027/files/H3.pdb"
trajectory_url = "https://zenodo.org/records/16753027/files/H3_1.xtc"
topology_filename = "./data/H3.pdb"
trajectory_filename = "./data/H3_1.xtc"

# Download the dataset files
download_dataset(topology_url, topology_filename)
download_dataset(trajectory_url, trajectory_filename)


topology_url = "https://zenodo.org/records/16753027/files/WT_rw_run1.pdb"
trajectory_url = "https://zenodo.org/records/16753027/files/WT_rw_run1_2000ns_40ps.xtc"
topology_filename = "./data/WT_rw_run1.pdb"
trajectory_filename = "./data/WT_rw_run1_2000ns_40ps.xtc"

# Download the dataset files
download_dataset(topology_url, topology_filename)
download_dataset(trajectory_url, trajectory_filename)
