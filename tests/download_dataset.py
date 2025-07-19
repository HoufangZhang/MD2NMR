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
topology_url = "https://zenodo.org/record/7806382/files/H2A.pdb"
trajectory_url = "https://zenodo.org/record/7806382/files/H2A_1.xtc"
topology_filename = "../data/H2A.pdb"
trajectory_filename = "../data/H2A_1.xtc"

# Download the dataset files
download_dataset(topology_url, topology_filename)
download_dataset(trajectory_url, trajectory_filename)
"""
# Define the dataset URL and filenames
topology_url = "https://zenodo.org/record/7806382/files/vanilla_run1.pdb"
trajectory_url = "https://zenodo.org/record/7806382/files/vanilla_run1.xtc"
topology_filename = "../data/vanilla_run1.pdb"
trajectory_filename = "../data/vanilla_run1.xtc"

# Download the dataset files
download_dataset(topology_url, topology_filename)
download_dataset(trajectory_url, trajectory_filename)


# Define the dataset URL and filenames
topology_url = "https://zenodo.org/record/7806382/files/H2A.pdb"
trajectory_url = "https://zenodo.org/record/7806382/files/H2A_1.xtc"
topology_filename = "../data/H2A.pdb"
trajectory_filename = "../data/H2A_1.xtc"

# Download the dataset files
download_dataset(topology_url, topology_filename)
download_dataset(trajectory_url, trajectory_filename)
