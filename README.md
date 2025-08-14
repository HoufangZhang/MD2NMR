# MD2NMR
Written by TW in Jan 2023, Queen's University Panchenko's Lab.
For any enquiries, please contact anna.panchenko@queensu.ca.
Backup email: t.wei@queensu.ca

Subsequently modified and updated by Houfang Zhang and Yunhui Peng in July 2025, Central China Normal University.
For any enquiries please contact yunhuipeng@ccnu.edu.cn

## Description
MD2NMR is a tool for calculating NMR relaxation observables (R1/R2/NOE/T1/T2/Tau_c) directly from MD trajectories. It was initially developed for nucleosome simulation analysis but can be extended to other proteins or complexes. This software is distributed under the MIT license. Some functions are imported from other repositories, as noted individually in the script comments.

## Dependencies
The required packages are listed in 'requirements.txt'.
To create a new environment using Anaconda (replace myenv as appropriate):
conda create --name myenv --file requirements.txt

Or using pip:
pip install virtualenv #create virtual environment
virtualenv test_env
source test_env/bin/activate
pip install numpy==1.23.4 pandas==1.5.1 scikit-learn==1.1.3 scipy==1.9.3 mdtraj==1.9.9

## Usage

For single-file mode (basic usage):

Before running, check the `config.py` file to ensure the parameters are suitable for your calculation.  
You can modify the `prefix_list` and working/output directories as needed.  
Double-check that the magnetic field strength (`B0` in `config.py`) matches the experimental value.


### Data Download Protocol
To download the test MD trajectory, run the following command to place it under the `./data` directory:

python ./tests/download_data.py

### Testing Protocol
After downloading the test MD trajectory, run:
python ./src/md2nmr.py -t H3.pdb -y H3_1.xtc

The result should be in the output directory. Then, change the `$use_chain_ids$` in `config.py` to `"use_chain_ids = False"` and run the following command to test on the other trajectory:

python3 ./src/md2nmr.py -t WT_rw_run1.pdb -y WT_rw_run1_2000ns_40ps.xtc 


### Test on your own trajectory.
Follow the steps below:
1. Check the magnetic field parameters — In config.py, the magnetic field unit is Tesla.
2. Ensure the trajectory is long enough — For accurate calculations, at least 200 ns is recommended. In config.py, the trajectory length should be longer than n_split * tau_max * 2. Example: For a 1000 ns trajectory, a recommended setting is n_split = 50 and tau_max = 1.8.
3. Check your topology file — Usually a .pdb file. Different MD packages (AMBER/GROMACS/NAMD) may use different output formats. If your PDB file contains multiple chains, set use_chain_ids = True; otherwise keep it as False (default).
4. Set the working directory (wd) and output directory (od) — The program will prompt when the trajectory is loaded.

### Batch Mode Protocol
For batch mode, check the `config.py` file and modify the `prefix_list` variable.  
Batch mode will generate results for all trajectory/topology files in the working directory that match the specified prefix.  
Then, run the following command in the terminal:

python ./src/md2nmr.py --batch=True

## Example Run

You can see a step-by-step run example in the Jupyter Notebook:
[Run_nmr_relaxation.ipynb](./Run_nmr_relaxation.ipynb)
