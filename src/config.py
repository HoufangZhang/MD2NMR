#this is a config file for NMR relaxation rate calculation
#Originally written by TW in January 2023, and later modified and updated by Houfang Zhang and Yunhui Peng in July 2025.
#This script and the main software is subject to the MIT licence

### DEFINE PARAMETERS ###


"""
1. Number of splits to truncate the MD trajectory. The averaged N–H vector and its standard deviation will be used to assist scipy.curve_fit in achieving a better fit for the correlation function.
2. tau_max in nanoseconds. This is the threshold value to truncate the correlation function. Normally select 1.8 ns or 6.5 ns based on results from Musselman et al., 2010.
3. The magnetic field strength. This corresponds to the experimental setup. Field Strength = 18.8 for 800 MHz [Tesla], 19.9 for 850 MHz.
4. Residue start and end indices in list format: [[start_resid, end_resid], [start_resid_2, end_resid_2]]. Calculations will be based on residues in these ranges (inclusive), excluding PRO residues since they lack a backbone N–H.
   Note: this is a 0-indexed residue ID. In most PDB files, residues are 1-indexed.
   For nucleosomes, there are two copies of the tails; input format example: [[0, 37], [136, 172]].
5. Working directory — where the trajectories and topology files are located. Change as needed. Default: '../data/'.
6. Output directory. Change as needed. Default: './'.
7. Timestep per frame (ps). This is the actual time interval between frames in the input trajectory, not the MD integration timestep.
8. Number of exponentials. These are used in the multi-exponential fitting of the autocorrelation function.
9. Align — whether to remove global tumbling of the nucleosome core particle (NCP) by aligning frames to the histone fold domains.
10. Tumbling time constant in nanoseconds. For nucleosomes, the default is 163.4 ns based on Rabdano et al., 2021.
11. traj_segment — start and end snapshot indices to be used from the trajectory, in the format [start_snap, end_snap]. Set to None to use the entire trajectory.
12. traj_stride — the interval between consecutive frames read from the trajectory (e.g., 1 processes every frame, 2 skips 1 frame per read, N skips N-1 frames).
"""
n_split = 90 
tau_max = 6
B0 = 14.1 
resid = [[0, 37]]#[[0, 135], [0, 135]]
use_chain_ids = False #True or False
chain_ids = ['H']
wd = './data/'  
od = './data/'
dt = 40
n_exp = 3
traj_segment = [5000, 50000]  #[start_snap:end_snap] or None
traj_stride = 1
use_align = True
tumbling_time = None 

prefix_list = ['WT_rw_run1']#for batch mode.

