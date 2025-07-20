#this is a config file for NMR relaxation rate calculation
#Written by TW Jan 2023
#This script and the main software is subject to the MIT liscence

### DEFINE PARAMETERS ###

""" 
1. Number of split to truncate the MD trajectory. The averaged N-H vector and its standard deviation will be used to assist scipy.curve_fit to gain a better fitting for the correlation function.
2. Tau_max in nanosecond. This is the threshold value to truncate the correlation function. Normally select 1.8ns or 6.5ns based on result from Musselman et al. 2010
3. The magnetic field strength. This is corresponds to the experiment result. Field Strength = 18.8 for 800 MHz [Teslas], 19.9 for 850 MHz
4. residue start and end index in list format. [[start_resid, end_resid], [start_resid_2, end_resid_2]]. Calculations will be based on residues in this range (including the start and end residues), excluding PRO residues since they don't have N-H on the backbone.
note this is 0-indexed resid. usually in pdb files residues are 1-indexed.
note for nucleosome we have two copies of tails, input like: [[0, 37], [136,172]]
5. Working Directory, where the trajectories and topology files are loaded # change as need. default is '../data/'
6. Output Directory # change as need. default is './'
7. Timestep per frame (ps), note this is the actual time interval between frames in the input trajectory not the simulation step.
8. Number of exponentials; note these are used in the multi-exponential fitting of the autocorrelation function.
9. Align, note whether to remove global tumbling of the nucleosome core particle (NCP) by aligning frames to the histone fold domains.
10. Tumbling time constant in nanoseconds. If the system is a nucleosome, the default is 163.4 ns based on Rabdano et al., 2021
"""

n_split = 20 
tau_max = 1.8
B0 = 14.1 
resid = [[0, 10]]#[[0, 135], [0, 135]]
use_chain_ids = False #True or False
chain_ids = ['C','G']
wd = '../data/'  
od = '../data/'
dt = 100
n_exp = 3

use_align = True
tumbling_time = None 

prefix_list = ['vanilla_run1']#for batch mode.

