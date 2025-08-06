import numpy as np
import mdtraj as md
import pandas as pd
import config
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from math import sqrt

def reindex_residues(topology):
    #helps re-index the residue from 0 for each chain if use chain_id selection mode.
    for chain in topology.chains:
        for i, residue in enumerate(chain.residues):
            residue.index = i
    return topology

def gen_res_selection(resid, element, use_chain_ids=False, chain_ids=None, chainid_dict=None):
    if chainid_dict is None:
        chainid_dict = {'A': '0',
                        'B': '1', 
                        'C': '2', 
                        'D': '3', 
                        'E': '4', 
                        'F': '5', 
                        'G': '6', 
                        'H': '7', 
                        'I': '8', 
                        'J': '9', 
                        'K': '10'}

    res_range_list = []  # list to store residue range.
    for start_end_list in resid:
        res_range_list.append('(resid {0} to {1})'.format(start_end_list[0], start_end_list[1]))

    if use_chain_ids and chain_ids:
        chain_selection = ' or '.join('chainid ' + str(chainid_dict[chain_id]) for chain_id in chain_ids)
        return "({0}) and (name {1} and not resname PRO and {2})".format(chain_selection, element, ' or '.join(res_range_list))
    else:
        return "name {0} and not resname PRO and {1}".format(element, ' or '.join(res_range_list))


def calc_NHVecs(traj_file, top_file, resid, start_snap=0, end_snap=None):
    """
    Uses mdtraj to load the trajectory and get the atomic indices and coordinates to calculate the correlation functions.
    For each, trajectory load the trajectory using mdtraj, get the atomic index for the the N-H atoms and calculate the vector between the two.
    Append the vector to the NHVecs list for all the trajectories. 
    NHVecs should return a list of shape: (# Trajectories, # Snapshots, # Residues w/N-H Vectors, 3)
    resid is a list of list. providing the residue start/end id.
    MIT License

    Copyright (c) 2019 Alan Hicks

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
    
    
    """
    if config.traj_segment is None:
        start_frame, end_frame = 0, None
    else:
        start_frame = config.traj_segment[0] if config.traj_segment[0] is not None else 0
        end_frame = config.traj_segment[1] if len(config.traj_segment) > 1 and config.traj_segment[1] is not None else None

    stride = getattr(config, "traj_stride", 1)

    print(f"Loading trajectory frames from {start_frame} to {end_frame if end_frame else 'end'} with stride={stride}...")
    
    traj_full = md.load(traj_file, top=top_file)
    traj = traj_full[start_frame:end_frame][::stride]

    print(f"Trajectory loaded. Total frames after slicing: {traj.n_frames}")
    
        
    
    #print("Loading trajectory...")
    #traj = md.load(traj_file, top=top_file)

    #print("Trajectory loaded.")
  
    
    if config.use_align:
    
        # Extract the first frame (unaligned initial structure)
        first_frame = traj[0]

        # Compute DSSP for the first frame
        ss = md.compute_dssp(first_frame)  # Output shape: (1, n_residues)
    
        # Select Cα atoms based on DSSP of the first frame
        align_indices = []
        for i, residue in enumerate(traj.topology.residues):
            if i >= len(ss[0]):
                break
            if ss[0][i] in ['H', 'G', 'I', 'E', 'B']:
                align_indices.extend([atom.index for atom in residue.atoms if atom.name == 'CA'])
    
        if not align_indices:
            raise ValueError("No suitable Cα atoms found for alignment based on secondary structure.")
    
        # Align the trajectory
        print("Aligning trajectory...")
        traj.superpose(traj, frame=0, atom_indices=align_indices)
        print("Alignment completed.")
        top = traj.topology
    else:
        top = traj.topology
        
    
    if config.use_chain_ids:
        reidx_top = reindex_residues(top)
        Nit = reidx_top.select(gen_res_selection(resid=config.resid, element='N', use_chain_ids=config.use_chain_ids, chain_ids=config.chain_ids))
        Hyd = reidx_top.select(gen_res_selection(resid=config.resid, element='H', use_chain_ids=config.use_chain_ids, chain_ids=config.chain_ids))
    
    else:
        Nit = top.select(gen_res_selection(resid=config.resid, element='N', use_chain_ids=config.use_chain_ids, chain_ids=config.chain_ids))
        Hyd = top.select(gen_res_selection(resid=config.resid, element='H', use_chain_ids=config.use_chain_ids, chain_ids=config.chain_ids))
    
    NH_Pair = [[i,j] for i,j in zip(Nit,Hyd)]
    NH_Pair_Name = [[top.atom(i),top.atom(j)] for i,j in NH_Pair]
    #NH_Res = ["{}-{}{}".format(str(i).split('-')[0],str(i).split('-')[1], str(j).split('-')[1]) for i,j in NH_Pair_Name]
    NH_Res = ["{}".format(str(i).split('-')[0]) for i,j in NH_Pair_Name]

    ##Generate the N-H vectors in Laboratory Frame
    NHVecs_tmp = np.take(traj.xyz, Hyd, axis=1) - np.take(traj.xyz, Nit, axis=1)
    sh = list(NHVecs_tmp.shape)
    sh[2] = 1
    NHVecs_tmp = NHVecs_tmp / np.linalg.norm(NHVecs_tmp, axis=2).reshape(sh)
    
    return NHVecs_tmp[start_snap:end_snap],NH_Res

def split_NHVecs(nhvecs, dt, tau):
    """
    This function will split the trajectory in chunks defined by tau. 
    nhvecs = array of N-H bond vectors,
    dt = timestep of the simulation
    tau = length of chunks
    
    MIT License

    Copyright (c) 2019 Alan Hicks

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
    
    """
    nFiles = len(nhvecs) ## number of trajectories
    nFramesPerChunk = int(tau/dt) ###tau/timestep 
    used_frames = np.zeros(nFiles,dtype=int)
    remainingFrames = np.zeros(nFiles,dtype=int)
    for i in range(nFiles):
        nFrames = nhvecs[i].shape[0]
        used_frames[i] = int(nFrames/nFramesPerChunk)*nFramesPerChunk
        remainingFrames[i] = nFrames % nFramesPerChunk
    
    nFramesTot=int(used_frames.sum())
    out = np.zeros((nFramesTot,nhvecs[0].shape[1],nhvecs[0].shape[2]), dtype=nhvecs[0].dtype)
    start = 0
    for i in range(nFiles):
        end = int(start+used_frames[i])
        endv = int(used_frames[i])
        out[start:end,...] = nhvecs[i][0:endv,...]
        start = end
        
    sh = out.shape
    vecs = out.reshape((int(nFramesTot/nFramesPerChunk), nFramesPerChunk, sh[-2], sh[-1]))
    print('vec shape after split:', vecs.shape)
    
    return vecs

def calc_Ct(nhvecs):
    """
    Calculates the correlation function of the N-H bond vectors found in nhvecs. 
    Direct space calculation. This could be changed to Fourier space calculation for increased speed. 
    
    LICENSE INFO:
    
    MIT License

    Copyright (c) 2017 Po-chia Chen

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

    """
    sh = nhvecs.shape
    nReplicates=sh[0] ; nDeltas=int(sh[1]/2) ; nResidues=sh[2]
    Ct  = np.zeros( (nDeltas, nResidues), dtype=nhvecs.dtype )
    dCt = np.zeros( (nDeltas, nResidues), dtype=nhvecs.dtype )
    
    for delta in range(1,1+nDeltas):
        nVals=sh[1]-delta
        # = = Create < vi.v'i > with dimensions (nRep, nFr, nRes, 3) -> (nRep, nFr, nRes) -> ( nRep, nRes ), then average across replicates with SEM.
        tmp = -0.5 + 1.5 * np.square( np.einsum( 'ijkl,ijkl->ijk', nhvecs[:,:-delta,...] , nhvecs[:,delta:,...] ) )
        tmp  = np.einsum( 'ijk->ik', tmp ) / nVals
        Ct[delta-1]  = np.mean( tmp, axis=0 )
        dCt[delta-1] = np.std( tmp, axis=0 ) / ( np.sqrt(nReplicates) - 1.0 ) #note if replica is 1, this will raise error and no NH_dCt given. (simply because no variance for 1 copy data.)
    
    return Ct, dCt

def _bound_check(func, params):
    """
    
    Checks if the fit returns a sum of the amplitudes greater than 1.
    
    MIT License

    Copyright (c) 2017 Po-chia Chen

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

    """
    if len(params) == 1:
        return False
    elif len(params) %2 == 0 :
        s = sum(params[0::2])
        return (s>1)
    else:
        s = params[0]+sum(params[1::2])
        return (s>1)

def calc_chi(y1, y2, dy=[]):
    """
    Calculates the chi^2 difference between the predicted model and the actual data. 
    
    LICENSE INFO:
    MIT License

    Copyright (c) 2017 Po-chia Chen

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

    """
    #if dy != []:
    if dy is not None and len(dy) > 0:
        return np.sum( (y1-y2)**2.0/dy )/len(y1)
    else:
        return np.sum( (y1-y2)**2.0 )/len(y1)

def calc_r2(y1, y2):
    ss_res = np.sum((y1 - y2)**2)
    ss_tot = np.sum((y1 - np.mean(y1))**2)
    return 1 - (ss_res/ss_tot)

def func_exp_decay(t, *params):

    n_exp = config.n_exp
    #print(params)
    if n_exp == 2:
       A, tau_a, tau_b = params[0], params[1], params[2]
       return ( A * np.exp(-t / tau_a) + (1 - A) * np.exp(-t / tau_b) ) 
    elif n_exp == 3:
       A, tau_a, B, tau_b, tau_c = params[0], params[1], params[2], params[3], params[4]
       return (A * np.exp(-t / tau_a) + B * np.exp(-t / tau_b) + (1 - A - B) * np.exp(-t / tau_c) ) 
    elif n_exp == 4:
         A, tau_a, B, tau_b, C, tau_c, tau_d = params[0], params[1], params[2], params[3], params[4], params[5], params[6]
         return (A * np.exp(-t / tau_a) + B * np.exp(-t / tau_b) + C * np.exp(-t / tau_c) + (1 - A - B - C) * np.exp(-t / tau_d) ) 

def exponential_fitting(x, y, dy=np.empty([]), tau_mem=50.0):
    """
    Fitting using exponential decay functions using scipy curve_fit.
    Return the calculated chi-squared error with calc_chi function.
    Return the calculated r-squared error with calc_r2 function.
    """
    n_exp = config.n_exp 
    if n_exp < 2 or n_exp > 4:
        raise ValueError("Only 2, 3, or 4 exponential components are supported.")
    
    if config.tumbling_time == None:
       tumbling_time = 163.4
    else:
       tumbling_time = config.tumbling_time
    #print(tumbling_time)
    b1_guess = [1.0 / n_exp] * n_exp
    t1_guess = np.logspace(np.log10(0.001), np.log10(round(tumbling_time, 0)), n_exp)
    func = func_exp_decay
    
    if config.use_align:
           #print("align")
           if n_exp == 2:
               guess = (b1_guess[0], t1_guess[0], t1_guess[1])
               bound = ([0.0, 0.0, 0.0], [1.0, tumbling_time, tumbling_time])
    
           elif n_exp == 3:
               guess = (b1_guess[0], t1_guess[0], b1_guess[1], t1_guess[1], t1_guess[2])
               bound = ([0.0, 0.0, 0.0, 0.0, 0.0], [1.0, tumbling_time, 1.0, tumbling_time, tumbling_time])
    
           elif n_exp == 4:
               guess = (b1_guess[0], t1_guess[0], b1_guess[1], t1_guess[1], b1_guess[2], t1_guess[2], t1_guess[3])
               bound = ([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [1.0, tumbling_time, 1.0, tumbling_time, 1.0, tumbling_time, tumbling_time])

    else:
           #print("without_align")
           if n_exp == 2:
               guess = (b1_guess[0], t1_guess[0], t1_guess[1])
               bound = ([0.0, 0.0, 0.0], [1.0, np.inf, np.inf])
    
           elif n_exp == 3:
               guess = (b1_guess[0], t1_guess[0], b1_guess[1], t1_guess[1], t1_guess[2])
               bound = ([0.0, 0.0, 0.0, 0.0, 0.0], [1.0, np.inf, 1.0, np.inf, np.inf])
    
           elif n_exp == 4:
               guess = (b1_guess[0], t1_guess[0], b1_guess[1], t1_guess[1], b1_guess[2], t1_guess[2], t1_guess[3])
               bound = ([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [1.0, np.inf, 1.0, np.inf, 1.0, np.inf, np.inf])
    
    #if dy != []:
    if dy is not None and len(dy) > 0:
        popt, popv = curve_fit(func, x, y, p0=guess, sigma=dy, bounds=bound, method='trf', loss='soft_l1') 
    else:
        popt, popv = curve_fit(func, x, y, p0=guess,  bounds=bound, loss='soft_l1')
    
    #calculate the approximated y using the fitted functions.
    #note this is for error calculations.
    ymodel=[ func(x[i], *popt) for i in range(len(x)) ] 
    
    bExceed=_bound_check(func, popt)
    if bExceed: #disable the sum(weight) = 1, use the global scaling, S^2 ???
        #print("!!! WARNING, curve fitting in do_LSstyle_fit returns a sum>1. !!!")
        #return 9999.99, popt, np.sqrt(np.diag(popv)), ymodel 
        return calc_chi(y, ymodel, dy), calc_r2(y, ymodel), popt, popv, ymodel
    else:
        return calc_chi(y, ymodel, dy), calc_r2(y, ymodel), popt, popv, ymodel


def fitCorrF(CorrDF, dCorrDF, tau_mem, pars_l, threshold=1.0):
    """
        Input Variables:
            CorrDF: Dataframe containing the correlation functions. Columns are the NH-bond vectors, rows are timesteps. 
            dCorrDF: Error in the correlation function at time t
            tau_mem: Cut-Off time to remove noise at the tail of the correlation function 
            pars_l : parameters list. 
            
        Main function to fit the correlation function. 
        Loops over all residues with N-H vectors and calculates the fit, appends the best fit from findbest_Expstyle_fits2.
        Passes the set of lists to fitstoDF to return a data frame of the best fits for each residue. 
        
        Takes the correlation function CorrDF and errors in the correlation function, maximum tau mem to cut correlation
    """
    NH_Res = CorrDF.columns
    chi_list=[]
    r2_list = []
    names_list=[]
    pars_list=[]
    errs_list=[]
    ymodel_list=[]
    covarMat_list = []
    
    #loop on all residues
    for i in CorrDF.columns:
               
        tstop = np.where(CorrDF.index.values==tau_mem)[0][0]
        #print('The tau_max (ns) for truncating the correlation function is:', tstop)
        x = CorrDF.index.values[:tstop]
        y = CorrDF[i].values[:tstop]
        dy = dCorrDF[i].values[:tstop]
        
        ## If there is no error provided then set dy to empty set
        if np.all(np.isnan(dy)):
            dy = []
            
        #print("Now performing a fitting.")    
        chi, r2, pars, covarMat, ymodel = exponential_fitting(x, y, dy, tau_mem)
        
        A_tau_pars = pars
        
        names = ['C_a', 'tau_a', 'C_b', 'tau_b', 'C_g', 'tau_g', 'C_d', 'tau_d']
        errs = np.sqrt(np.diag(covarMat)) # this is the standard covariance error.
        
        #append the calculated params to list
        chi_list.append(chi)
        r2_list.append(r2)
        names_list.append(names)
        pars_list.append(pars)
        errs_list.append(errs)
        ymodel_list.append(ymodel)
        
    FitDF = fitstoDF(NH_Res, chi_list, r2_list, pars_list, errs_list, names_list)
    
    return FitDF

def fitstoDF(resnames, chi_list, r2_list, pars_list, errs_list, names_list):
    ## Set Up columns indices and names for the data frame
    """
    Function that takes the residue names, chi^2, r2, parameters, errors and names of the fits and returns a data frame
    of the parameters.
    """
    mparnames = ['C_a', 'tau_a', 'C_b', 'tau_b', 'C_g', 'tau_g', 'C_d', 'tau_d']
    mtau_names = np.array(mparnames)[1::2]
    mc_names = np.array(mparnames)[::2]
    colnames = np.array(['Resname','NumExp'])
    tau_errnames = np.array([[c,"{}_err".format(c)] for c in mtau_names]).flatten()
    mc_errnames = np.array([[c, "{}_err".format(c)] for c in mc_names]).flatten()
    colnames = np.hstack([colnames,mc_errnames])
    colnames = np.hstack([colnames,tau_errnames])
    colnames = np.hstack([colnames,np.array(['Chi_Fit'])])
    colnames = np.hstack([colnames,np.array(['R2_Fit'])])
    FitDF = pd.DataFrame(index=np.arange(len(pars_list)), columns=colnames).fillna(0.0)
    FitDF['Resname'] = resnames
    FitDF['Chi_Fit'] = chi_list
    FitDF['R2_Fit'] = r2_list
            
    
    for i in range(len(pars_list)):
        npar = len(pars_list[i])
        if (npar%2)==1:
            ccut = npar-2
            tau_f, terr = pars_list[i][1:ccut+1:2], errs_list[i][1:ccut+1:2]
            tau_f = np.hstack([tau_f, pars_list[i][-1]])
            terr = np.hstack([terr, errs_list[i][-1]])
            sort_tau = np.argsort(tau_f)
            coeff, cerr= pars_list[i][0:ccut:2], errs_list[i][0:ccut:2]
            Clast = 1; Clasterr = 0.0;
            for n,m in zip(coeff, cerr):
                Clast -= n
                Clasterr += m
            
            coeff = np.hstack([coeff, np.array(Clast)])
            cerr = np.hstack([cerr, np.array(Clasterr)])
    
            tne = np.array([[c,"{}_err".format(c)] for c in mparnames[1:npar+1:2]]).flatten()
            cne = np.array([[c, "{}_err".format(c)] for c in mparnames[0:npar:2]]).flatten()
                
        else:
            tau_f, terr = pars_list[i][1::2], errs_list[i][1::2] 
            coeff, cerr= pars_list[i][0::2], errs_list[i][0::2]
            sort_tau = np.argsort(tau_f)[::-1]
            tne = np.array([[c,"{}_err".format(c)] for c in names_list[i][1::2]]).flatten()
            cne = np.array([[c, "{}_err".format(c)] for c in names_list[i][0::2]]).flatten()
    
        NumExp=np.array(len(tau_f))
        tau_err = np.array([[t,e] for t,e in zip(tau_f[sort_tau],terr[sort_tau])]).flatten()
        c_err = np.array([[c,e] for c,e in zip(coeff[sort_tau], cerr[sort_tau])]).flatten()
        namesarr = np.hstack([np.array('NumExp'),cne,tne])
        valarr = np.hstack([NumExp,c_err,tau_err])
    
        FitDF.loc[i,namesarr] = valarr
        
    FitDF['AUC_a'] = FitDF.C_a*FitDF.tau_a; FitDF['AUC_b'] = FitDF.C_b*FitDF.tau_b; 
    FitDF['AUC_g'] = FitDF.C_g*FitDF.tau_g; FitDF['AUC_d'] = FitDF.C_d*FitDF.tau_d;
    FitDF['AUC_Total'] = FitDF[['AUC_a','AUC_b','AUC_g','AUC_d']].sum(axis=1)
    
    
    return FitDF            
            
            
def J_direct_transform(om, consts, taus):
    
    """
        Calculation of the spectral density from the parameters of the fit by direct fourier transform
        
        return the sum of estimated spectral density on point: x = om. 
        
        MIT License

    Copyright (c) 2019 Alan Hicks

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
        
    """
    ## Calculation for the direct spectral density 
    if config.tumbling_time == None:
       tumbling_time = 163.4
    else:
       tumbling_time = config.tumbling_time
    
    ndecay=len(consts)-1 ; noms=1;###lnden(om)
    Jmat = np.zeros( (ndecay, noms ) ) #Jmat is in shape [2,1] 2 is the number of fitted decay func
    for i in range(ndecay):
        if config.use_align:
            #print(taus[i])
            #print((taus[i]*163.4*1e-9)/(taus[i]+163.4))
            Jmat[i] = consts[i]*(taus[i]*tumbling_time*1e-9)/(taus[i]+tumbling_time)/(1 + np.power((taus[i]*1e-9*tumbling_time/(taus[i]+tumbling_time))*(om),2.))
        else:
            Jmat[i] = consts[i]*(taus[i]*1e-9)/(1 + np.power((taus[i]*1e-9)*(om),2.))
    
    return (0.89)*Jmat.sum(axis=0)

def calc_NMR_Relax(J, fdd, fcsa, gammaH, gammaN):
    """
        Function to calculate the R1, R2 and NOE from the spectral densities and the physical parameters for the 
        dipole-dipole and csa contributions, fdd and fcsa. 
    """
    R1 = fdd * (J['Diff'] + 3*J['15N'] + 6*J['Sum']) + fcsa * J['15N']
    
    R2 = (1./2. * fdd * (4*J['0'] + J['Diff'] + 3*J['15N'] + 6*J['1H']) 
          + (1./6.) * fcsa*(4*J['0'] + 3*J['15N']) )
    
    NOE = 1 + ((fdd*gammaH)/(gammaN*R1))*(6*J['Sum'] - J['Diff']) #changed constant NOE = 1 + ((fdd*gammaH)/(*gammaN*R1))*(6*J['Sum'] - J['Diff'])
    
    return R1, R2, NOE
