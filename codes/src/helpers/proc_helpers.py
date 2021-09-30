#!/usr/bin/env python3
## Example Python script calling DMRG.

# load the pyten module
import pyten as ptn
import numpy as np
import sys, csv
from tabulate import tabulate
from ..utils import SimpleNamespace

def params(Lx=6,
           Ly=4,
           Nphi=2.0,
           U=8.0,
           N=2,
           S=None,                 #SU2 DMRG
           Sz=None,                  #U1 DMRG
           pbc=True,
           bond=100,
           ind=0,
           sweep = [38, 2],
           chis = [],
           source="",            #in which correpoding MPS file is stored
           lat=ptn.mp.lat.nil.genBosonLattice(4,2.),  #simple lattice for initialization purposes
           pin=False,
           g=0
          ):
    return SimpleNamespace(Lx=Lx, Ly=Ly, Nphi=Nphi, U=U, N=N, S=S, Sz=Sz, pbc=pbc, bond=bond, ind=ind, sweep=sweep, chis=chis, source=source, lat=lat, pin=pin, g=g)

def mps_load(params):
    '''
    params: SimpleNamespace object containing necessary parameters
    
    Returns: MPS object in which MPS specified above is stored
    
    '''
    Lx=params.Lx; Ly=params.Ly;
    loc_dr = params.source
    mps_fl = mps_nm(params)
    
    mps_obj = ptn.mp.MPS(loc_dr + "/" + mps_fl)  #load the mps object from the specified directory
    
    return mps_obj

def mps_nm(par):
    '''
    returns: frequently used name format for parameter object par in this project, for MPS,txt,dat files etc
    '''
    
    if par.S is not None and par.Sz is None:
        mps_str = "Lx"+str(par.Lx)+"_Ly"+str(par.Ly)+"_Nphi"+str(par.Nphi)+"_U"+str(par.U)+"_N"+str(par.N)+"_S"+str(par.S)+"_PBC"+\
                    str(par.pbc)+str(par.bond)+ "_"+str(par.ind)
        
    elif par.S is None and par.Sz is not None:
        mps_str = "Lx"+str(par.Lx)+"_Ly"+str(par.Ly)+"_Nphi"+str(par.Nphi)+"_U"+str(par.U)+"_N"+str(par.N)+"_Sz"+str(par.Sz)+"_PBC"+\
                    str(par.pbc)+str(par.bond)+ "_"+str(par.ind)
        
    elif par.S is None and par.Sz is None:
        raise ValueError("None of Spin numbers is initiated, cannot decide whether SU2 or U1")
    else:
        raise ValueError("Both Spin numbers is initiated, cannot decide whether SU2 or U1")
    if par.pin:
        mps_str += "_g"+str(par.g) 
    
    return mps_str + ".mps"


def name_to_S(name):
    '''
    name:string representing the name of mps file of which spin value to be returned
    returns: the total spin number of the mps file
    '''
    idx = name.find("_S")
    rest = name[idx+2:]
    idx2 = rest.find("_")
    S = float(rest[:idx2])

    return S

def name_to_Nphi(name):
    '''
    name:string representing the name of mps file of which spin value to be returned
    returns: the total spin number of the mps file
    '''
    idx = name.find("_Nphi")
    rest = name[idx+5:]
    idx2 = rest.find("_")
    Nphi = float(rest[:idx2])

    return Nphi

def name_to_N(name):
    '''
    name:string representing the name of mps file of which spin value to be returned
    returns: the total spin number of the mps file
    '''
    idx = name.rfind("_N")
    rest = name[idx+2:]
    idx2 = rest.find("_")
    N = float(rest[:idx2])

    return N


def n_arr_save(mps_obj, tar_loc, params):
    '''
    mps_obj: MPS object
    tar_loc: location in which the resultant particle density array to be stored, target location
    params: SimpleNamespace object containing necessary parameters
    
    Stores numpy array for particle density in the target location
    Returns: the name of file in which the array is stored, convenient for extracting the file later
    '''
    
    #produce the correct lattice for mps
    lat=params.lat
    Lx=params.Lx; Ly=params.Ly
    
    #particle density distribution as a function of x-coord
    n_arr = np.zeros((Lx, Ly))               #use the same Lx and Ly used above

    for i in range(Lx):
        for j in range(Ly):
            ind = Ly*i + j

            n_arr[i,j] = np.real(ptn.mp.expectation(mps_obj, lat.get("n", ind)))
            
    mps_fl=mps_nm(params)
    # Saving the 2D array in a text file
    full_name = tar_loc+mps_fl[:-4]+"_nArr.txt"
    np.savetxt(full_name, n_arr)
    
    return full_name        


def cur_arr_save(mps_obj, tar_loc, params):
    '''
    mps_obj: MPS object
    tar_loc: location in which the resultant particle density array to be stored, target location
    params: SimpleNamespace object containing necessary parameters
    
    Stores numpy array for current operator in the target location
    Returns: the name of file in which the array is stored, convenient for extracting the file later
    '''
    #produce the correct lattice for mps
    lat=params.lat
    Lx=params.Lx; Ly=params.Ly; Nphi=params.Nphi; pbc=params.pbc
    
    #array for current
    cur_arr = np.zeros((Lx,Ly,2), dtype=complex)    #last index: x,y component of current
    
    for n in range(Lx):
        for m in range(Ly):
            if pbc:
                alpha = Nphi/(Lx-1)/Ly
            else:
                alpha = Nphi/(Lx-1)/(Ly-1)
            i = n*Ly + m   #site index of MPS, starting from 0
            t_g = np.exp(2j*np.pi*alpha*n)           #t=1
            if m!=Ly-1:
                #y-current
                cur_op=ptn.mp.dot(lat.get("c",i), lat.get("c",i+1))
                cur_arr[n,m,1]=ptn.mp.expectation(mps_obj, cur_op)*t_g
                
            elif params.pbc:
                #periodic y-current at the end
                cur_op = ptn.mp.dot(lat.get("c",i), lat.get("c",i-Ly+1))
                cur_arr[n,m,1]=ptn.mp.expectation(mps_obj, cur_op)*t_g
                
            if n!=Lx-1:
                #x-current
                cur_op= ptn.mp.dot(lat.get("c",i), lat.get("c",i+Ly)) 
                cur_arr[n,m,0]=ptn.mp.expectation(mps_obj, cur_op)
            else:
                cur_arr[n,m,0]=0*1j
                

    mps_fl=mps_nm(params)
    cur_arr_reshaped = cur_arr.reshape(cur_arr.shape[0], -1)
    # Saving the 2D array in a text file
    full_name = tar_loc+mps_fl[:-4]+"_cur_Arr.txt"
    np.savetxt(full_name, cur_arr_reshaped)
    
    return full_name
   
    
def nCorr_arr_save(mps_obj, tar_loc, params):
    '''
    mps_obj: MPS object
    tar_loc: location in which the resultant particle density array to be stored, target location
    params: SimpleNamespace object containing necessary parameters
    
    Stores numpy array for density-density correlation in the target location
    Returns: the name of file in which the array is stored, convenient for extracting the file later
    '''
    #produce the correct lattice for mps
    lat=params.lat
    Lx=params.Lx; Ly=params.Ly
    
    #density-density correlation
    y0=0
    x0_l = [0, Lx-1, Lx//2, Lx//2+1, Lx//2+2, Lx//2+3, Lx//2+4]
    nCorr_arr = np.zeros((len(x0_l),Lx,2))    #first index: 0, Lx, or center, #second index: correlaton functoin
                                      #third index: value of correlation function and distance, to be used for plots

    for k in range(nCorr_arr.shape[0]):
        x0 = x0_l[k]
        for i in range(nCorr_arr.shape[1]):
            ind_1 = Ly*x0+y0
            ind_2 = Ly*i+ y0

            #density-density
            n0_nr = ptn.mp.dot(lat.get("n", ind_1), lat.get("n", ind_2))
            nCorr_arr[k,i] = np.real(ptn.mp.expectation(mps_obj, n0_nr)), np.abs(x0-i)

    mps_fl=mps_nm(params)
    nCorr_arr_reshaped = nCorr_arr.reshape(nCorr_arr.shape[0], -1)
    # Saving the 2D array in a text file
    full_name = tar_loc+mps_fl[:-4]+"_nCorr_Arr.txt"
    np.savetxt(full_name, nCorr_arr_reshaped)
    
    return full_name

def sCorr_arr_save(mps_obj, tar_loc, params):
    '''
    mps_obj: MPS object
    tar_loc: location in which the resultant particle density array to be stored, target location
    params: SimpleNamespace object containing necessary parameters
    
    Stores numpy array for spin-spin correlation in the target location
    Returns: the name of file in which the array is stored, convenient for extracting the file later
    '''
    
    #produce the correct lattice for mps
    lat=params.lat
    Lx=params.Lx; Ly=params.Ly
    
    #spin-spin correlation functions at y=0
    y0=0
    x0_l = [0, Lx-1, Lx//2, Lx//2+1, Lx//2+2, Lx//2+3, Lx//2+4]
    sCorr_arr = np.zeros((len(x0_l),Lx,2))    #first index: 0, Lx, or center, #second index: correlaton functoin
                                      #third index: value of correlation function and distance, to be used for plots

    for k in range(sCorr_arr.shape[0]):
        x0 = x0_l[k]
        for i in range(sCorr_arr.shape[1]):
            ind_1 = Ly*x0+y0
            ind_2 = Ly*i+ y0

            #spin-spin
            s0_sr = ptn.mp.dot(lat.get("s", ind_1), lat.get("s", ind_2))
            sCorr_arr[k,i] = np.real(ptn.mp.expectation(mps_obj, s0_sr)), np.abs(x0-i)    #keep the distance for plotting later
    
    mps_fl=mps_nm(params)
    sCorr_arr_reshaped = sCorr_arr.reshape(sCorr_arr.shape[0], -1)           #change to 2d array, necessary for savetxt function.
    # Saving the 2D array in a text file
    full_name = tar_loc+mps_fl[:-4]+"_sCorr_Arr.txt"
    np.savetxt(full_name, sCorr_arr_reshaped)
    
    return full_name


def load_arr_high_bond(p1, tar_loc, chis, c, num_locs=3):
    '''
    returns arrays of physical quantities with highest bond dimension from a given folder
    c: choice of physical quantity, see the dictionary below
    chis:bond dimension sequence
    tar_loc: target folder to which the files to be stored
    '''
    
    funcs = {"n": ["_nArr.txt", p1.Lx, p1.Ly] , "nn": ["_nCorr_Arr.txt", num_locs, p1.Lx, 2], "ss": ["_sCorr_Arr.txt", num_locs, p1.Lx, 2]}
    
    choose = funcs[c]
    for i in range(1, len(chis)):
        p1.ind=len(chis)-i; p1.bond = chis[-i]
        fl_name = tar_loc+mps_nm(p1)[:-4] + choose[0]

        try:
            Corr_arr_rs = np.loadtxt(fl_name)
            Corr_Arr = Corr_arr_rs.reshape(choose[1:])
            break
            
        except:
            pass
            
    return Corr_Arr