import pyten as ptn
import numpy as np
import sys, csv, os
from tabulate import tabulate
from .proc_helpers import name_to_S, name_to_N, name_to_Nphi
from ..utils import getListOfFiles, read_last_line

def print_DMRG_table(params):
    '''
    prints: table for summary of DMRG computation, indicating several properties, see the list headers, below
    Files from DMRG calculation is gathered in loc_dir
    '''
    Lx=params.Lx; Ly=params.Ly; Nphi=params.Nphi; U=params.U; N=params.N; S=params.S
    #extract energy variance from txt files of DMRG
    target = params.target 
    if not params.pin:
        var_fl = "Lx"+str(Lx)+"_Ly"+str(Ly)+"_Nphi"+str(Nphi)+"_U"+str(U)+"_N"+str(N)+"_S"+str(S)+"_variance_FHH_SU2.dat"
    else:
        var_fl = "Lx"+str(Lx)+"_Ly"+str(Ly)+"_Nphi"+str(Nphi)+"_U"+str(U)+"_N"+str(N)+"_S"+str(S)+"_g" + str(params.g)+"_variance_FHH_SU2.dat"

    print("DMRG Table for " + var_fl[:-21] + " - 40 sweeps for each bond dim - " + "80 sweeps for 2 stages" + "\n")
    
    if not params.pin:
        num_list= [8, 7, -2, 9, -1]
    else:
        num_list= [9, 8, -2, 10, -1]

    with open(target+var_fl, 'r') as DMRGtxt:

        content = np.array(list(csv.reader(DMRGtxt, delimiter=',', quotechar='|')))

        print(tabulate(content[:,num_list], headers=['Variance', 'Bond Dimension', 'Energy Difference', 'Energy', 'Time for one stage'], tablefmt='orgtbl', floatfmt=".10f"))
   

    return content[:, num_list]


def load_energies(params, tar_loc):
    '''
    Lx: The length of the system along x-direction
    Ly: The length in y-direction
    tar_loc: where to look for the .dat files

    returns: an array of energy values for each particle filling, magnetic filling and spin polarization
    This is to be used when plottin phase diagram
    '''

    fileList = getListOfFiles(tar_loc)
    energy_list = []
    energy=0
    name = "Lx" + str(params.Lx) + "_Ly" + str(params.Ly)
    for file in fileList:
        if os.path.basename(file)[:len(name)] == name and file[-3:] == "dat":
            last_line_list = read_last_line(file)
            idx = 10 if params.pin else 9
            energy = last_line_list[idx]
            dnmr = (params.Lx-1)*params.Ly
            energy = float(energy)
            S = name_to_S(file); N = name_to_N(file); Nphi = name_to_Nphi(file)
            energy_list.append([energy/dnmr, Nphi, N/dnmr, S])
            
    return energy_list

    #Lx16_Ly5_Nphi34.0_U8.0_N34_S0.0_variance_FHH_SU2.dat

def phase_Diag(Lx, Ly):


    return