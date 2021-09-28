import pyten as ptn
import numpy as np
import sys, csv
from tabulate import tabulate

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