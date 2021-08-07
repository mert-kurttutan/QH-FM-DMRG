#!/usr/bin/env python3
import pyten as ptn
import sys,os
sys.path.append(os.path.abspath(os.path.dirname("../src"))) 
from src import helpers, utils, DMRG 

'''
Calculates the current density with parameter given below,
Used for ensuring covergence
'''
Lx = int(float(sys.argv[1]))
Ly = int(float(sys.argv[2]))   #number of sites along y-direction
Nphi = float(sys.argv[3])
U = float(sys.argv[4])
N = int(float(sys.argv[5]))
S = float(sys. argv[6])
pbc = bool(float(sys.argv[7]))

#load the particle density to txt file, and then to numpy array, then plot as a function of x
chis = [100, 100, 200, 200, 400, 400, 400, 800, 800, 800, 800, 1600, 1600, 1600, 2000, 2000, 2000, 4000, 4000, 4000, 4000, 4000, 4000]

p1=helpers.params()
p1.folder="QH-FM-02"; p1.pbc=True; 
#tar_loc="/project/th-scratch/m/Mert.Kurttutan/QH-FM-02/corr-func/"     #necessary when submitting a job
#tar_loc="../data/dat-files-02-expc-arr/cur-arr/"
tar_loc="../temp/"
src_folder = "/project/th-scratch/m/Mert.Kurttutan/QH-FM-01/"

p1.Lx=Lx; p1.Ly=Ly; p1.Nphi=Nphi; p1.U=U; p1.N=int(N); p1.S=S;
p1.lat=ptn.mp.lat.su2u1.genFermiHubbardSpinCharge(p1.Lx*p1.Ly)
p1.source = src_folder + "Lx" + str(Lx) + "_Ly" + str(Ly) + "/"
for i in range(1,len(chis)+1):
    p1.bond=chis[-i]; p1.ind=len(chis)-i
                    
    try:
        mps_obj = helpers.mps_load(p1)
        file_nm=helpers.sCorr_arr_save(mps_obj, tar_loc, p1)        #save it in the local dir
        break                                                  #calculate it in the highest bond dim
                    
    except:
        print("State with " + "m_B=" + str(p1.bond) + " is not produced")
                        
                    
                