#!/usr/bin/env python3
import pyten as ptn
import sys,os
pp=os.path.dirname(os.path.abspath(__file__))
pp = os.path.dirname(pp)
sys.path.append(pp)          #used this, be careful
from src import helpers, utils, DMRG   

par=helpers.params()

par.Lx = int(float(sys.argv[1]))
par.Ly = int(float(sys.argv[2]))   #number of sites along y-direction
par.Nphi = float(sys.argv[3])
par.U = float(sys.argv[4])
par.N = int(float(sys.argv[5]))
par.S = float(sys. argv[6])
par.pbc = bool(float(sys.argv[7]))
par.chis = [100, 100, 200, 200, 400, 400, 400, 800, 800, 800, 800, 1600, 1600, 1600, 2000, 2000, 2000, 4000, 4000, 4000, 4000, 4000, 4000]
par.lat=ptn.mp.lat.su2u1.genFermiHubbardSpinCharge(par.Lx*par.Ly)

#tar_loc = "/project/th-scratch/m/Mert.Kurttutan/QH-FM-01/convergence/n-arr/"
#tar_loc = "../data/dat-files-01-expc-arr/n-arr/"
tar_loc= pp+"/temp_trial/"
src_folder = "/project/th-scratch/m/Mert.Kurttutan/QH-FM-01/"

if __name__ == "__main__":            #make sure it runs only when it is executed
    DMRG.conv_FHH_SU2_cur(par, tar_loc, src_folder)