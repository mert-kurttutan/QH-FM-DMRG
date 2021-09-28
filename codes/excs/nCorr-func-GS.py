#!/usr/bin/env python3
import pyten as ptn
import sys,os
import traceback
pp=os.path.dirname(os.path.abspath(__file__))
pp = os.path.dirname(pp)
sys.path.append(pp)          #used this, be careful
from src import helpers, utils, DMRG 

'''
Calculates the current density with parameter given below,
Used for ensuring covergence
'''
par=helpers.params()

par.Lx = int(float(sys.argv[1]))
par.Ly = int(float(sys.argv[2]))   #number of sites along y-direction
par.Nphi = float(sys.argv[3])
par.U = float(sys.argv[4])
par.N = int(float(sys.argv[5]))
par.S = float(sys.argv[6])
par.pbc = bool(float(sys.argv[7]))
par.g = float(sys.argv[8])
par.x0 = int(float(sys.argv[9])); par.y0 = int(float(sys.argv[10]))
par.chis = [100, 100, 200, 200, 400, 400, 400, 800, 800, 800, 800, 1600, 1600, 1600, 2000, 2000, 2000, 4000, 4000, 4000, 4000, 4000, 4000]
lat=DMRG.FHH_Ham_SU2(par.Ly, par.Lx, par.Nphi, 1.0, par.pbc)  
par.lat = DMRG.add_pin_SU2(lat, par.g, par.x0, par.y0, par)
par.pin = True

#tar_loc = "/project/th-scratch/m/Mert.Kurttutan/QH-FM-01/Lx" + str(par.Lx) + "_Ly" + str(par.Ly) + "/"
#tar_loc = "/project/th-scratch/m/Mert.Kurttutan/QH-FM-02/Lx" + str(par.Lx) + "_Ly" + str(par.Ly) + "/"
tar_loc = "/project/th-scratch/m/Mert.Kurttutan/QH-FM-03/Lx" + str(par.Lx) + "_Ly" + str(par.Ly) + "/"
#tar_loc="../temp/"
#tar_loc= pp+"/temp_trial/"
src_folder = "/project/th-scratch/m/Mert.Kurttutan/QH-FM-03/"

par.lat=ptn.mp.lat.su2u1.genFermiHubbardSpinCharge(par.Lx*par.Ly)
par.source = src_folder + "Lx" + str(par.Lx) + "_Ly" + str(par.Ly) + "/"
for i in range(1,len(par.chis)+1):
    par.bond=par.chis[-i]; par.ind=len(par.chis)-i
                    
    try:
        mps_obj = helpers.mps_load(par)
        file_nm=helpers.nCorr_arr_save(mps_obj, tar_loc, par)        #save it in the local dir
        break                                                  #calculate it in the highest bond dim
                    
    except Exception:
        #traceback.print_exc()
        print("State with " + "m_B=" + str(par.bond) + " is not produced")
                        
                    
                