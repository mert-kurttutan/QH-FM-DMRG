#!/usr/bin/env python3
## Example Python script calling DMRG.

import pyten as ptn
from .DMRG_lat import FHH_Ham_SU2, FHH_Ham_U1
from ..helpers import mps_nm, mps_load, n_arr_save, cur_arr_save
                
import numpy as np
import sys, time, csv

def run_dmrg_FHH_SU2(par, tar_folder):
    Lx = par.Lx; Ly = par.Ly   #Number of sites along x and y directions
    Nphi = par.Nphi
    U = par.U
    N = par.N
    S = par.S
    pbc = par.pbc
    #where files stored, e.g. tar_folder="/project/th-scratch/m/Mert.Kurttutan/QH_FM_02/Lx"+str(Lx)+"_Ly"+str(Ly) + "/"
    chis = par.chis        #bond dimensions for each stage
    sweep = par.sweep

    Q_nums = str(N) + " " + str(S)
    ##################
    ###### MAIN ######
    ##################

    ## path prefix
    pref = "/log-files"
    pref = tar_folder + pref  

    try:
        os.system("mkdir "+pref)
        os.system("mkdir "+pref+"log")
    except:
        pass


    print("Lx"+str(Lx)+"_Ly"+str(Ly)+"_Nphi"+str(Nphi)+"_U"+str(U)+"_N"+str(N)+"_S"+str(S)+"_PBC"+str(pbc))
    print("Generating lattice…")
    ## the lattice to be used
    lat = FHH_Ham_SU2(Ly, Lx, Nphi, 1.0, pbc)          #tperp=1.0

    print("Generating random state…")
    ## our initial random state, here generated with keyword arguments
    rnd = ptn.mp.generateCompleteState(lat, Q_nums)

    ## define Hamiltonians
    H = lat.get("Hj") + U*lat.get("Hu")
    lat.add("H", "full Hamiltonian", H)
    
    ## dmrg config object
    dmrgconf = ptn.dmrg.DMRGConfig()

    pre_str = "log/Lx"+str(Lx)+"_Ly"+str(Ly)+"_Nphi"+str(Nphi)+"_U"+str(U)+"_N"+str(N)+"_S"+str(S)+"_PBC"+str(pbc)
    pre_str = tar_folder + "log-files/" + pre_str

    ## prefix to be used for log files
    dmrgconf.prefix = pre_str

    for chi in chis:
        ## (m 100 x sweep[0])
        dmrgconf.stages += [ptn.dmrg.DMRGStage("(m "+str(chi)+" x "+ str(sweep[0]) +")")]
        dmrgconf.stages += [ptn.dmrg.DMRGStage("(m "+str(chi)+" x "+ str(sweep[1]) +" l 2 eb 0)")]

    ## set multi-threading
    ptn.threading.setTensorNum(4)

    ## set log-output
    ptn.setLogGLvl(0)
    ptn.setLogTLvl(0)

    ## PDMRG management object. Initialised with our random state, a list
    # of the desired Hamiltonians, the config object and a list of the
    # to-be-orthogonal states
    pdmrg = ptn.mp.dmrg.PDMRG(rnd, [lat.get("H")], dmrgconf)

    out_variance = "Lx"+str(Lx)+"_Ly"+str(Ly)+"_Nphi"+str(Nphi)+"_U"+str(U)+"_N"+str(N)+"_S"+str(S)+"_variance_FHH_SU2.dat"
    #mps_file = "Lx"+str(Lx)+"_Ly"+str(Ly)+"_Nphi"+str(Nphi)+"_U"+str(U)+"_N"+str(N)+"_S"+str(S)+"_PBC"+str(pbc)

    #mps_file = tar_folder  + mps_file
    out_variance = tar_folder + out_variance       #location if submitted via job

    e_new = 0
    for i in range(len(chis)):
        e_old = e_new
        par.ind = i; par.bond = chis[i]

        starttime = time.time()

        mps_0 = pdmrg.run()
        mps_tmp = pdmrg.run()

        if i > 6:
            mps_tmp.save(tar_folder + mps_nm(par))

        endtime = time.time()
        timediff = endtime - starttime

        e_new = ptn.mp.expectation(mps_tmp, lat.get("H"))
        esq = ptn.mp.expectation(mps_tmp, lat.get("H")*lat.get("H"))
        var = abs(esq - e_new**2)
        print("E = ", e_new)
        print("Δ = ", e_new - e_old)
        print("Var = ", var)


        f = open(out_variance, 'a')
        writer = csv.writer(f, delimiter=',')
        writer.writerow([Lx, Ly, Nphi, U, N, S, str(pbc), chis[i], var, np.real(e_new), np.real(e_new - e_old), timediff])
        f.close()



def conv_FHH_SU2_n(par, tar_loc, src_folder):
    '''
    Calculates the particle density and current density for states of parameter object par,
    Used for ensuring covergence
    '''
    
    source = src_folder + "Lx" + str(par.Lx) + "_Ly" + str(par.Ly) + "/"
    
    for i in range(7,len(par.chis)):
        par.bond=par.chis[i]; par.ind=i
        par.source = source

        try:
            mps_obj = mps_load(par)
            file_nm=n_arr_save(mps_obj, tar_loc, par)        #save it in the local dir

        except:
            print("State with " + "m_B=" + str(par.bond) + " is not produced")
                        
                                     
def conv_FHH_SU2_cur(par, tar_loc, src_folder):
    '''
    Calculates the particle density and current density for states of parameter object par,
    Used for ensuring covergence
    '''

    source = src_folder + "Lx" + str(par.Lx) + "_Ly" + str(par.Ly) + "/"
    par.source = source
    for i in range(7,len(par.chis)):
        par.bond=par.chis[i]; par.ind=i
        
        try:
            mps_obj = mps_load(par)
            file_nm=cur_arr_save(mps_obj, tar_loc, par)        #save it in the local dir

        except:
            print("State with " + "m_B=" + str(par.bond) + " is not produced")
                        
                        
                