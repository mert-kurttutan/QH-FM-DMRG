#!/usr/bin/env python3

import pyten as ptn
import numpy as np

def FHH_Ham_SU2(legs, length, Nphi, tperp, PBC=False):
    '''
    
    legs: number of sites in y-direction,
    length: number of sites in x-drection
    Nphi: total flux [in units of flux quanta]
    tperp: in-plane hopping parameter
    PBC: periodic boundar condition imposed along y-direction, i.e. cylinder geometry
    
    returns: SyTeN lattice object with hamiltonian terms attached to it.
    Associated MPS of output has the snake-like shape used for 2D systems, see the notes
    
    '''
    #width = legs +1  if PBC=False
    #width = legs  if PBC=True
    if PBC:
        alpha = Nphi/(legs*(length-1))

    if not PBC:
        alpha = Nphi/((legs-1)*(length-1))
    
    
    lat = ptn.mp.lat.su2u1.genFermiHubbardSpinCharge(legs*length)      #using su2u1
    
    
    onsite = [lat.get("nb",i) for i in range(length*legs)]  #Hubbard onsite Interaction term, double occupancy
    hopping = []
    
    for n in range(length):
        t_g = tperp*np.exp(2j*np.pi*alpha*n)           #hopping term with Landau Gauge
        t_gc = tperp*np.exp(-2j*np.pi*alpha*n)        #complex conjugate of it
        for m in range(legs):
            i = n*legs + m   #site index of MPS, starting from 0
            
            if m!=legs-1:
                #y-hopping
                hopping += [-t_g * ptn.mp.dot(lat.get("c",i), lat.get("c",i+1)) 
                            -t_gc * ptn.mp.dot(lat.get("c",i+1), lat.get("c",i))]
            elif PBC:
                #periodic y-hopping at the end
                hopping += [-t_g * ptn.mp.dot(lat.get("c",i), lat.get("c",i-legs+1)) 
                            -t_gc * ptn.mp.dot(lat.get("c",i-legs+1), lat.get("c",i))]
            if n!=length-1:
                #x-hopping
                hopping += [-tperp * ptn.mp.dot(lat.get("c",i), lat.get("c",i+legs)) 
                            -tperp * ptn.mp.dot(lat.get("c",i+legs), lat.get("c",i))]         
            
    

    hJ = ptn.mp.addLog(hopping)        #a way to add MPOS inside a list
    hU = ptn.mp.addLog(onsite)

    lat.add("Hj", "Hopping term (including -!)", hJ)
    lat.add("Hu", "Hubbard Interaction (with double occupancy)", hU)

    # print("Generated a MP lattice and Hamiltonian for a "+str(legs)+"-leg ladder of length "+str(length)
    #                     +" at alpha="+str(alpha)+" and rung hopping t_perp="+str(tperp))

    return lat



import pyten as ptn
import numpy as np

def FHH_Ham_U1(legs, length, Nphi, tperp, PBC=False):
    '''
    
    legs: number site along y-direction,
    length: number of sites in x-drection
    Nphi: total flux [in units of flux quanta]
    tperp: in-plane hopping parameter
    PBC: periodic boundar condition imposed along y-direction -> cylinder geometry
    
    returns: SyTeN lattice object with hamiltonian terms attached to it.
    Associated MPS of output has the snake-like shape used for 2D systems, see the notes
    
    '''
    #width = legs +1  if PBC=False
    #width = legs  if PBC=True
    if PBC:
        alpha = Nphi/(legs*(length-1))

    if not PBC:
        alpha = Nphi/((legs-1)*(length-1))
    
    lat = ptn.mp.lat.u1u1.genFermiHubbard(legs*length)      #using u1u1
    
    
    onsite = [lat.get("nd",i)*lat.get("nu",i) for i in range(length*legs)]  #Hubbard onsite Interaction term
    hopping = []
    
    for n in range(length):
        t_g = tperp*np.exp(2j*np.pi*alpha*n)           #hopping term with Landau Gauge, x-dependent
        t_gc = tperp*np.exp(-2j*np.pi*alpha*n)        #complex conjugate of it
        for m in range(legs):
            i = n*legs + m   #site index of MPS, starting from 0
            
            if m!=legs-1:
                #y-hopping
                hopping += [-t_g*lat.get("cu",i)*lat.get("chu",i+1) - t_gc*lat.get("cu",i+1)*lat.get("chu",i)   
                            - t_g*lat.get("cd",i)*lat.get("chd",i+1) - t_gc*lat.get("cd",i+1)*lat.get("chd",i)]
            elif PBC:
                #periodic y-hopping at the end
                hopping += [-t_g*lat.get("cu",i)*lat.get("chu",i-legs+1) - t_gc*lat.get("cu",i-legs+1)*lat.get("chu",i)  
                            - t_g*lat.get("cd",i)*lat.get("chd",i-legs+1) - t_gc*lat.get("cd",i-legs+1)*lat.get("chd",i)]
            if n!=length-1:
                #x-hopping
                hopping += [-tperp*lat.get("cu",i)*lat.get("chu",i+legs) - tperp*lat.get("cu",i+legs)*lat.get("chu",i)
                            - tperp*lat.get("cd",i)*lat.get("chd",i+legs) - tperp*lat.get("cd",i+legs)*lat.get("chd",i)]           
            
    

    hJ = ptn.mp.addLog(hopping)        #a way to add MPOS inside a list
    hU = ptn.mp.addLog(onsite)

    lat.add("Hj", "Hopping term (including -!)", hJ)
    lat.add("Hu", "Hubbard Interaction", hU)

    # print("Generated a MP lattice and Hamiltonian for a "+str(legs)+"-leg ladder of length "+str(length)
    #                     +" at alpha="+str(alpha)+" and rung hopping t_perp="+str(tperp))

    return lat


def add_pin_SU2(lat, g, x, y, par):
    '''
    Adds a pinning potential term to a given SU2-U1 fermionic lattice
    '''
    Lx=par.Lx; Ly=par.Ly
    
    
    i = x*(Ly) + y          #both x,y starts from 0
    h_pin = g*lat.get("n", i)
    
    lat.add("H_pin", "Pinning potential - for localized qp", h_pin)
    
    return lat
    
    