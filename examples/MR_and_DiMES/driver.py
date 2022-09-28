#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 13:39:30 2022

@author: de, bergstrom
"""
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, 'surface_model/')
sys.path.append('.')
import surface_model_functions
from surface_model_functions import *


# GITR = "GITR_except_efield"
GITR = "GITR_final"
wGITR = True

species = [74] # List of GITR impurity species
dt = 1e-10
nT = 1e4

print("GITR exe:",GITR)
print("Run surface model w/GITR:",wGITR)
print("GITR species:",species[0])
print("Number of GITR particles",N_GITR) # number of GITR particles
print("GITR timestep:",dt) # match this with the driver code
print("GITR steps",nT)
print("Stopping criteria",Stopping_criteria) # for C_C and C_W


# Make initial particle source for GITR and
# initial concentration file for the surface model

# rank = 0
# Conc,Zs,Surfaces,Flux_proportionality,chi_W_ero,chi_C_ero_1,chi_C_dep_1,\
#     chi_C_dep_2,Gamma_W_redep,Surface_time,prop_W,numP = InitParticleSource(rank=rank,wGITR=wGITR,plot=False)

# nP = numP
# runGITR(rank,GITR,"input/",nP,dt,nT,0)

# Loop through GITR and surface model
start = 3
finish = 4
for rank in range(start, finish):
    # Run below if NOT starting from rank 0
    Conc,Zs,Surfaces,Flux_proportionality,chi_W_ero,chi_C_ero_1,chi_C_dep_1,\
        chi_C_dep_2,Gamma_W_redep,Surface_time,prop_W,numP = \
            ParticleSourceFromSurfaceModel(rank=rank,wGITR=wGITR,plot=False)

    # Load results of GITR simulation and run surface model
    Conc, Delta_t_Stopping = runSurfaceModel(Conc,Zs,Surfaces,chi_W_ero,\
        chi_C_ero_1,chi_C_dep_1,chi_C_dep_2,Gamma_W_redep,Flux_proportionality,prop_W,Surface_time)

    # Make particle source from the output of the surface model
    Conc,Zs,Surfaces,Flux_proportionality,chi_W_ero,chi_C_ero_1,chi_C_dep_1,\
        chi_C_dep_2,Gamma_W_redep,Surface_time,prop_W,numP = \
            ParticleSourceFromSurfaceModel(rank=rank,wGITR=wGITR,plot=False)
    
    # Run GITR simulation
    nP = numP
    runGITR(rank,GITR,"input/",nP,dt,nT,0)


# for surface in Surfaces:
#     if surface == 84:
#         plotConc(surface,Surface_time,Conc)

