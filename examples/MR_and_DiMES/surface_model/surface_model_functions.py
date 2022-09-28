import electronvolt as units
import io, libconf
import numpy as np
import matplotlib.pyplot as plt
import pyGITR
from pyGITR.Particles import *
import netCDF4
from netCDF4 import Dataset
import os
import math
from pyGITR.process import *
from pyGITR.process_functions import *
from pyGITR.Physical_Sputtering import *
from pyGITR.make_particleSource import *

import sys
sys.path.insert(0, 'setup/')
sys.path.append('.')
import input
from input import make_input
import plot_tracks
from plot_tracks import plotTracks

N_GITR = 10000 # number of GITR particles
FileNameSurfaceConcentration = 'input/surface_evolution_C_W.nc'
Delta_t = 0.1 # in seconds
Delta_t_gitr = 1e-10 # match this with the driver code
Delta_implant = 1e-5 # enter parameter value and units
Stopping_criteria = 0.02 # for C_C and C_W

Flux_H = 1.0e20
alpha_c = 0.02       # Carbon concentration in the background plasma
Flux_C = alpha_c*Flux_H


amu_C = 12 #for carbon
amu_W = 184 #for tungsten

n_atom = 6e22 # average number density
weight_gitr = Delta_t/Delta_t_gitr

Sputtering_yield_H_to_C = 0.005  # max 0.005
Sputtering_yield_C_to_C = 0.22   # 0.22 treated almost constant

Sputtering_yield_H_to_W = 0.002  # max 0.002 
Sputtering_yield_C_to_W = 0.5  # 0.5 treated almost constant

Reflection_yield_C_to_C = 0.005    # max 0.005 min 0.001  steady-state :  0.9
Reflection_yield_C_to_W = 0.75  # max 0.95 min 0.67    steady-state :  0.005


def getE(amu,vx,vy,vz):
    return np.array(0.5*amu*1.66e-27*(vx**2 + vy**2 + vz**2)/1.602e-19)


def getGeom(File):
    with io.open(File) as f:
        config = libconf.load(f)

    area = np.array(config.geom.area)
    surf = np.array(config.geom.surface)
    Z = np.array(config.geom.Z)

    return area,surf,Z

def makeInitNC(Surfaces,area,Conc):
    print("Creating initial {}...".format(FileNameSurfaceConcentration))
    initial_token_flux = 1.0e19  # tunable parameter
    dim = len(Surfaces)

    Flux_proportionality = {}
    for Z in Conc.keys():
        if Z==74:
            Flux_proportionality[Z] = 0
            for k in Surfaces:
                Flux_proportionality[Z] += initial_token_flux*area[k]*Delta_t_gitr/N_GITR

    Surface_time = np.full((1,1),0.0)
    Surface_number = np.array(range(dim))

    ncFile = netCDF4.Dataset(FileNameSurfaceConcentration, 'w', format='NETCDF4')
    s_number_dim = ncFile.createDimension('surface_dim', dim) # surface number dimension
    s_time_dim = ncFile.createDimension('time_dim', len(Surface_time)) # time dimension

    s_number = ncFile.createVariable('surface_number', np.float32, ('surface_dim',))
    s_time = ncFile.createVariable('time', np.float32, ('time_dim',))
    s_concentration = {}
    flux_proportionality = {}
    for Z in Conc.keys():
        s_concentration[Z] = ncFile.createVariable('surface_concentration_{}'.format(Z), np.float64, ('surface_dim','time_dim'))
        if Z == 74:
            flux_proportionality[Z] = ncFile.createVariable('Flux_Conversion_{}'.format(Z),np.float64,('time_dim'))

    s_number[:] = np.linspace(1,dim,dim)
    s_time[:] = Surface_time
    for Z in Conc.keys():
        s_concentration[Z][:,:] = Conc[Z]
        if Z == 74:
            flux_proportionality[Z][:] = Flux_proportionality[Z]

    ncFile.close()


def loadSurfaceConcentration(Zs):
    # Reading the surface features from the surface evolution netcdf file
    print("Loading {}...".format(FileNameSurfaceConcentration))
    SurfaceConcentrationData = netCDF4.Dataset(FileNameSurfaceConcentration, "r", format="NETCDF4")

    # Record concentrations of all surface elements and their initial Z
    Flux_proportionality = {}
    Conc = {}
    for z in Zs:
        Conc[z] = np.array(SurfaceConcentrationData['surface_concentration_{}'.format(z)][:,:]) # [:,:]
        if z == 74:
            Flux_proportionality[z] = np.array(SurfaceConcentrationData['Flux_Conversion_{}'.format(z)][:])

    Surface_time = np.array(SurfaceConcentrationData['time'][:])
    Surface_number = np.array(SurfaceConcentrationData['surface_number'][:])
    counter = len(Surface_time)-1

    SurfaceConcentrationData.close()

    return Conc,Flux_proportionality,Surface_time


def saveSurfaceConcentration(Surfaces,Surface_time,Conc,Flux_proportionality):
    #Writing the surface features with time
    print("Saving {}...".format(FileNameSurfaceConcentration))

    ncFile = netCDF4.Dataset(FileNameSurfaceConcentration, 'w', format='NETCDF4')
    s_number_dim = ncFile.createDimension('surface_dim', len(Surfaces)) # surface number dimension
    s_time_dim = ncFile.createDimension('time_dim', len(Surface_time)) # time dimension

    s_number = ncFile.createVariable('surface_number', np.float32, ('surface_dim',))
    s_time = ncFile.createVariable('time', np.float32, ('time_dim',))

    s_concentration = {}
    flux_proportionality = {}
    for Z in Conc.keys():
        s_concentration[Z] = ncFile.createVariable('surface_concentration_{}'.format(Z),np.float64,('surface_dim','time_dim'))
        if Z == 74:
            flux_proportionality[Z] = ncFile.createVariable('Flux_Conversion_{}'.format(Z),np.float64,('time_dim'))

    s_number[:] = np.linspace(1,len(Surfaces),len(Surfaces))
    s_time[:] = Surface_time

    for Z in Conc.keys():
        s_concentration[Z][:,:] = Conc[Z]
        if Z == 74:
            flux_proportionality[Z][:] = Flux_proportionality[Z]

    ncFile.close()


def getZsandSurfaces(surf,Z):
    # Make Zs and Surfaces arrays
    Zs = []
    Surfaces = []
    idx = np.arange(0,len(surf))
    for surface,z,i in zip(surf,Z,idx):
        if surface:
            Zs.append(z)
            Surfaces.append(i)
    Zs = np.unique(Zs)

    return Zs, Surfaces


def initConc(Zs,Surfaces,area,Z):
    Conc = {}
    for z in Zs:
        Conc[z] = np.full((len(Surfaces),1), 0.0)
        # step through the surfaces, set to 1.0 for those surfaces with Z
        for k,surface in enumerate(Surfaces):
            if Z[surface] == z:
                Conc[z][k] = 1.0

    # Create initial ncFile
    makeInitNC(Surfaces,area,Conc)


def calcEroDep(Zs,Surfaces,area,wGITR,Energy,surfaceHit):
    Conc,Flux_proportionality,Surface_time = loadSurfaceConcentration(Zs)

    # Tungsten loop
    Gamma_W_redep = np.zeros((len(Surfaces),1))
    Y_WW_Gamma_W_redep = np.zeros((len(Surfaces),1))
    # Y_WC_Gamma_W_redep = np.zeros((len(Surfaces),1))

    print("Index   Energy   Y_WW_Gamma_W_redep")
    if wGITR:
        for i in range(len(Energy)):
            if surfaceHit[i] != -1:
                surface_index = int(surfaceHit[i])
                sr_object = Sputtering_and_reflection()

                for j in Surfaces:
                    if j == surface_index:
                        Flux_W_local = Flux_proportionality[74][-1]/(Delta_t_gitr*area[surface_index])
                        Gamma_W_redep[j] += Flux_W_local  # check this
                        Y_WW_Gamma_W_redep[j] += sr_object.Calculate_PhysicalSputteringParameters('W','W',Energy[i])*Flux_W_local
                        # Y_WC_Gamma_W_redep[j] += sr_object.Calculate_PhysicalSputteringParameters('W','C',Energy[i])*Flux_W_local

                        # print(j,Energy[i],Y_WW_Gamma_W_redep[j])

    chi_W_ero =  Y_WW_Gamma_W_redep + Sputtering_yield_H_to_W*Flux_H + Sputtering_yield_C_to_W*Flux_C

    Gamma_W_ero_global = np.reshape(Conc[74][:,-1],(len(Surfaces),1)) * chi_W_ero
    Gamma_W_dep_global = Gamma_W_redep

    chi_C_ero_1 =  Sputtering_yield_H_to_C*Flux_H + Sputtering_yield_C_to_C*Flux_C   

    chi_C_dep_1 = np.zeros((len(Surfaces),1)) + (1-Reflection_yield_C_to_C)*Flux_C
    chi_C_dep_2 = np.zeros((len(Surfaces),1)) + (1-Reflection_yield_C_to_W)*Flux_C
    
    Gamma_C_ero_global = np.reshape(Conc[6][:,-1],(len(Surfaces),1)) * chi_C_ero_1
    Gamma_C_dep_global = np.reshape(Conc[6][:,-1],(len(Surfaces),1)) * chi_C_dep_1 +\
        np.reshape(Conc[74][:,-1],(len(Surfaces),1)) * chi_C_dep_2

    prop_W = 0
    for i,surface in enumerate(Surfaces):
        prop_W = prop_W + Gamma_W_ero_global[i]*area[surface]*Delta_t_gitr
    prop_W = prop_W/N_GITR

    return Conc, Flux_proportionality, Surface_time, chi_W_ero, chi_C_ero_1, chi_C_dep_1,chi_C_dep_2,Gamma_W_redep, Gamma_W_ero_global, prop_W


def particleSource(Surfaces,Z,Gamma_W_ero_global,area,prop_W,rank,GeomFile,plot):
    particleSourceDict = {}
    for i,surface in enumerate(Surfaces):
        if Z[i] == 74: # == specie
            num_particles = round(np.array(Gamma_W_ero_global[i]*Delta_t_gitr*area[surface]/prop_W).item())
            if num_particles!=0: 
                # print("Surface:",surface,"particles:",num_particles)
                particleSourceDict[surface] = num_particles

    number = 0
    for key in particleSourceDict.keys():
        number += particleSourceDict[key]

    print("Making particleSource_{}.nc".format(rank))
    makeParticleSource(particleSourceDict, GeomFile, "input/particleSource_{}.nc".format(rank),plot=plot)

    return number


def loadPositions(Surfaces):
    # print("Loading positions...")
    PositionsFile = "output/positions.nc"
    PositionData = netCDF4.Dataset(PositionsFile, "r", format="NETCDF4")
    dims, vars = getPositionsData(PositionData)

    surfaceHit = np.array(PositionData['surfaceHit'])
    vx = np.array(PositionData['vx'])
    vy = np.array(PositionData['vy'])
    vz = np.array(PositionData['vz'])
    Energy = getE(amu_W,vx,vy,vz)
    PositionData.close()

    W_elem = {}
    for surface_element in np.unique(surfaceHit):
        if surface_element in Surfaces:
            W_elem[surface_element] = 0

    for surface_element in surfaceHit:
        if surface_element in Surfaces:
            W_elem[surface_element] += 1
    print("Surface element: number of W hit")
    print(W_elem)

    return Energy, surfaceHit


##########################################################
#----------------- INIT PARTICLE SOURCE -----------------#
##########################################################

def InitParticleSource(rank,wGITR,plot=False):

    # Reading geomtry file
    GeomFile = "input/gitrGeom.cfg"
    area,surf,Z = getGeom(GeomFile)

    Zs, Surfaces = getZsandSurfaces(surf,Z)

    initConc(Zs,Surfaces,area,Z)

    Conc, Flux_proportionality, Surface_time, \
        chi_W_ero, chi_C_ero_1, chi_C_dep_1, chi_C_dep_2, \
            Gamma_W_redep, Gamma_W_ero_global, prop_W = calcEroDep(Zs,Surfaces,area,False,0,0)

    numParticles = particleSource(Surfaces,Z,Gamma_W_ero_global,area,prop_W,rank,GeomFile,plot)

    return Conc,Zs,Surfaces,Flux_proportionality,chi_W_ero,chi_C_ero_1,\
        chi_C_dep_1,chi_C_dep_2,Gamma_W_redep,Surface_time,prop_W,numParticles


##########################################################
#-------------------- PARTICLE SOURCE -------------------#
##########################################################

def ParticleSourceFromSurfaceModel(rank,wGITR,plot=False):

    # Reading geomtry file
    GeomFile = "input/gitrGeom.cfg"
    area,surf,Z = getGeom(GeomFile)

    Zs, Surfaces = getZsandSurfaces(surf,Z)

    # Reading position files of Tungsten
    if wGITR:
        Energy, surfaceHit = loadPositions(Surfaces)

        Conc, Flux_proportionality, Surface_time, \
            chi_W_ero, chi_C_ero_1, chi_C_dep_1, chi_C_dep_2, \
                Gamma_W_redep, Gamma_W_ero_global, prop_W = calcEroDep(Zs,Surfaces,area,wGITR,Energy,surfaceHit)

    else:
        Conc, Flux_proportionality, Surface_time, \
            chi_W_ero, chi_C_ero_1, chi_C_dep_1, chi_C_dep_2, \
                Gamma_W_redep, Gamma_W_ero_global, prop_W = calcEroDep(Zs,Surfaces,area,wGITR,0,0)

    numParticles = particleSource(Surfaces,Z,Gamma_W_ero_global,area,prop_W,rank,GeomFile,plot)

    return Conc,Zs,Surfaces,Flux_proportionality,chi_W_ero,chi_C_ero_1,\
        chi_C_dep_1,chi_C_dep_2,Gamma_W_redep,Surface_time,prop_W,numParticles


##########################################################
#--------------------- SURFACE MODEL --------------------#
##########################################################

def runSurfaceModel(Conc,Zs,Surfaces,chi_W_ero,chi_C_ero_1,\
    chi_C_dep_1,chi_C_dep_2,Gamma_W_redep,Flux_proportionality,prop_W,Surface_time):

    # Estimating the total time evolution for the surface model
    last_entry_C = np.reshape(Conc[6][:,-1],(len(Surfaces),1))
    last_entry_W = np.reshape(Conc[74][:,-1],(len(Surfaces),1))
    # print(last_entry_W)

    Gamma_W_ero = last_entry_W*chi_W_ero
    Gamma_C_ero = last_entry_C*chi_C_ero_1

    Gamma_C_dep = last_entry_C*chi_C_dep_1 + last_entry_W*chi_C_dep_2
    # Gamma_W_dep = Gamma_W_redep

    Gamma_C_net = Gamma_C_dep - Gamma_C_ero
    Gamma_W_net = -Gamma_W_ero

    Gamma_C_bulk = np.zeros((len(Surfaces),1))
    Gamma_W_bulk = np.zeros((len(Surfaces),1))

    for surface_index in range(len(Surfaces)):
        if (Gamma_C_net[surface_index] + Gamma_W_net[surface_index]) > 0: # deposition regime
            # print("deposition")
            Gamma_C_bulk[surface_index] = last_entry_C[surface_index,0]*(Gamma_C_net[surface_index]+Gamma_W_net[surface_index])
            Gamma_W_bulk[surface_index] = last_entry_W[surface_index,0]*(Gamma_C_net[surface_index]+Gamma_W_net[surface_index])
        else:  #  erosion regime
            # print("erosion")
            Gamma_C_bulk[surface_index] = 0
            Gamma_W_bulk[surface_index] = (Gamma_C_net[surface_index]+Gamma_W_net[surface_index])

    RHS_C = Gamma_C_net - Gamma_C_bulk
    RHS_W = Gamma_W_net - Gamma_W_bulk

    Delta_t_surface_estimate_C = (Delta_implant*n_atom*Stopping_criteria)/RHS_C
    Delta_t_surface_estimate_W = (Delta_implant*n_atom*Stopping_criteria)/RHS_W
    Delta_t_surface = min(np.amin(Delta_t_surface_estimate_C),np.amin(Delta_t_surface_estimate_C))
    print("Delta_t_surface:", Delta_t_surface)

    # Run surface model
    Time = Delta_t_surface
    Time_steps = 1e4
    Delta_Time = Delta_t/Time_steps # THIS IS THE TIMESTEP?
    Delta_t_Stopping = 0

    newConc = {}
    for z in Zs:
        newConc[z] = np.reshape(Conc[z][:,-1],(len(Surfaces),1))
            
    print("Running surface model...")
    for t in range(1,int(Time_steps)):

        # Incomplete        
        # Gamma_ero = {}
        # Gamma_dep = {}
        # for z in Zs:
        #     Gamma_ero[z] = newConc[z] * chi_ero[z] # C has 2 different chis
            
        #     if z == 74: Gamma_dep[z] = Gamma_redep[z]
        #     if z == 6: Gamma_dep[z] = newConc[z] * chi_dep[z]

        Gamma_W_ero = newConc[74] * chi_W_ero
        Gamma_C_ero = newConc[6] * chi_C_ero_1 
        Gamma_C_dep = newConc[6] * chi_C_dep_1 +newConc[74] * chi_C_dep_2
        Gamma_W_dep = Gamma_W_redep
        
        # determining erosion or deposition
        Gamma_C_net = Gamma_C_dep - Gamma_C_ero
        Gamma_W_net = -Gamma_W_ero
        
        Gamma_C_bulk = np.zeros((len(Surfaces),1))
        Gamma_W_bulk = np.zeros((len(Surfaces),1))
        
        for surface_index in range(len(Surfaces)):
            if (Gamma_C_net[surface_index] + Gamma_W_net[surface_index]) > 0: # deposition regime
                # print("deposition")
                Gamma_C_bulk[surface_index] = newConc[6][surface_index,0]*(Gamma_C_net[surface_index]+Gamma_W_net[surface_index])
                Gamma_W_bulk[surface_index] = newConc[74][surface_index,0]*(Gamma_C_net[surface_index]+Gamma_W_net[surface_index])
            else:  #  erosion regime
                # print("erosion")
                Gamma_C_bulk[surface_index] = 0
                Gamma_W_bulk[surface_index] = (Gamma_C_net[surface_index]+Gamma_W_net[surface_index])
        
        newConc[6] = newConc[6] + Delta_Time*(Gamma_C_net - Gamma_C_bulk)/(Delta_implant*n_atom)
        newConc[74] = newConc[74] + Delta_Time*(Gamma_W_net - Gamma_W_bulk)/(Delta_implant*n_atom)
        
        Delta_t_Stopping += Delta_Time
            
        if (np.abs(newConc[6]-last_entry_C)>Stopping_criteria).any() or (np.abs(newConc[74]-last_entry_W)>Stopping_criteria).any():
            print(Delta_t_Stopping," Delta_t_Stopping ", t)
            break

    # print(np.shape(Conc[74]),np.shape(newConc[74]))
    for z in Zs:
        Conc[z] = np.concatenate((Conc[z],newConc[z]),axis=1) # ADDS THE NEW ENTRY TO THE RHS ON CONC DICT

    for z in Zs:
        if z == 74:
            Flux_proportionality[z] = np.append(Flux_proportionality[z],(1/prop_W))

    Surface_time = np.append(Surface_time,Surface_time[-1]+Delta_t_Stopping)

    saveSurfaceConcentration(Surfaces,Surface_time,Conc,Flux_proportionality)

    return Conc, Delta_t_Stopping


def runGITR(rank,executable,folder,nP,dt,nT,plot):

    print("- GITR -")
    print(nP,dt,nT)
    make_input(nP,dt,nT,ParticleFile="particleSource_{}.nc".format(rank),folder=folder)
    print("\nRunning GITR, step {}".format(rank))
    os.system("./exe/{} > log{}".format(executable,rank))

    plotTracks(filename='output/history.nc',geomFile='input/gitrGeom.cfg',plasmaFile='input/profiles.nc',plot=plot)

    # save results
    # os.system("cp output/positions.nc output/positions_{}.nc".format(rank))
    os.system("cp output/surface.nc output/surface_{}.nc".format(rank))
    os.system("cp output/history.nc output/history_{}.nc".format(rank))
    # os.system("cp output/spec.nc output/spec_{}.nc".format(rank))
    os.system("cp {} input/surface_evolution_C_W_{}.nc".format(FileNameSurfaceConcentration,rank))

def plotConc(surface_index,Surface_time,Conc):
    plt.figure()
    lineStyle = {74:'d-k',6:'d-r'}
    for Z in Conc.keys():
        plt.plot(Surface_time, Conc[Z][surface_index,:], lineStyle[Z], label='Concentration of {}'.format(Z))
    plt.legend()
    plt.title("Surface Element {}".format(surface_index))
    plt.show()