#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 17:53:18 2021

@author: jguterl
"""
import matplotlib.pyplot as plt

import electronvolt as units
import os
os.system("clear")

import numpy as np
FileNameSurface='output/surface.nc'
FileNameParticle='output/particleSource.nc'
FileNameHistory='output/history.nc'

import netCDF4
from netCDF4 import Dataset

SurfaceData = Dataset(FileNameSurface, "r", format="NETCDF4")
ParticleData = Dataset(FileNameParticle, "r", format="NETCDF4")
HistoryData = Dataset(FileNameHistory, "r", format="NETCDF4")

# print("SurfaceData")
# for dim in SurfaceData.dimensions.values():
#     print(dim)
# print("ParticleData")
# for dim in ParticleData.dimensions.values():
#     print(dim)
# print("HistoryData")
# for dim in HistoryData.dimensions.values():
#     print(dim)

print("SurfaceData")
for var in SurfaceData.variables.values():
    print(var)
# print("ParticleData")
# for var in ParticleData.variables.values():
#     print(var)
# print("HistoryData")
# for var in HistoryData.variables.values():
#     print(var)

D={}
vx = np.array(ParticleData['vx'])
vy = np.array(ParticleData['vy'])
vz = np.array(ParticleData['vz'])
#amu = np.array(ParticleData.variables['amu'])
amu = np.array([2])*vx.size
E0=1/2*amu*(vz**2+vy**2+vx**2)*units.mp/units.eV*(units.V)**2

Distrib = np.array(SurfaceData['surfEDist'])
GridE = np.array(SurfaceData['gridE'])
GridA = np.array(SurfaceData['gridA'])
print('Total particle={}'.format(np.sum(Distrib,(0,1,2))))

E=np.sum(Distrib[:,:,:],(0,2))
print('Total particle={}'.format(np.sum(Distrib,(0,1,2))))
plt.figure()
plt.plot(GridE,E)
plt.show()
# plt.figure
# plt.hist(E0,GridE)
# plt.show()
