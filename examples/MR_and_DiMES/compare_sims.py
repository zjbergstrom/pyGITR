import matplotlib.pyplot as plt
import netCDF4
from netCDF4 import Dataset
import numpy as np
# import sys
# sys.path.insert(0, 'surface_model/')
# sys.path.append('.')
# import surface_model_functions
# from surface_model_functions import *

noGITR = "input/surface_evolution_C_W_noGITR_1to10.nc"
wGITR = "input/surface_evolution_C_W.nc"
Zs = [6,74]


plt.figure()
lineStyle1 = {74:'d-k',6:'d-r'}
lineStyle2 = {74:'o--k',6:'o--r'}
surface_index = 90

print("Loading {}".format(noGITR))
noGITRData = netCDF4.Dataset(noGITR, "r", format="NETCDF4")
Flux_proportionality = {}
Conc = {}
for z in Zs:
    Conc[z] = np.array(noGITRData['surface_concentration_{}'.format(z)][:,:]) # [:,:]
    if z == 74:
        Flux_proportionality[z] = np.array(noGITRData['Flux_Conversion_{}'.format(z)][:])

Surface_time = np.array(noGITRData['time'][:])
Surface_number = np.array(noGITRData['surface_number'][:])
counter = len(Surface_time)-1
noGITRData.close()

for Z in Conc.keys():
    plt.plot(Surface_time, Conc[Z][surface_index,:], lineStyle1[Z], label='C$_{}$'.format(Z))

print(Conc[74][82:-1])

print("Loading {}".format(wGITR))
wGITRData = netCDF4.Dataset(wGITR, "r", format="NETCDF4")
Flux_proportionality = {}
Conc = {}
for z in Zs:
    Conc[z] = np.array(wGITRData['surface_concentration_{}'.format(z)][:,:]) # [:,:]
    if z == 74:
        Flux_proportionality[z] = np.array(wGITRData['Flux_Conversion_{}'.format(z)][:])

Surface_time = np.array(wGITRData['time'][:])
Surface_number = np.array(wGITRData['surface_number'][:])
counter = len(Surface_time)-1
wGITRData.close()

for Z in Conc.keys():
    plt.plot(Surface_time, Conc[Z][surface_index,:], lineStyle2[Z], label='C$_{}$'.format(Z))

print(Conc[74][82:-1])

plt.legend()
plt.title("Surface Element {}".format(surface_index))
plt.show()