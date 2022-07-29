#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generation of particles distribution for GITR.
@author: guterl
"""
import numpy as np
from numpy import random

from pyGITR.Particles import *
ParticleFile='particleConf.nc'
p = ParticleDistribution()

# Attributes of particles x,y,z,v,x,vy,vz are stored in the dictionary p.Particles
# Display the list of attributes of particles.
# p.ShowAttr()

# First, set numbers of particles. This is not affecting exisiting attributes
# but only the generation of those attributes
NP=100
p.SetAttr('Np', NP)

# Set distribution for attributes of particles
# First, show list of available distribution pdf
# p.ShowAvailablePdfs()

# Set positions of particles
xld = []
yld = []
while len(xld) != NP*0.5:
    xr = random.uniform(-0.005,0.005)
    yr = random.uniform(-0.005,0.005)
    if np.sqrt(xr*xr+yr*yr) <= 0.005:
        xld.append(xr)
        yld.append(yr)
xld = np.array(xld)
yld = np.array(yld)

xsd = []
ysd = []
while len(xsd) != NP*0.5:
    xr = random.uniform(-0.0105,-0.0095)
    yr = random.uniform(-0.0005,0.0005)
    if np.sqrt((xr+0.01)**2+yr*yr) <= 0.0005:
        xsd.append(xr)
        ysd.append(yr)
xsd = np.array(xsd)
ysd = np.array(ysd)

xs = np.concatenate((xld,xsd),axis=0)
ys = np.concatenate((yld,ysd),axis=0)

p.SetAttr('x',xs)
p.SetAttr('y',ys)
p.SetAttr('z',1e-5)

# Set velocities of particles
p.SetAttr(['vx','vy'],'Gaussian')
p.SetAttr(['vz'],'Gaussian')

# Rescale velocity by characteristic velocity
vpara = 1e4
vperp = 1e5
p.ScaleAttr(['vy','vz'],vperp)
p.ScaleAttr('vx',vpara)

# Rotate velocity of particles to put in the reference frame of the magnetic field

# Define rotation
thetaB = 2 #degree by default, set Degree=False in RotateVector for radian
AxisrotB = [0,1,0]

# Check results of the rotation
B0 = [1,0,0]
B = RotateVector(B0, AxisrotB, thetaB)
#PlotVector(B)

# Rotate velocity vector
p.Rotate('v', AxisrotB, thetaB) #rotate (vx,vy,vz) along AxisrotB of angle thetaB

# Write particle distribution in netcdf file
p.WriteParticleFile('particleConf.nc')

# Distribution can be also defined with a user-defined pdf. Must be formatted as
# f(x, **kwargs)
def LevyDistrib(x, c=1, mu=0):
    return np.sqrt(c/2/np.pi)*np.exp(-c/(x-mu))/((x-mu)**1.5)

p.SetAttr(['vz'],LevyDistrib, x=np.linspace(0.001,10,1000), c=2, mu=0)
p.ScaleAttr(['vz'],vperp)
p.PlotAttr('vz')
p.WriteParticleFile('particleConf2.nc')

os.system("mv particleConf* input")

FileNameParticleConf='input/particleConf2.nc'

import netCDF4
from netCDF4 import Dataset
ParticleConfData = Dataset(FileNameParticleConf, "r", format="NETCDF4")

import matplotlib.pyplot as plt
x = np.array(ParticleConfData['x'])
y = np.array(ParticleConfData['y'])
z = np.array(ParticleConfData['z'])
vx = np.array(ParticleConfData['vx'])
vy = np.array(ParticleConfData['vy'])
vz = np.array(ParticleConfData['vz'])

plt.close()
plt.scatter(x,y)
plt.xlabel('x')
plt.ylabel('y')
plt.title('x-y scatter')
plt.savefig('plots/x-y_scatter')

#plot Gaussian v dist
plt.close()
plt.hist(vx,bins=100)
plt.xlabel('Energy Bins [eV]')
plt.ylabel('Histogram')
plt.title('Thomson Energy Distribution')
plt.savefig('plots/vx_hist')

plt.close()
plt.hist(vz,bins=100)
plt.xlabel('Energy Bins [eV]')
plt.ylabel('Histogram')
plt.title('Thomson Energy Distribution')
plt.savefig('plots/vz_hist')
