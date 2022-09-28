# -*- coding: utf-8 -*-

import os
from pyGITR import Input

ParticleFile='particleConf2.nc'
GeometryFile='gitrGeom.cfg'

B0 = 2.25
thetaB = 2
phiB = 0
nP = 100

# Initiallize input object
i = Input()

# Add structures to configuration file
i.SetBField(B0=2.25, theta = thetaB, phi = phiB)
i.SetTimeStep(dt=1e-9,nT=1e4)
i.SetGeometryFile(GeometryFile)
i.SetParticleSource(ParticleFile, nP=nP, Zmax=74, M=183, Z=4)
i.SetSurfaces()
i.SetDiagnostics()
i.SetBackgroundPlasmaProfiles()
i.SetSurfaceModel()
# i.SetGeomHash()
# i.SetConnectionLength()
i.SetForceEvaluation()
# i.SetGeomSheath()

# i.Input['flags']['GEOM_HASH_SHEATH'] = 1
# i.Input['flags']['GEOM_HASH'] = 1

i.Input['backgroundPlasmaProfiles']['Bfield']['interpolation'] = 0 # in fields.cpp, only needs to be >0
i.Input['flags']['USESHEATHEFIELD'] = 0
i.Input['flags']['USEPERPDIFFUSION'] = 0
i.Input['flags']['BFIELD_INTERP'] = 0
i.Input['flags']['GENERATE_LC'] = 0
i.Input['flags']['LC_INTERP'] = 0 # 0, 0 no Lc, 2 get from file, 3 consider flowV
i.Input['flags']['EFIELD_INTERP'] = 0
i.Input['flags']['PRESHEATH_INTERP'] = 0 # 3
i.Input['flags']['DENSITY_INTERP'] = 0
i.Input['flags']['TEMP_INTERP'] = 0
i.Input['flags']['FLOWV_INTERP'] = 0 # 3, 0 const from file, 1 Lc based, 2 2D, 3 3D
i.Input['flags']['GRADT_INTERP'] = 0 # 3, 1 R, 2 R+Z, 3 R+Z+Y
i.Input['flags']['USECYLSYMM'] = 0 # rotates the plasma with cylindrical geometry



# Set the standard options
i.Input['flags']['USESURFACEMODEL'] = 1
i.Input['flags']['USECOULOMBCOLLISIONS'] = 1
i.Input['flags']['USEFRICTION'] = 1
i.Input['flags']['USEANGLESCATTERING'] = 1
i.Input['flags']['USEHEATING'] = 1

i.Input['impurityParticleSource']['Z'] = 74
i.Input['impurityParticleSource']['source_material_Z'] = 74
# i.Input['backgroundPlasmaProfiles']['FlowVelocity']['flowVz'] = -2000
i.Input['backgroundPlasmaProfiles']['Bfield']['rString'] = 'br'
i.Input['backgroundPlasmaProfiles']['Bfield']['zString'] = 'bz'
i.Input['backgroundPlasmaProfiles']['Bfield']['yString'] = 'bt'

# Write input file
i.WriteInputFile(Folder='input',OverWrite=True)

