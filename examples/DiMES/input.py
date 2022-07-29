# -*- coding: utf-8 -*-

import os
from pyGITR import Input

ParticleFile='particleConf2.nc'
GeometryFile='gitrGeom.cfg'
B0 = 2.25
thetaB = 2
phiB = 0
NP = 100

# Initiallize input object
i = Input()

# Add structures to configuration file
i.SetBField(B0=2.25, theta = thetaB, phi = phiB)
i.SetTimeStep(dt=1e-11,nT=1e5)
i.SetGeometryFile(GeometryFile)
i.SetParticleSource(ParticleFile, nP=NP, Zmax=74, M=183, Z=4)
i.SetSurfaces()
i.SetDiagnostics()
i.SetBackgroundPlasmaProfiles()
i.SetSurfaceModel()

# i.SetGeomHash()
# i.Input['flags']['GEOM_HASH'] = 1
# i.SetGeomSheath()
# i.Input['flags']['GEOM_HASH_SHEATH'] = 1

i.SetConnectionLength()
i.SetForceEvaluation()

# Set the standard options
i.Input['flags']['USESURFACEMODEL'] = 1
i.Input['flags']['USECOULOMBCOLLISIONS'] = 1
i.Input['flags']['USEFRICTION'] = 1
i.Input['flags']['USEANGLESCATTERING'] = 1
i.Input['flags']['USEHEATING'] = 1
i.Input['impurityParticleSource']['Z'] = 74
i.Input['impurityParticleSource']['source_material_Z'] = 74
i.Input['backgroundPlasmaProfiles']['FlowVelocity']['flowVz'] = -2000

# Write input file
i.WriteInputFile()

os.system("mv gitrInput.cfg input")