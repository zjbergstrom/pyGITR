# -*- coding: utf-8 -*-

import os
from pyGITR import Input

def make_input(nP,dt,nT,ParticleFile='particleConf2.nc',folder='input/'):
    GeometryFile='gitrGeom.cfg'

    B0 = -2.25
    thetaB = 2
    phiB = 0

    # Initiallize input object
    i = Input()

    # Add structures to configuration file
    i.SetBField(B0=B0, theta = thetaB, phi = phiB)
    i.SetTimeStep(dt=dt,nT=nT)
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

    # Set the standard flags
    i.Input['flags']['BIASED_SURFACE'] = 1
    i.Input['flags']['USE_SURFACE_POTENTIAL'] = 1
    i.Input['flags']['USETHERMALFORCE'] = 1
    i.Input['flags']['USEPERPDIFFUSION'] = 1
    i.Input['flags']['USESURFACEMODEL'] = 1
    i.Input['flags']['USESHEATHEFIELD'] = 1 # 1 or 0, on or off
    i.Input['flags']['USEPRESHEATHEFIELD'] = 1  # 1 or 0, on or off
    i.Input['flags']['USECOULOMBCOLLISIONS'] = 1
    i.Input['flags']['USEFRICTION'] = 1
    i.Input['flags']['USEANGLESCATTERING'] = 1
    i.Input['flags']['USEHEATING'] = 1
    i.Input['flags']['FORCE_EVAL'] = 0
    i.Input['flags']['USECYLSYMM'] = 1 # rotates the plasma with cylindrical geometry

    # Set the INTERP flags
    i.Input['flags']['BFIELD_INTERP'] = 2
    i.Input['flags']['PRESHEATH_INTERP'] = 2 # 3
    i.Input['flags']['DENSITY_INTERP'] = 2
    i.Input['flags']['TEMP_INTERP'] = 2
    i.Input['flags']['FLOWV_INTERP'] = 2 # 3, 0 const from file, 1 Lc based, 2 2D, 3 3D
    i.Input['flags']['GRADT_INTERP'] = 2 # 3, 1 R, 2 R+Z, 3 R+Z+Y



    # Set other options
    i.Input['backgroundPlasmaProfiles']['Bfield']['interpolation'] = 1 # in fields.cpp, only needs to be >0
    i.Input['impurityParticleSource']['Z'] = 74
    i.Input['impurityParticleSource']['source_material_Z'] = 74
    i.Input['backgroundPlasmaProfiles']['FlowVelocity']['flowVz'] = -2000
    i.Input['backgroundPlasmaProfiles']['Bfield']['rString'] = 'br'
    i.Input['backgroundPlasmaProfiles']['Bfield']['zString'] = 'bz'
    i.Input['backgroundPlasmaProfiles']['Bfield']['yString'] = 'bt'
    i.Input['backgroundPlasmaProfiles']['biasPotential'] = -3*25
    i.Input['backgroundPlasmaProfiles']['Efield']['Er'] = 0.02
    i.Input['backgroundPlasmaProfiles']['Efield']['Et'] = 40
    i.Input['backgroundPlasmaProfiles']['Efield']['Ez'] = 2
    i.Input['diagnostics']['nettx0'] = 1.38
    i.Input['diagnostics']['nettx1'] = 1.58
    i.Input['diagnostics']['nX'] = 100

    # Optimization
    i.Input['flags']['GEOM_HASH_SHEATH'] = 0
    i.Input['flags']['GEOM_HASH'] = 0
    i.Input['flags']['GENERATE_LC'] = 0

    # Write input file
    i.WriteInputFile(Folder=folder,OverWrite=True)

if __name__ == '__main__':
    make_input(100)