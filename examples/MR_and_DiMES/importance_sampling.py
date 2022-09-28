import os
from pyGITR.Run import *
from pyGITR.process_functions import *

import sys
sys.path.insert(0, 'surface_model/')
sys.path.append('.')
import surface_model_functions
from surface_model_functions import *

nP = 1000
dt = 1e-11
nT = 1e4



# Run = Run()
# Run.Verbose = True
# Run.SetReferenceDirectory('.')
# Run.SetSimRootPath('~/GeneralAtomics/simulations/')
# Run.ModifParam('input/gitrInput.cfg','impurityParticleSource.nP',nP)
# Run.Clean()
# Run.LaunchBatch()

# from pyGITR.PostProcess import *
# Post = PostProcess(Run.CurrentSimu)

# Post.PlotArray(PlotEStartEnd ,alpha=0.2)
# Post.PlotArray(SurfaceAngle)

# print('Tot:',[S.Data['SurfaceData']['Data']['Tot'] for S in Post.Simulations])


# Post.GetSurfaceData()

# plt.figure()
# plt.scatter(Post.Simulations[0].Data['ParticleStartData']['Data']['x'],Post.Simulations[0].Data['ParticleStartData']['Data']['y'],Post.Simulations[0].Data['ParticleStartData']['Data']['z'])
# plt.scatter(Post.Simulations[0].Data['ParticleEndData']['Data']['x'],Post.Simulations[0].Data['ParticleEndData']['Data']['y'],Post.Simulations[0].Data['ParticleEndData']['Data']['z'])

# print('Tot:',[S.Data['SurfaceData']['Data']['Tot'] for S in Post.Simulations])



# For importance sampling, we want to find which particles hit the surfaces we care about
# Reading geomtry file
GeomFile = "input/gitrGeom.cfg"
area,surface,Z = getGeom(GeomFile)
Zs, Surfaces = getZsandSurfaces(surface,Z)
Energy, surfaceHit = loadPositions(Surfaces)

surfaceW = []
print("The W surface elements:")
for i in range(len(surface)):
    if surface[i] == 1 and Z[i] == 74:
        print(i)
        surfaceW.append(i)
print(surfaceW)

# SURFACE
SurfaceFile = "output/surface.nc"
surface = netCDF4.Dataset(SurfaceFile)
printInfo(surface)
dims, vars = getSurfaceData(surface)
# plotSurface(dims,vars,["spylCounts","sumParticlesStrike","sumWeightStrike"])
plotSurface(dims,vars,["surfEDist","surfReflDist","surfSputtDist"])
    


# for surface_element in np.unique(surfaceHit):
#         if surface_element in Surfaces:
#             W_elem[surface_element] = 0