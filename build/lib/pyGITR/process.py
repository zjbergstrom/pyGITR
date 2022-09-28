from pyGITR.process_functions import *

# sys.path.insert(0, '../particleSource/')
# sys.path.append('.')
# from particleSource_functions import *


# HISTORY
# HistoryFile = "gitrm-history.nc"
# history = netCDF4.Dataset(HistoryFile)
# printInfo(history)
# dims, vars = getHistoryData(history)
# plotHistory2D(dims,vars,10)
# plotHistory3D(dims,vars,20)
# hits = findHits(dims,vars,1,plotLongTracks=False)


# SURFACE
# SurfaceFile = "gitrm-surface.nc"
# surface = netCDF4.Dataset(SurfaceFile)
# printInfo(surface)
# dims, vars = getSurfaceData(surface)
# plotSurface(dims,vars,["spylCounts","sumParticlesStrike","sumWeightStrike"])


# SPECTROSCOPY
# SpecFile = "gitrm-spec.nc"
# spec = netCDF4.Dataset(SpecFile)
# printInfo(spec)
# dims, vars = getSpecData(spec)
# plotSpec(dims,vars)

# PLASMA
# PlasmaFile = "bField.nc"
# plasma = netCDF4.Dataset(PlasmaFile)
# printInfo(plasma)
# dims, vars = getPlasmaData(plasma)
# plotPlasma(dims,vars)

# POSITIONS
# PositionsFile = "positions.nc"
# positions = netCDF4.Dataset(PositionsFile)
# printInfo(positions)
# dims, vars = getPositionsData(positions)
# plotPositions(dims,vars)


# PARTICLESOURCE
# ParticlesFile = "particleSourceForGITRm.nc"
# particles = netCDF4.Dataset(ParticlesFile)
# # printInfo(particles)
# dims, vars = getParticlesData(particles)


# x1,x2,x3,y1,y2,y3,z1,z2,z3,a,b,c,d,area,plane_norm,surface,indir = loadCFG()

# # The bad particle is indexed at 405
# data = {}
# for i in range(len(surface)):
#     if surface[i]==1:
#         data[i] = 1
# print(len(data))

# counter = 0 # the bad particle is at counter 406
# for num,i in enumerate(data.keys()):
#     if num%10==0:
#         counter+=1
#         if counter==406:

#             # plot corresponding element
#             surface_x, surface_y, surface_z = [x1[i],x2[i],x3[i],x1[i]], \
#                     [y1[i],y2[i],y3[i],y1[i]], [z1[i],z2[i],z3[i],z1[i]]

#             # plot corresponding particle
#             x,y,z = vars['x'][counter-1], vars['y'][counter-1], vars['z'][counter-1]

#             plotPointsAndElement(surface_x,surface_y,surface_z,x,y,z,a[i],b[i],c[i],indir[i])