import os
from pyGITR.Run import *


os.system("python metalrings_DIMES.py")
os.system("python setup_geom_DIMES.py")
os.system("python particles_distribution_generator.py")
os.system("python input.py")
# os.system("./GITR_no_geom_hash")


''' I've commented below because I just want to do the above, the below is good though
# ############################
Run = Run()
Run.Verbose = True

Ti=25

Run.SetReferenceDirectory('.')
Run.SetSimRootPath('~/GeneralAtomics/simulations/')
#Run.AddParamScan('input/gitrGeom.cfg',{'geom.lambda':[ 0.00005, 0.0001, 0.00025, 0.0005, 0.00075, 0.001, 0.002, 0.005]})
#Run.AddParamScan('input/gitrGeom.cfg',{'geom.lambda':[ 0.0001, 0.005]})
Run.AddParamScan('input/gitrInput.cfg',{'backgroundPlasmaProfiles.biasPotential':[-3*Ti, 0]})
#Run.ModifParam('input/gitrInput.cfg','impurityParticleSource.nP',1000)

Run.SetupScan(OverWrite=True)

Run.Clean()
Run.LaunchBatch()

# #%%
# from pyGITR.PostProcess import *
# Post = PostProcess(Run.CurrentSimu)
# Post.PlotArray(PlotEStartEnd ,alpha=0.2)
# Post.PlotArray(SurfaceAngle)

# print('Tot:',[S.Data['SurfaceData']['Data']['Tot'] for S in Post.Simulations])







# Data = PostProcess(Run.CurrentSimu)
# Data.GetSurfaceData()









# ################################
# Run = pyGITR.Run()
# Run.Verbose = True

# Run.SetReferenceDirectory('.')
# Run.SetSimRootPath('~/simulations/{}_{}_{}__{}_single_microtrench/'.format(Elem,charge, Ti,nP))
# Run.AddParamScan('input/gitrGeom.cfg',{'geom.lambda':[0.003]})
# #Run.AddParamScan('input/gitrGeom.cfg',{'geom.lambda':[ 0.0001, 0.005]})
# #Run.AddParamScan('input/gitrInput.cfg',{'backgroundPlasmaProfiles.biasPotential':[-3*Ti, 0]})
# Run.ModifParam('input/gitrInput.cfg','impurityParticleSource.nP',nP)

# Run.SetupScan(OverWrite=True)


# Run.Clean()
# Run.LaunchBatch()

# #%%
# from pyGITR.PostProcess import *
# Post = PostProcess(Run.CurrentSimu)
# Post.PlotArray(PlotEStartEnd ,alpha=0.2)
# Post.PlotArray(SurfaceAngle)
# plt.figure()
# plt.scatter(Post.Simulations[0].Data['ParticleStartData']['Data']['x'],Post.Simulations[0].Data['ParticleStartData']['Data']['y'],Post.Simulations[0].Data['ParticleStartData']['Data']['z'])
# plt.scatter(Post.Simulations[0].Data['ParticleEndData']['Data']['x'],Post.Simulations[0].Data['ParticleEndData']['Data']['y'],Post.Simulations[0].Data['ParticleEndData']['Data']['z'])

# print('Tot:',[S.Data['SurfaceData']['Data']['Tot'] for S in Post.Simulations])
'''