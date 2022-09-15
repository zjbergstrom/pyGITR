#from pyGITR import PostProcess
from pyGITR.PostProcess import *
from pyGITR.Run import *

'''
Information about how to use PostProcess.py can be found in the following:
./flat_surface/read_output.py:Data = PostProcess(Run.CurrentSimu)
./flat_surface/input.py:from pyGITR.PostProcess import *
./flat_surface/input.py:Post = PostProcess(Run.CurrentSimu)
./micro-trenches/input_single_microtrench.py:from pyGITR.PostProcess import *
./micro-trenches/input_single_microtrench.py:Post = PostProcess(Run.CurrentSimu)
'''

Post = PostProcess(Run.CurrentSimu)
print(Post)
#pp.GetHistoryData()
#pp.GetParticleData()
#pp.GetSurfaceData()
