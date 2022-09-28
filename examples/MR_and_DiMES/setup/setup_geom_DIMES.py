# -*- coding: utf-8 -*-

import os
import matplotlib.pyplot as plt
from pyGITR.Geom import GeomSetup

# Load geometry
g = GeomSetup('metal_rings_and_DiMES.msh', Verbose=True)

# Show existing groups
# g.ShowGroups()

# Plot mesh for some groups
# g.Plot(["DiMES","Metal Ring"])
# g.Plot(["BoundBox"])
g.Plot([])

# g.SetAxisLim(-1.25, 1.25)
border = 0.01
g.SetAxisLim3D((1.37875-border, 1.59125+border), (-0.10625-border, 0.10625+border), (-1.25-border,-1.0375+border))
# plt.show()

# Show Centroids
# g.ShowCentroids()
# plt.show()

# Show Normals
# g.ShowNormals(["Metal Ring","DiMES"])
# g.ShowNormals(["BoundBox"])
# plt.show()

# Set properties for each group
g.SetElemAttr([], 'surface', 1)
g.SetElemAttr(["BoundBox"], 'surface', 0)

# set inDir. Empty list = all elements
g.SetElemAttr([],'inDir',1)
g.SetElemAttr(["BoundBox"], 'inDir',-1) # changed this from +1

g.ShowInDir(["BoundBox"])
# g.ShowInDir(["DiMES","Metal Ring"])
plt.show()

# Set Z for material
g.SetElemAttr(["DiMES"], 'Z', 6)
g.SetElemAttr(["Metal Ring"], 'Z', 74)
g.SetElemAttr(["BoundBox"], 'Z', 0)

# Set potential for biasing. Empty list = all elements
g.SetElemAttr([],'potential',0)

# Plot geometry showing values of Z with color
# g.Plot(ElemAttr='Z', Alpha=0.1)
# border = 0.01
# g.SetAxisLim3D((1.37875-border, 1.59125+border), (-0.10625-border, 0.10625+border), (-1.25-border,-1.0375+border))
# plt.show()

# Write the geometry file
g.WriteGeomFile(Folder="../input",OverWrite=True)




