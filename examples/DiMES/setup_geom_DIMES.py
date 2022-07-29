# -*- coding: utf-8 -*-

import os
from pyGITR.Geom import GeomSetup

# Load geometry
g = GeomSetup('small_large_dots_DiMES.msh', Verbose=True)

# Show existing groups
#g.ShowGroups()

# Plot mesh for some groups
#g.Plot(["DiMES"])
#g.Plot(["Small Dot","Large Dot"])
g.Plot(["BoundBox"])

# Zoom on the plot (the zoom function for 3D axis in matplotlib 3.3.3 is not working so it has to be done manually)
g.SetAxisLim(-0.08, 0.08)

# Show Centroids
#g.ShowCentroids()

# Show normals for the DiMES
# Must initiallize before with Plot
#g.ShowNormals("DiMES")
#g.ShowNormals(["Small Dot","Large Dot"])
#g.ShowNormals(["BoundBox"])

# Set properties for each group
# set 'surface' property. Default is 0.
g.SetElemAttr([], 'surface', 1)
g.SetElemAttr(["DiMES","Large Dot","Small Dot"], 'surface', 1)

# set inDir. Empty list = all elements
g.SetElemAttr([],'inDir',1)
g.SetElemAttr(["BoundBox"], 'inDir',-1) # changed this from +1
g.ShowInDir(["BoundBox"])

# Set Z for material
g.SetElemAttr(["DiMES"], 'Z', 6)
g.SetElemAttr(["Large Dot","Small Dot"], 'Z', 74)
g.SetElemAttr(["BoundBox"], 'Z', 1)

# Set potential for biasing. Empty list = all elements
g.SetElemAttr([],'potential',0)

# Plot geometry showing values of Z with color
#g.Plot(ElemAttr='Z', Alpha=0.1)

# Write the geometry file
g.WriteGeomFile()

os.system("mv gitrGeom.cfg input")



