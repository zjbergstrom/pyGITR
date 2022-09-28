#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
# Geometry and meshing for small/large dots experiments on DiMES
# Reference paper: https://iopscience.iop.org/article/10.1088/1361-6587/ab5144/meta
@author: Jerome Guterl (guterlj@fusion.gat.com)
"""

from tokenize import Pointfloat
from uuid import getnode
import gmsh
from pyGITR.gmsh_helper import SetGroups
from pyGITR.Geom import GeomSetup

# Define DiMES cap surface
xDiMES=0.0
yDiMES=0.0
zDiMES=0.0
rDiMES=0.025

# Define small dot on DiMES cap surface
xsDot=-0.01
ysDot=0.0
zsDot=0.0
rsDot=0.0005

# Define large dot on DiMES cap surface
xlDot=0.0
ylDot=0.0
zlDot=0.0
rlDot=0.005

# Define bounding box
xBox = xDiMES
yBox = yDiMES
zBox = zDiMES
dxBox = rDiMES*2.5
dyBox = rDiMES*2.5
dzBox = rDiMES*2.5

# Tags
TagDiMES0=0
TagDiMES=10
TaglDot=20
TagsDot=30
TagBox0=40



# elements = ["DiMES","Small Dot","Large Dot", "BoundaryBox"]
elements = ["BoundaryBox"]
#elements = ["Large Dot"]
for element in elements:

    # Initialize gmsh session
    gmsh.initialize()

    print(element)

    if element == "DiMES":
        # Create 2D surfaces
        gmsh.model.occ.addDisk(xDiMES,yDiMES,zDiMES,rDiMES,rDiMES,TagDiMES0) # the second r is used for elipse definition

        # Create dots
        gmsh.model.occ.addDisk(xsDot,ysDot,zsDot,rsDot,rsDot,TagsDot)
        gmsh.model.occ.addDisk(xlDot,ylDot,zlDot,rlDot,rlDot,TaglDot)
        DiMESTag = gmsh.model.occ.cut( [(2,TagDiMES0)],[(2,TagsDot),(2,TaglDot)], TagDiMES,removeTool=False) # cuts the end from the first in the list, assigns new tag

    if element == "Large Dot":
        gmsh.model.occ.addDisk(xlDot,ylDot,zlDot,rlDot,rlDot,TaglDot)

    if element == "BoundaryBox":
        # Create simulation bounding box
        s0 = gmsh.model.occ.getEntities(2)
        TagBox0 = gmsh.model.occ.addBox(xBox-dxBox/2, yBox-dyBox/2, zBox, dxBox, dyBox, dzBox, TagBox0)
        gmsh.model.occ.remove([(3, TagBox0)])
        #gmsh.model.occ.cut([(2, 6)], [(2, TagDiMES0)], -1, removeTool=False)
        #s1 = gmsh.model.occ.getEntities(2)
        #TagBoundBox = list(set(s1) - set(s0))

    # Synchronize necessary before mesh setup and generation
    gmsh.model.occ.synchronize()

    mesh = gmsh.model.mesh.generate(2)

    # Launch the GUI to see the results:
    gmsh.fltk.run()

    # Write mesh into a meshio format
    gmsh.write("{}.msh".format(element))

    g = GeomSetup('{}.msh'.format(element), Verbose=True)
    g.SetElemAttr([], 'surface', 0)
    g.SetElemAttr([], 'inDir',1)
    g.SetElemAttr([], 'Z', 6)
    g.SetElemAttr([],'potential',0)
    g.WriteGeomFile(FileName='{}.cfg'.format(element))

# Close gmsh session
gmsh.finalize()