#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
# Geometry and meshing for small/large dots experiments on DiMES
# Reference paper: https://iopscience.iop.org/article/10.1088/1361-6587/ab5144/meta
@author: Jerome Guterl (guterlj@fusion.gat.com)
"""

import gmsh
from pyGITR.gmsh_helper import SetGroups

# Initialize gmsh session
gmsh.initialize()

# shift entire mesh to align with center of DiMES
r_shift = 1.485
z_shift = -1.250

# Define DiMES cap surface
xDiMES=0.0+r_shift
yDiMES=0.0
# zDiMES=0.0
zDiMES=z_shift
rDiMES=0.025

# Define small dot on DiMES cap surface
xsDot=xDiMES
ysDot=0.01
zsDot=zDiMES
rsDot=0.0005

# Define large dot on DiMES cap surface
xlDot=xDiMES
ylDot=0.0
zlDot=zDiMES
rlDot=0.005

# Define bounding box
xBox = xDiMES
yBox = yDiMES
zBox = zDiMES
dxBox = rDiMES*8.5
dyBox = rDiMES*8.5
dzBox = rDiMES*8.5

# Define metal rings
dxBlock = 0.05
dyBlock = dyBox
xBlock=xDiMES-0.01-rDiMES-dxBlock
yBlock=yDiMES
zBlock=zDiMES

# Tags
TagDiMES0=0
TagDiMES=10
TaglDot=20
TagsDot=30
TagBox0=40
TagMR=50

# Box dimensions
print("xmin",xBox-dxBox/2)
print("xmax",xBox+dxBox/2)
print("ymin",yBox-dyBox/2)
print("ymax",yBox+dyBox/2)
print("zmin",zBox)
print("zmax",zBox+dzBox)

# MR dimensions
# print('xmin',xBlock-dxBlock/2)
# print('xmax',dxBlock)
# print('ymin',yBlock)
# print('ymax',dyBlock)

# When add a shape, return an int tag, either automatically, or the tag given
# Get entities returns the dimensionality and tag of the current entities
# Cut removes a list of specified entities (via their tags), and returns ???

# there are 3d entities, and 2d entities
# when the 3d box is made, 1 3d entitiy with the assigned tag is made
# along with the 2d entities that comprise the 3d shape

# Create metal ring
r1 = gmsh.model.occ.addRectangle(xBlock, yBlock-dyBlock/2, zBlock, dxBlock, dyBlock)
print(r1)
# Create 2D surfaces
gmsh.model.occ.addDisk(xDiMES,yDiMES,zDiMES,rDiMES,rDiMES,TagDiMES0) # the second r is used for elipse definition

# Create simulation bounding box
s0 = gmsh.model.occ.getEntities(2)
gmsh.model.occ.addBox(xBox-dxBox/2, yBox-dyBox/2, zBox, dxBox, dyBox, dzBox, TagBox0)
gmsh.model.occ.remove([(3, TagBox0)]) # removes the 3d shape from the 3d shape list, leaving only 2d shapes in the 2d shape list
gmsh.model.occ.cut([(2, 6)], [(2, TagDiMES0),(2,r1)], -1, removeTool=False)
#gmsh.model.occ.cut([(2, 6)], [(2, TagDiMES0)], -1, removeTool=False)
s1 = gmsh.model.occ.getEntities(2)
TagBoundBox = list(set(s1) - set(s0))
print(type(s1))
print('s1',s1)
print(TagBoundBox)

# Create dots
gmsh.model.occ.addDisk(xsDot,ysDot,zsDot,rsDot,rsDot,TagsDot)
gmsh.model.occ.addDisk(xlDot,ylDot,zlDot,rlDot,rlDot,TaglDot)
DiMESTag = gmsh.model.occ.cut( [(2,TagDiMES0)],[(2,TagsDot),(2,TaglDot)], TagDiMES,removeTool=False) # cuts the end from the first in the list, assigns new tag
print('DiMESTag',DiMESTag)
# Synchronize necessary before mesh setup and generation
gmsh.model.occ.synchronize()

# Set number of elements on the boundary of each dots and DiMES cap
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 1)
gmsh.option.setNumber("Mesh.MinimumElementsPerTwoPi", 50)

# Prevent very small elements in small dots
gmsh.option.setNumber("Mesh.MeshSizeMin", rsDot/3)

# Define groups to allow setting properties of elements when generating geometry input with GeomSetup
TaglDot = [(2,20)]
TagsDot = [(2,30)]
TagMR = [(2,r1)]
SetGroups(gmsh.model, DiMESTag, "DiMES", [255, 0, 0])
SetGroups(gmsh.model, TaglDot, "Large Dot", [0, 255, 0])
SetGroups(gmsh.model, TagsDot, "Small Dot", [0, 0, 255])
SetGroups(gmsh.model, TagBoundBox, "BoundBox", [230, 230, 250])
SetGroups(gmsh.model, TagMR, "Metal Ring", [120,120,120])

# Generate 2D mesh
mesh = gmsh.model.mesh.generate(2)

# Launch the GUI to see the results:
gmsh.fltk.run()

# Write mesh into a meshio format
gmsh.write("metal_rings_and_DiMES.msh")

# Close gmsh session
gmsh.finalize()
