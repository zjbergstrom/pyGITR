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
zDiMES=z_shift
rDiMES=0.025

# Define bounding box
xBox = xDiMES
yBox = yDiMES
zBox = zDiMES
dxBox = rDiMES*8.5
dyBox = rDiMES*8.5
dzBox = rDiMES*8.5

# Define metal ring
xRing = 0
yRing = 0
zRing = zDiMES
innerR = 1.4
outerR = 1.45

# Tags
TagDiMES = 10
TagOuterRing = 20
TagInnerRing = 30
TagBox = 40

# Box dimensions
print("xmin",xBox-dxBox/2)
print("xmax",xBox+dxBox/2)
print("ymin",yBox-dyBox/2)
print("ymax",yBox+dyBox/2)
print("zmin",zBox)
print("zmax",zBox+dzBox)


# Create DiMES cap
s1 = gmsh.model.occ.addDisk(xDiMES,yDiMES,zDiMES,rDiMES,rDiMES,TagDiMES) # the second r is used for elipse definition, (2,0)

# Disks
s2 = gmsh.model.occ.addDisk(xRing,yRing,zRing,outerR,outerR,TagOuterRing,)
s3 = gmsh.model.occ.addDisk(xRing,yRing,zRing,innerR,innerR,TagInnerRing)

# Bounding Box
s4 = gmsh.model.occ.addBox(xBox-dxBox/2, yBox-dyBox/2, zBox, dxBox, dyBox, dzBox, TagBox) # (3,40), (2,2), (2,3), (2,4), (2,5), (2,6), (2,7)
gmsh.model.occ.remove([(3, TagBox)]) # removes the 3d shape from the 3d shape list, leaving only 2d shapes in the 2d shape list

e1 = gmsh.model.occ.getEntities(2) # get 2D entities
print(e1)

# Create Annulus
annulus = gmsh.model.occ.cut([(2, s2)], [(2, s3)], -1, removeTool=False) # cuts out the MR and DiMES cap from the bottom surface
print(annulus)

e2 = gmsh.model.occ.getEntities(2) # get 2D entities
print(e2)

# Remove Annulus and DiMES from Box
box = gmsh.model.occ.cut([(2, 35)], [(2, s1), (2, s2), (2,s3)], -1, removeTool=False) # cuts out the MR and DiMES cap from the bottom surface
print(box)



# Synchronize necessary before mesh setup and generation
gmsh.model.occ.synchronize()

# Set number of elements on the boundary of each dots and DiMES cap
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 1)
gmsh.option.setNumber("Mesh.MinimumElementsPerTwoPi", 20)

# Define groups to allow setting properties of elements when generating geometry input with GeomSetup
# SetGroups(gmsh.model, TagDiMES, "DiMES", [255, 0, 0])
# SetGroups(gmsh.model, TagBoundBox, "BoundBox", [0, 255, 0])
# SetGroups(gmsh.model, TagMR, "Metal Ring", [0,0,255])

# Generate 2D mesh
mesh = gmsh.model.mesh.generate(2)

# Launch the GUI to see the results:
gmsh.fltk.run()

# Write mesh into a meshio format
# gmsh.write("metal_rings_and_DiMES.msh")

# Close gmsh session
gmsh.finalize()











# Create 2D surfaces
# New metal ring
# gmsh.model.occ.addDisk(0,0,0,outerR,outerR,TagOuterRing)
# gmsh.model.occ.addDisk(0,0,0,innerR,innerR,TagInnerRing)
# gmsh.model.occ.cut([(2, TagOuterRing)], [(2, TagInnerRing)], -1, removeTool=False) # cuts out the MR and DiMES cap from the bottom surface
# TagMR = gmsh.model.occ.get_entities(2)
# print(TagMR)




# Create metal ring
# r1 = gmsh.model.occ.addRectangle(xBlock, yBlock-dyBlock/2, zBlock, dxBlock, dyBlock) # (2,1)
# TagMR = gmsh.model.occ.get_entities(2)
# print(TagMR)

# # Create DiMES cap
# r4 = gmsh.model.occ.addDisk(xDiMES,yDiMES,zDiMES,rDiMES,rDiMES,TagDiMES0) # the second r is used for elipse definition, (2,0)
# s0 = gmsh.model.occ.getEntities(2) # get 2D entities

# TagDiMES = list(set(s0) - set(TagMR)) # intersects the 2 sets
# print(TagDiMES)

# # Create simulation bounding box
# s0 = gmsh.model.occ.getEntities(2) # get 2D entities
# gmsh.model.occ.addBox(xBox-dxBox/2, yBox-dyBox/2, zBox, dxBox, dyBox, dzBox, TagBox0) # (3,40), (2,2), (2,3), (2,4), (2,5), (2,6), (2,7)
# gmsh.model.occ.remove([(3, TagBox0)]) # removes the 3d shape from the 3d shape list, leaving only 2d shapes in the 2d shape list
# gmsh.model.occ.cut([(2, 6)], [(2, TagDiMES0), (2,r1)], -1, removeTool=False) # cuts out the MR and DiMES cap from the bottom surface
# # gmsh.model.occ.cut([(2, 6)], [(2, TagDiMES0), (2,TagMR)], -1, removeTool=False) # cuts out the MR and DiMES cap from the bottom surface

# s1 = gmsh.model.occ.getEntities(2)
# TagBoundBox = list(set(s1) - set(s0)) # intersects the 2 sets
# print(TagBoundBox)

