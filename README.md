# Workflow of pyGITR

## make_geom.py
Define the geometry and mesh of the simulation.

## setup_geom.py
Converts the geometry mesh to gitrGeom.cfg. Many plot-checks, such as inDir, centroids, norms...

## make_particle_source.py
Define the particle source. Currently, this is static. Eventually want to make this a function of impinging flux (such as C flux).

## input.py
Script that makes gitrInput.cfg. Set physics and simulation options here.

## setup.py
Script that runs all of the above

