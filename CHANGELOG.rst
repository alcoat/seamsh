Changes
=======

seamsh 0.4.9
------------
* add transifinite curves (physical lines with prescribed number of elements)

seamsh 0.4.8
------------
* numpy 1.24 compatibility (np.bool -> bool)

seamsh 0.4.7
------------
* write WKT projection as gmsh model attribute
* default msh output version set to 4.1

seamsh 0.4.6
------------
* License files in source package

seamsh 0.4.5
------------
* Inpoly field (requires shapely)

seamsh 0.4.4
------------
* fix unrefine when no identical points
* ugrid conversion (requires xarray)
* fix domain physical name

seamsh 0.4.3
------------
* add output_srs parameter to gmsh.mesh
* add gmsh.reproject function

seamsh 0.4.2
------------
* do not crash on features without geometry in shapefiles

seamsh 0.4.1
------------
* do not call gmsh.initialize() if already initialized
* use a temporary gmsh model for meshing

seamsh 0.4.0
------------

* precompute and smooth mesh size field on a quadtree-like mesh
* build library directly from setup.py, no cmake anymore
* discretize curves in their own projection, not in the model projection
* handle intersection of more than 2 interior curves on the same point
* source distribution on pypi
