Seamsh_ is a python library wrapping gmsh_, gdal_  and scipy_ to simplify the generation of unstructured meshes. It is primarily developed for coastal ocean simulations but can be used in other GIS contexts.

Seamsh allows to :

- import ESRI shapefiles to define tagged domain boundaries, interior lines and interior points
- define arbitrary mesh elements size fields based on distances from lines or raster files
- create a low-resolution valid topology from high-resolution non conformal (i.e. intersecting) data

.. _gmsh : https://www.gmsh.info
.. _gdal : https://gdal.org
.. _scipy : https://www.scipy.org
.. _Seamsh : https://git.immc.ucl.ac.be/jlambrechts/seamsh

Seamesh is distributed under the GPL_. See the gitlab page of the project for the `source code`__ , documentation__ and `bug reports`__
Binary packages for 64 bits linux, windows and OSX are available on pypi_.

__ https://git.immc.ucl.ac.be/jlambrechts/seamsh
__ http://jlambrechts.git-page.immc.ucl.ac.be/seamsh
__ https://git.immc.ucl.ac.be/jlambrechts/seamsh/-/issues
.. _GPL : https://www.gnu.org/licenses/gpl-3.0.html
.. _pypi : https://test.pypi.org/project/seamsh
