Seamsh_ is a Python library wrapping gmsh_, gdal_  and scipy_ to simplify the generation of unstructured meshes. While primarily developed for coastal ocean simulations, it can be used in other GIS contexts.

Main Features :

- Import ESRI shapefiles to define tagged domain boundaries, interior lines and interior points.
- Define arbitrary mesh elements size fields based on distances from lines or raster files.
- Create a low-resolution valid topology from high-resolution non-conformal (i.e. intersecting) data.

Seamesh is distributed under the GPL_. See the project gitlab page for the `source code`_  and `bug reports`_. The documentation_ contains examples, Python API reference and installation instructions.
Binary packages for 64 bits linux, windows and OSX are available on pypi_.

.. _gmsh : https://www.gmsh.info
.. _gdal : https://gdal.org
.. _scipy : https://www.scipy.org
.. _Seamsh : https://git.immc.ucl.ac.be/jlambrechts/seamsh
.. _source code : https://git.immc.ucl.ac.be/jlambrechts/seamsh
.. _documentation : http://jlambrechts.git-page.immc.ucl.ac.be/seamsh
.. _bug reports : https://git.immc.ucl.ac.be/jlambrechts/seamsh/-/issues
.. _GPL : https://www.gnu.org/licenses/gpl-3.0.html
.. _pypi : https://test.pypi.org/project/seamsh
