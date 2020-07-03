Installation
============
Dependencies
------------
Seamsh depends on the following python3 packages:

- `gmsh-dev <https://pypi.org/project/gmsh-dev>`_
- `gdal <https://pypi.org/project/gdal>`_
- `scipy <https://pypi.org/project/scipy>`_

Depending on your platform, gdal and scipy can probably be installed directly from your package manager, otherwise pip can be used (follow the links above for installation instructions).

To install gmsh-dev, follow the instuctions on the gmsh-dev_ pypi page. Do not forget, to uninstall other previous gmsh installations. Check that gmsh-dev is correctly installed with the following command.

.. code-block:: bash

   python3 -c "import gmsh;gmsh.initialize();print(gmsh.option.getString('General.Version'))"

The output should be a version number >= 3.6.1.
On some platforms, it is necessary to set the python path manually after installation e.g. :

.. code-block:: bash

   export "PYTONPATH=$HOME/.local/lib/python3.8/site-packages/gmsh-git-Linux64-sdk/lib/:$PYTHONPATH"


Installation with pip
---------------------
Pre-compiled seamsh libraries are provided for 64bit linux, windows and OSX.

.. code-block:: bash 

   pip3 install seammsh

Installation from sources
-------------------------
To build the package from sources, it is necessary to first compile the (small) seamsh C library. For example on linux/OSX :

.. code-block:: bash 
  
   git clone https://git.immc.ucl.ac.be/jlambrechts/seamsh
   cd seamsh/seamshlib
   mkdir build
   cd build
   cmake .. -DCMAKE_BUILD_TYPE=Release
   make
   cp *.so ../../seamsh
   cd ../../
   python3 setup.py build
