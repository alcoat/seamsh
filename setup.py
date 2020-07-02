# seamsh - Copyright (C) <2010-2020>
# <Universite catholique de Louvain (UCL), Belgium
# 	
# List of the contributors to the development of seamsh: see AUTHORS file.
# Description and complete License: see LICENSE file.
# 	
# This program (seamsh) is free software: 
# you can redistribute it and/or modify it under the terms of the GNU Lesser General 
# Public License as published by the Free Software Foundation, either version
# 3 of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License
# along with this program (see COPYING file).  If not, 
# see <http://www.gnu.org/licenses/>.

import setuptools 
import sys
import os

#with open("README.md", "r") as fh:
#    long_description = fh.read()

long_description =''

version = "0.0.1"
commit_tag = os.environ.get("CI_COMMIT_TAG")
if commit_tag and (commit_tag.startswith("v-") or commit_tag.startswith("w-")):
    version = commit_tag[2:]




setuptools.setup(
    name="seamsh",
    version=version,
    author="Jonathan Lambrechts",
    author_email="jonathan.lambrechts@uclouvain.be",
    description="Ocean mesh generation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    include_package_data=True,
    url="https://www.migflow.be",
    packages=["seamsh"],
    package_dir={"seamsh":"seamsh"},
    #ext_modules=[CMakeExtension("seamshlib")],
    #ext_modules=[setuptools.Extension("seamshlib",["seamsh.c"])],
    package_data={"seamsh":["*.so","*.dll","*.dll.a","*.dylib","COPYING.txt","AUTHORS.txt","LICENSE.txt"]},
    classifiers=[
        "Environment :: Console",
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: POSIX :: Linux",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: MacOS :: MacOS X",
        "Programming Language :: C",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering"
        ],
    install_requires=["scipy","numpy","gdal","gmsh-dev"],
    python_requires='>=3.6'
)




