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

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _build_ext
from wheel.bdist_wheel import bdist_wheel as _bdist_wheel
import pkg_resources
import sys
import os

version = "0.1"
commit_tag = os.environ.get("CI_COMMIT_TAG")
if commit_tag and (commit_tag.startswith("v-") or commit_tag.startswith("w-")):
    version = commit_tag[2:]

lib = Extension("seamsh.libseamsh", sources = ["seamshlib/seamsh.c","seamshlib/polymesh.c","seamshlib/robustPredicates.c"])

with open("README.rst", "r") as fh:
    long_description = fh.read()
with open("CHANGELOG.rst", "r") as fh:
    long_description += "\n" + fh.read()

class bdist_wheel(_bdist_wheel):
    def get_tag(self):
        otag = _bdist_wheel.get_tag(self)
        return ("py3", "none", otag[2])

class build_ext(_build_ext):
    def get_ext_filename(self, fullname):
        filename = os.path.join(*fullname.split('.')) + ".so"
        return filename

setup(
    name="seamsh",
    version=version,
    author="Jonathan Lambrechts",
    author_email="jonathan.lambrechts@uclouvain.be",
    description="Ocean mesh generation",
    long_description=long_description,
    long_description_content_type = "text/x-rst",
    include_package_data=True,
    url="https://git.immc.ucl.ac.be/jlambrechts/seamsh",
    packages=["seamsh"],
    ext_modules = [lib],
    package_dir={"seamsh":"seamsh"},
    cmdclass = {'bdist_wheel':bdist_wheel,"build_ext":build_ext},
    package_data={"seamsh":["*.so"]},
    license_files=["COPYING.txt","AUTHORS.txt","LICENSE.txt"],
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
    install_requires=["scipy","numpy","gdal","gmsh"],
    python_requires='>=3.6'
)
