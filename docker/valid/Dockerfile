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

FROM ubuntu:22.04

RUN apt update
ENV DEBIAN_FRONTEND=noninteractive
RUN apt update &&  apt install -y git python3 python3-numpy python3-scipy python3-petsc4py python3-pip gdal-bin python3-shapely

VOLUME ["/etc/gitlab-runner"]

RUN useradd -ms /bin/bash validator
USER validator
RUN pip3 install numpy
RUN pip3 install -i https://gmsh.info/python-packages-dev-nox --force-reinstall --no-cache-dir gmsh
RUN pip3 install gdal
RUN pip3 install sphinx sphinx-rtd-theme sphinx-autodoc-typehints

WORKDIR /home/validator
