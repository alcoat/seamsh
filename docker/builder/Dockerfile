#seamsh - Copyright (C) <2010-2020>
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

FROM ubuntu:20.04

env DEBIAN_FRONTEND=noninteractive 
#linux
RUN apt update && apt install -y git python3 make cmake gcc
#windows
RUN apt update && apt install -y mingw-w64 mingw-w64-tools
run echo "\
SET(CMAKE_SYSTEM_NAME Windows)\n\
SET(CMAKE_C_COMPILER x86_64-w64-mingw32-gcc)\n\
SET(CMAKE_CXX_COMPILER x86_64-w64-mingw32-g++)\n\
SET(CMAKE_RC_COMPILER x86_64-w64-mingw32-windres)\n\
SET(CMAKE_Fortran_COMPILER x86_64-w64-mingw32-gfortran)\n\
SET(CMAKE_FIND_ROOT_PATH /mingw64)\n\
SET(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)\n\
SET(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)\n\
SET(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)\n"\
> /cmake-mingw

RUN apt update && apt install -y git python3-setuptools python3-wheel twine


### osxcross

RUN apt-get update && apt-get install -y clang patch libssl-dev liblzma-dev libxml2-dev llvm-dev uuid-dev zlib1g-dev

RUN useradd -ms /bin/bash validator
RUN mkdir -p /osxcross && chown validator /osxcross 
USER validator
RUN git clone https://github.com/tpoechtrager/osxcross.git /osxcross
COPY MacOSX10.15.sdk.tar.xz /osxcross/tarballs/MacOSX10.15.sdk.tar.xz
RUN cd /osxcross && SDK_VERSION=10.15 UNATTENDED=1 ./build.sh
ENV PATH=/osxcross/target/bin/:$PATH
ENV OSXCROSS_HOST=x86_64-apple-darwin19


### doc
USER root
RUN apt-get update && apt-get install -y python3-pip gmsh gdal-bin python3-scipy && pip3 install gmsh-dev sphinx sphinx-rtd-theme sphinx-autodoc-typehints gdal
ENV PYTHONPATH=/usr/local/lib/python3/dist-packages/gmsh-git-Linux64-sdk/lib/

USER validator

VOLUME ["/etc/gitlab-runner"]

WORKDIR /home/validator
