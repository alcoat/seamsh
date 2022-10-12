# seamsh - Copyright (C) <2010-2020>
# <Universite catholique de Louvain (UCL), Belgium
#
# List of the contributors to the development of seamsh: see AUTHORS file.
# Description and complete License: see LICENSE file.
#
# This program (seamsh) is free software:
# you can redistribute it and/or modify it under the terms of the GNU
# Lesser General Public License as published by the Free Software Foundation,
# either version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program (see COPYING file).  If not,
# see <http://www.gnu.org/licenses/>.

import sys
import uuid
import struct
import atexit
import os
import time
from itertools import chain
import numpy as np
from typing import List, Tuple, Callable
from osgeo import gdal, ogr, osr
from scipy.spatial import cKDTree, Delaunay
from enum import Enum

import platform
import ctypes as c

try :
    import xarray as xr
    xarray_available = True
except:
    xarray_available = False

try :
    import shapely.geometry
    shapely_available = True
except:
    shapely_available = False


def log(txt, title=False):
    decorator = " * " if title else ""
    print(decorator+txt+decorator)


class ProgressLog:

    def __init__(self, msg, title=False):
        self._tic = time.time()
        self._last_msg = ""
        log(msg, title)
        if sys.stdout.isatty():
            print("")

    def log(self, msg):
        if sys.stdout.isatty():
            print("\033[F\033[K", end="")  # Cursor up one line and clear
            print(msg + " (%.1fs)" % (time.time()-self._tic))
        self._last_msg = msg

    def end(self):
        if not sys.stdout.isatty():
            print(self._last_msg)


_transform_cache = {}


def project_points(x, pfrom, pto):
    x = np.asarray(x)[:, :2]
    if not pfrom.IsSame(pto):
        trans = _transform_cache.get((pfrom, pto), None)
        if trans is None:
            trans = osr.CoordinateTransformation(pfrom, pto)
            _transform_cache[(pfrom, pto)] = trans
        x = trans.TransformPoints(x)
        x = np.asarray(x)[:, :2]
    return x


libdir = os.path.dirname(os.path.realpath(__file__))
libpath = os.path.join(libdir, "libseamsh.so")

lib = c.CDLL(libpath)


def np2c(a, dtype=np.float64):
    tmp = np.require(a, dtype, "C")
    r = c.c_void_p(tmp.ctypes.data)
    r.tmp = tmp
    return r


class PolyMesh:
    def __init__(self, bbmin, bbmax):
        lib.polymesh_new.restype = c.c_void_p
        self._ptr = c.c_void_p(lib.polymesh_new(np2c(bbmin), np2c(bbmax)))
        self._n = 0

    def add_points(self, v):
        tags = np.arange(self._n, self._n+v.shape[0], dtype=np.int32)
        self._n += v.shape[0]
        lib.polymesh_add_points(
                self._ptr,
                c.c_int(v.shape[0]),
                np2c(v), np2c(tags, np.int32))

    def get_triangles(self):
        lib.polymesh_n_faces.restype = c.c_int
        n = lib.polymesh_n_faces(self._ptr)
        tri = np.empty((n, 3), np.int32)
        lib.polymesh_faces(self._ptr, np2c(tri, np.int32))
        return tri
