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

from .geometry import Domain
import numpy as np
from osgeo import osr, gdal
from scipy.spatial import cKDTree
import logging
from typing import List
import itertools
from .gmsh import _curve_sample


def _ensure_valid_points(x, pfrom, pto):
    x = np.asarray(x)[:, :2]
    if not pfrom.IsSame(pto):
        x = osr.CoordinateTransformation(pfrom, pto).TransformPoints(x)
        x = np.asarray(x)
    return x


class Distance:
    """Callable evaluating the distance to a set of discretized curves.

    The curves are discretized as sets of points, then the distance to the
    closest point is computed.
    """

    def __init__(self, domain: Domain, sampling: float,
                 tags: List[str] = None):
        """
        Args:
            domain: a Domain object containing the set of curves
            sampling: the interval between two consecutive sampling points.
            tags: List of physical tags specifying the curve from the domain.
                if None, all curves are taken into account.
        """
        points = []
        for curve in itertools.chain(domain._curves, domain._interior_curves):
            if (tags is None) or (curve.tag in tags):
                points.append(_curve_sample(curve, sampling))
        points = np.vstack(points)
        self._tree = cKDTree(points)
        self._projection = domain._projection

    def __call__(self, x: np.ndarray, projection: osr.SpatialReference
                 ) -> np.ndarray:
        """Computes the distance between each point of x and the curves.

        Args:
            x: the points [n,2]
            projection: the coordinate system of the points, should be
                the same coordinate system as the domain, otherwise no
                conversion is done and an exception is raised.
        Returns:
            The distance expressed in the domain unit. [n]
        """
        if not projection.IsSame(self._projection):
            raise ValueError("incompatible projection")
        x = x[:, :2]
        return self._tree.query(x)[0]


class Raster:
    """Callable to evaluate a raster field loaded from a file."""

    def __init__(self, filename: str):
        """
        Args:
            filename: A geotiff file or any other raster supported by gdal.
        """
        src_ds = gdal.Open(filename)
        self._geo_matrix = src_ds.GetGeoTransform()
        self._data = src_ds.GetRasterBand(1).ReadAsArray()
        assert(self._geo_matrix[2] == 0.)
        assert(self._geo_matrix[4] == 0.)
        self._projection = osr.SpatialReference()
        self._projection.ImportFromWkt(src_ds.GetProjection())

    def __call__(self, x: np.ndarray, projection: osr.SpatialReference
                 ) -> np.ndarray:
        """Evaluates the field value on each point of x.

        Keyword arguments:
            x: the points [n,2]
            projection: the coordinate system of the points
        Returns:
            The field value on points x. [n]
        """
        x = _ensure_valid_points(x, projection, self._projection)
        gm = self._geo_matrix
        xi = (x[:, 0]-gm[3])/gm[5]
        eta = (x[:, 1]-gm[0])/gm[1]
        return self._data[xi.astype(int), eta.astype(int)]
