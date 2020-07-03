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

from osgeo import ogr, osr
import platform
import numpy as np
from scipy.spatial import cKDTree, Delaunay
import logging
import os
import itertools
import ctypes as c
from enum import Enum
import typing


libdir = os.path.dirname(os.path.realpath(__file__))
if platform.system() == "Windows":
    libpath = os.path.join(libdir, "seamsh.dll")
elif platform.system() == "Darwin":
    libpath = os.path.join(libdir, "libseamsh.dylib")
else:
    libpath = os.path.join(libdir, "libseamsh.so")

lib = c.CDLL(libpath)


def _ensure_valid_points(x, pfrom, pto):
    x = np.asarray(x)[:, :2]
    if not pfrom.IsSame(pto):
        x = osr.CoordinateTransformation(pfrom, pto).TransformPoints(x)
        x = np.asarray(x)[:, :2]
    return x


MeshSizeCallback = typing.Callable[[np.ndarray, osr.SpatialReference],
                                   np.ndarray]


class CurveType(Enum):
    """ Determine how curves control points are interpolated. """
    POLYLINE = 1
    """ linear interpolation """
    STRICTPOLYLINE = 2
    """ linear interpolation with a  mesh point forced on each curve point """
    SPLINE = 3
    """ cubic spline interpolation """
    BSPLINE = 4
    """ cubic b-spline interpolation (smoother than spline) """


def _generate_unique_points(x):
    tree = cKDTree(x)
    unique_id = np.full([x.shape[0]], -1, np.int32)
    bbmin = np.min(x, axis=0)
    bbmax = np.max(x, axis=0)
    eps = np.linalg.norm(bbmax-bbmin)*1e-12
    eps = 1e-12
    cid = 0
    for p0, p1 in tree.query_pairs(eps):
        uid = max(unique_id[p0], unique_id[p1])
        if uid == -1:
            uid = cid
            cid += 1
        unique_id[p0] = uid
        unique_id[p1] = uid
    nbreakpoints = cid
    for i, u in enumerate(unique_id):
        if unique_id[i] == -1:
            unique_id[i] = cid
            cid += 1
    xo = np.ndarray([cid, 2])
    xo[unique_id, :] = x
    return xo, unique_id, nbreakpoints


class _Curve:

    def __init__(self, points, tag, curve_type):
        self.points = points[:, :2]
        self.mesh_size = None
        self.tag = tag
        self.curve_type = curve_type


class Domain:
    """ List the domain boundaries, forced mesh points and
    forced mesh lines. """

    def __init__(self, projection: osr.SpatialReference):
        """
        Args:
            projection: coordinate system of the geometry and mesh
        """

        self._projection = projection
        self._interior_curves = []
        self._interior_points = []
        self._curves = []

    def _build_topology(self):
        curvesiter = itertools.chain(self._curves, self._interior_curves)
        allpoints = np.row_stack(list(itertools.chain(
            (curve.points for curve in curvesiter), self._interior_points)))
        self._points, unique_id, nbreakpoints = _generate_unique_points(
            allpoints)
        # assign points id in curves
        cid = 0
        for curve in itertools.chain(self._curves, self._interior_curves):
            n = curve.points.shape[0]
            curve.pointsid = unique_id[cid:cid+n]
            cid += n
        # assign points id for interior points
        self._interior_points_id = unique_id[cid:]
        # break curves on repeated points

        def split_curves(curves, nbk):
            newcurves = []
            for curve in curves:
                breaks = np.where(curve.pointsid[1:-1] < nbk)[0]
                if len(breaks) == 0:
                    newcurves.append(curve)
                else:
                    breaks = [0] + list(breaks+1) + [len(curve.points)-1]
                    for i, j in zip(breaks[:-1], breaks[1:]):
                        ncurve = _Curve(
                            curve.points[i:j+1], curve.tag, curve.curve_type)
                        ncurve.pointsid = curve.pointsid[i:j+1]
                        newcurves.append(ncurve)
            return newcurves
        self._curves = split_curves(self._curves, nbreakpoints)
        self._interior_curves = split_curves(
            self._interior_curves, nbreakpoints)
        p2curve = np.full([self._points.shape[0], 2, 2], -1, int)
        curve2p = np.ndarray([len(self._curves), 2], int)
        for icurve, curve in enumerate(self._curves):
            curve2p[icurve, 0] = curve.pointsid[0]
            curve2p[icurve, 1] = curve.pointsid[-1]
            p0 = p2curve[curve.pointsid[0]]
            i0 = 0 if p0[0][0] == -1 else 1
            p0[i0][0] = icurve
            p0[i0][1] = 0
            p1 = p2curve[curve.pointsid[-1]]
            i1 = 0 if p1[0][0] == -1 else 1
            p1[i1][0] = icurve
            p1[i1][1] = 1
        touched = np.full((len(self._curves)), False, np.bool)
        self._curveloops = []
        count = 0
        try:
            for i, _ in enumerate(self._curves):
                if touched[i]:
                    continue
                self._curveloops.append([])
                j = 0
                while ((not self._curveloops[-1]) or
                        self._curveloops[-1][0][0] != i):
                    self._curveloops[-1].append((i, j == 0))
                    p = curve2p[i, (j+1) % 2]
                    if p2curve[p][0][0] != i:
                        i, j = p2curve[p][0]
                    else:
                        i, j = p2curve[p][1]
                    touched[i] = True
                    assert(i != -1)
                    count += 1
        except AssertionError:
            raise ValueError("Invalid topology")
        loopbboxarea = np.zeros([len(self._curveloops)])
        for i, loop in enumerate(self._curveloops):
            lpts = np.row_stack([self._curves[j].points for j, o in loop])
            bbox = np.max(lpts, axis=0)-np.min(lpts, axis=0)
            loopbboxarea[i] = bbox[0]*bbox[1]
        self._curveloops.insert(
            0, self._curveloops.pop(np.argmax(loopbboxarea)))

    def _add_geometry(self, geometry, tag, projection, curve_type, interior,
                      onlypoints=False):
        if geometry.GetGeometryCount() != 0:
            for subg in geometry:
                self._add_geometry(subg, tag, projection,
                                   curve_type, interior, onlypoints)
            return
        if onlypoints:
            self.add_interior_points(geometry.GetPoints(), projection)
        else:
            assert(geometry.GetGeometryType() == 2)
            if interior:
                self.add_interior_curve(geometry.GetPoints(), tag,
                                        projection, curve_type)
            else:
                self.add_boundary_curve(geometry.GetPoints(), tag,
                                        projection, curve_type)

    def _add_shapefile(self, filename, physical_name_field,
                       interior, points, curve_type):
        driver = ogr.GetDriverByName('ESRI Shapefile')
        data = driver.Open(filename, 0)
        layer = data.GetLayer()
        layerdef = layer.GetLayerDefn()
        physfield = None
        if physical_name_field is not None:
            for i in range(layerdef.GetFieldCount()):
                field_name = layerdef.GetFieldDefn(i).GetName()
                if field_name == physical_name_field:
                    physfield = i
            if physfield is None:
                raise ValueError("field '"+physical_name_field +
                                 "' not found in shapefile")
        layerproj = layer.GetSpatialRef()
        for i in layer:
            phys = i.GetField(physfield) if physfield else "boundary"
            self._add_geometry(i.geometry(), phys, layerproj, curve_type,
                               interior, points)

    def add_interior_points(self, points: np.ndarray,
                            projection: osr.SpatialReference) -> None:
        """ Add forced interior mesh points

        Args:
            x: the points [n,2]
            projection: the points coordinate system
        """
        points = _ensure_valid_points(points, projection, self._projection)
        self._interior_points.append(points)

    def add_boundary_curve(self, points: np.ndarray, physical_tag: str,
                           projection: osr.SpatialReference,
                           curve_type: CurveType = CurveType.POLYLINE) -> None:
        """ Add a tagged curve to the domain boundary.

        Args:
            points: the curve control points [n,2]
            physical_tag: the curve physical_tag
            projection: the points coordinate system
            curve_type: curve interpolation
        """
        points = _ensure_valid_points(points, projection, self._projection)
        curve = _Curve(points, physical_tag, curve_type)
        self._curves.append(curve)

    def add_interior_curve(self, points: np.ndarray, physical_tag: str,
                           projection: str,
                           curve_type: CurveType = CurveType.POLYLINE) -> None:
        """ Adds a tagged curve inside the domain. The curve
        is not part of the domain boundary and will be meshed on both sides.

        Args:
            points: the curve control points [n,2]
            physical_tag: the curve physical tag
            projection: the points coordinate system
            curve_type: curve interpolation
        """
        points = _ensure_valid_points(points, projection, self._projection)
        curve = _Curve(points, physical_tag, curve_type)
        self._interior_curves.append(curve)

    def add_interior_points_shp(self, filename: str,
                                physical_name_field: str = None) -> None:
        """ Adds all points of a shape file as forced mesh points.

        Args:
            filename: path to a shapefile.
            physical_name_field: name of an attribute string field with the
                curves physical tags
        """
        self._add_shapefile(filename, physical_name_field, None, True, None)

    def add_interior_curves_shp(self, filename: str,
                                physical_name_field: str = None,
                                curve_type: CurveType = CurveType.POLYLINE
                                ) -> None:
        """Adds all lines, polylines and polygons of a shape file as forced
        interior mesh lines

        Args:
            filename: path to a shapefile.
            physical_name_field: name of and attribute string field with the
                 curves physical tags
            curve_type: curves interpolation
        """
        self._add_shapefile(filename, physical_name_field, True, False,
                            curve_type)

    def add_boundary_curves_shp(self, filename: str,
                                physical_name_field: str = None,
                                curve_type: CurveType = CurveType.POLYLINE
                                ) -> None:
        """ Adds all lines, polylines and polygons of a shapefile as domain
        boundaries

        Args:
            filename: path to a shapefile.
            physical_name_field: name of an attribute string field containing
                the curves physical tags
            curve_type: curves interpolation
        """
        self._add_shapefile(filename, physical_name_field, False, False,
                            curve_type)


from .gmsh import _curve_sample


def coarsen_boundaries(domain: Domain, x0: typing.Tuple[float, float],
                       x0projection: osr.SpatialReference,
                       mesh_size: MeshSizeCallback,
                       sampling: float) -> Domain:
    """ Creates a new Domain with the same projection and coarsened
    boundaries.

    For the moment, the physical tags are not transferred to the new Domain,
    this will be implemented in the future.

    Args:
        domain: the domain to coarsen
        x0: the coordinates of one point inside the domain.
        x0projection: the coordinates system of x0. For now, this parameter
            is ignored and the x0 should be in the same coordinate system
            as the domain.
        mesh_size: a function returning the desired mesh element size for given
            coordinates
        sampling: should be (at least) two times smaller than the smallest
            target mesh element size, otherwise the coarsening algorithm will
            fail. This parameter will be removed (and determined automatically)
            in the future.
    """

    sampled = [_curve_sample(curve, sampling) for curve in domain._curves]
    x = np.vstack(list(sampled))
    x, _, _ = _generate_unique_points(x)
    x = np.copy(x[:, :2])
    # avoid cocircular points
    eps = (np.max(x, axis=0, keepdims=True) -
           np.min(x, axis=0, keepdims=True))*1e-8
    x = x + np.random.random(x.shape)*eps
    tag = np.hstack(list(
        np.full(curve.points.shape[0], i, np.int32)
        for i, curve in enumerate(domain._curves)))

    def np2c(a, dtype=np.float64):
        tmp = np.require(a, dtype, "C")
        r = c.c_void_p(tmp.ctypes.data)
        r.tmp = tmp
        return r
    tri = Delaunay(x).simplices
    n_ll = c.c_int()
    n_l = c.c_int()
    n_xo = c.c_int()
    p_xo = c.POINTER(c.c_double)()
    p_l = c.POINTER(c.c_int)()
    p_ll = c.POINTER(c.c_int)()
    ms = mesh_size(x, domain._projection)
    lib.gen_boundaries_from_points(
        c.c_int(x.shape[0]), np2c(x), np2c(tag, np.int32),
        c.c_int(tri.shape[0]), np2c(tri, np.int32),
        c.c_double(x0[0]), c.c_double(x0[1]),
        np2c(ms),
        c.byref(p_xo), c.byref(n_xo),
        c.byref(p_l), c.byref(n_l),
        c.byref(p_ll), c.byref(n_ll))
    xbuf = c.cast(p_xo, c.POINTER(n_xo.value*2*c.c_double)).contents
    xo = np.ctypeslib.frombuffer(xbuf, dtype=np.float64)
    xo = xo.reshape([-1, 2]).copy()
    lib.libcfree(p_xo)
    linesbuf = c.cast(p_l, c.POINTER(n_l.value*3*c.c_int))
    lines = np.ctypeslib.frombuffer(linesbuf.contents, dtype=np.int32)
    lines = lines.reshape([-1, 3]).copy()
    lib.libcfree(p_l)
    # lib.libcfree(p_ll)
    odomain = Domain(domain._projection)
    for i, (lfrom, lto, llast) in enumerate(lines):
        pts = list(range(lfrom, lto+1))+[llast]
        odomain.add_boundary_curve(xo[pts], "boundaries",
                                   domain._projection,
                                   curve_type=CurveType.POLYLINE)
    return odomain


__all__ = ["Domain", "CurveType", "coarsen_boundaries"]
