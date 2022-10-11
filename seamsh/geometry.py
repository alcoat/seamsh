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

from . import _tools

__all__ = ["Domain", "CurveType", "coarsen_boundaries"]

MeshSizeCallback = _tools.Callable[[_tools.np.ndarray,
                                    _tools.osr.SpatialReference],
                                   _tools.np.ndarray]


class CurveType(_tools.Enum):
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
    tree = _tools.cKDTree(x)
    unique_id = _tools.np.full([x.shape[0]], -1, _tools.np.int32)
    bbmin = _tools.np.min(x, axis=0)
    bbmax = _tools.np.max(x, axis=0)
    eps = _tools.np.linalg.norm(bbmax-bbmin)*1e-12
    eps = 1e-12
    cid = 0
    pairs = _tools.np.array(list(tree.query_pairs(eps)))
    if pairs.shape[0] == 0 : pairs = pairs.reshape(-1,2)
    #sorting the pairs is necessary to handle points connecting more than 2 lines
    pairs.sort(axis=1)
    pairsint64 = pairs[:,0].astype(_tools.np.int64)*2**32+pairs[:,1].astype(_tools.np.int64)
    pairs = pairs[pairsint64.argsort(),:]
    for p0, p1 in pairs:
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
    xo = _tools.np.ndarray([cid, 2])
    xo[unique_id, :] = x
    breakpoints = _tools.np.zeros([xo.shape[0]], dtype=bool)
    breakpoints[:nbreakpoints] = True
    return xo, unique_id, breakpoints


class _Curve:

    def __init__(self, points, tag, curve_type, projection):
        self.points = _tools.np.array(points)[:, :2]
        self.mesh_size = None
        self.tag = tag
        self.curve_type = curve_type
        self.projection = projection


class _Point:

    def __init__(self, x, tag):
        self.tag = tag
        self.x = x

class Domain:
    """ List the domain boundaries, forced mesh points and
    forced mesh lines. """

    def __init__(self, projection: _tools.osr.SpatialReference):
        """
        Args:
            projection: coordinate system of the geometry and mesh
        """

        self._projection = projection
        self._interior_curves = []
        self._interior_points = []
        self._curves = []

    def _build_topology(self):
        _tools.log("Build topology")
        curvesiter = _tools.chain(self._curves, self._interior_curves)
        allpoints = _tools.np.row_stack(list(_tools.chain(
            (_tools.project_points(curve.points, curve.projection,
                                   self._projection)
             for curve in curvesiter),
            (point.x for point in self._interior_points))))
        self._points, unique_id, breakpts = _generate_unique_points(
            allpoints)
        # assign points id in curves
        cid = 0
        touched_interior = _tools.np.zeros((self._points.shape[0],), bool)
        for curve in _tools.chain(self._curves, self._interior_curves):
            n = curve.points.shape[0]
            curve.pointsid = unique_id[cid:cid+n]
            touched_interior[curve.pointsid] = True
            cid += n
        # assign points id for interior points
        for point in self._interior_points:
            point.pointid = unique_id[cid]
            touched_interior[point.pointid] = True
            cid += 1
        breakpts = _tools.np.logical_and(touched_interior, breakpts)

        # break curves on repeated points
        def split_curves(curves, breakpts, loops=[]):
            newcurves = []
            newmap = []
            for curve in curves:
                breaks = _tools.np.where(breakpts[curve.pointsid[1:-1]])[0]
                if len(breaks) == 0:
                    newcurves.append(curve)
                    newmap.append([len(newcurves)-1])
                else:
                    breaks = [0] + list(breaks+1) + [len(curve.points)-1]
                    nc = []
                    for i, j in zip(breaks[:-1], breaks[1:]):
                        ncurve = _Curve(
                            curve.points[i:j+1], curve.tag, curve.curve_type,
                            curve.projection)
                        ncurve.pointsid = curve.pointsid[i:j+1]
                        newcurves.append(ncurve)
                        nc.append(len(newcurves)-1)
                    newmap.append(nc)
            for iloop, loop in enumerate(loops):
                newloop = []
                for oi, ori in loop:
                    if ori:
                        for j in newmap[oi]:
                            newloop.append((j, ori))
                    else:
                        for j in reversed(newmap[oi]):
                            newloop.append((j, ori))
                loops[iloop] = newloop
            return newcurves

        p2curve = _tools.np.full([self._points.shape[0], 2, 2], -1, int)
        curve2p = _tools.np.ndarray([len(self._curves), 2], int)
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
        touched = _tools.np.full((len(self._curves)), False, _tools.np.bool)
        self._curveloops = []
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
                    assert(not(i != self._curveloops[-1][0][0] and touched[i]))
                    touched[i] = True
                    assert(i != -1)
        except AssertionError:
            raise ValueError("Invalid topology")
        self._curves = split_curves(self._curves, breakpts, self._curveloops)
        self._interior_curves = split_curves(self._interior_curves, breakpts)
        loopbboxarea = _tools.np.zeros([len(self._curveloops)])
        for i, l in enumerate(self._curveloops):
            lpts = _tools.np.row_stack([self._curves[j].points for j, o in l])
            bbox = _tools.np.max(lpts, axis=0)-_tools.np.min(lpts, axis=0)
            loopbboxarea[i] = bbox[0]*bbox[1]
        self._curveloops.insert(
            0, self._curveloops.pop(_tools.np.argmax(loopbboxarea)))

    def _add_geometry(self, geometry, tag, projection, curve_type, interior,
                      onlypoints=False):
        if geometry.GetGeometryCount() != 0:
            for subg in geometry:
                self._add_geometry(subg, tag, projection,
                                   curve_type, interior, onlypoints)
            return
        if onlypoints:
            self.add_interior_points(geometry.GetPoints(), tag, projection)
        else:
            if interior:
                self.add_interior_curve(geometry.GetPoints(), tag,
                                        projection, curve_type)
            else:
                self.add_boundary_curve(geometry.GetPoints(), tag,
                                        projection, curve_type)

    def _add_shapefile(self, filename, physical_name_field,
                       interior, points, curve_type):
        progress = _tools.ProgressLog(
                    "Import features from \"{}\"".format(filename), True)
        if filename[-5:] == ".gpkg":
            driver = _tools.ogr.GetDriverByName('GPKG')
        else:
            driver = _tools.ogr.GetDriverByName('ESRI Shapefile')
        data = driver.Open(filename, 0)
        layer = data.GetLayer()
        layerdef = layer.GetLayerDefn()
        physfield = None
        count = 0
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
            if i.geometry() is None :
                continue
            phys = (i.GetField(physfield)
                    if not (physfield is None) else "boundary")
            self._add_geometry(i.geometry(), phys, layerproj, curve_type,
                               interior, points)
            progress.log("{} features imported".format(count))
            count += 1
        progress.end()

    def add_interior_points(self, points: _tools.np.ndarray, physical_tag: str,
                            projection: _tools.osr.SpatialReference) -> None:
        """ Add forced interior mesh points

        Args:
            x: the points [n,2]
            physical_tag: the points physical tag
            projection: the points coordinate system
        """
        x = _tools.project_points(points, projection, self._projection)
        points = list(_Point(p, physical_tag) for p in x)
        self._interior_points += points

    def add_boundary_curve(self, points: _tools.np.ndarray, physical_tag: str,
                           projection: _tools.osr.SpatialReference,
                           curve_type: CurveType = CurveType.POLYLINE) -> None:
        """ Add a tagged curve to the domain boundary.

        Args:
            points: the curve control points [n,2]
            physical_tag: the curve physical tag
            projection: the points coordinate system
            curve_type: curve interpolation
        """
        curve = _Curve(points, physical_tag, curve_type, projection)
        self._curves.append(curve)

    def add_interior_curve(self, points: _tools.np.ndarray, physical_tag: str,
                           projection: _tools.osr.SpatialReference,
                           curve_type: CurveType = CurveType.POLYLINE) -> None:
        """ Adds a tagged curve inside the domain. The curve
        is not part of the domain boundary and will be meshed on both sides.

        Args:
            points: the curve control points [n,2]
            physical_tag: the curve physical tag
            projection: the points coordinate system
            curve_type: curve interpolation
        """
        curve = _Curve(points, physical_tag, curve_type, projection)
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

def coarsen_boundaries(domain: Domain, x0: _tools.Tuple[float, float],
                       x0_projection: _tools.osr.SpatialReference,
                       mesh_size: MeshSizeCallback) -> Domain:
    """ Creates a new Domain with the same projection and coarsened
    boundaries.

    Args:
        domain: the domain to coarsen
        x0: the coordinates of one point inside the domain.
        x0_projection: the coordinates system of x0.
        mesh_size: a function returning the desired mesh element size for given
            coordinates
    """
    _tools.log("Coarsen boundaries", True)
    x0 = _tools.project_points(_tools.np.array([x0]), x0_projection,
                               domain._projection)[0]
    sampled = []
    tags = []
    maxtag = 1
    str2tag = {}

    def mesh_size_half(x, p):
        return mesh_size(x, p)*0.5

    progress = _tools.ProgressLog("Sampling curves for coarsening")
    for icurve, curve in enumerate(domain._curves):
        cs = _curve_sample(curve, mesh_size_half, domain._projection)
        sampled.append(cs)
        if curve.tag not in str2tag:
            str2tag[curve.tag] = maxtag
            maxtag += 1
        tags.append(_tools.np.full(cs.shape[0], str2tag[curve.tag],
                                   dtype=_tools.np.int32))
        progress.log("{} curves sampled".format(icurve+1))
    progress.end()
    x = _tools.np.vstack(sampled)
    tags = _tools.np.concatenate(tags)
    x, unique_id, _ = _generate_unique_points(x)
    utags = list([-1, -1] for i in x)
    for i, tag in zip(unique_id, tags):
        if utags[i][0] == -1:
            utags[i][0] = tag
        found = False
        j = i
        while j != -1:
            found |= utags[j][0] == tag
            j = utags[j][1]
        if not found:
            utags[j][1] = len(utags)
        utags.append([tag, -1])
    tags = _tools.np.array(utags).reshape([-1])
    x = _tools.np.copy(x[:, :2])
    # avoid cocircular points
    eps = (_tools.np.max(x, axis=0, keepdims=True) -
           _tools.np.min(x, axis=0, keepdims=True))*1e-8
    x = x + _tools.np.random.random(x.shape)*eps
    _tools.log("Delaunay mesh of sampled points")
    tri = _tools.Delaunay(x)
    _tools.log("Extract boundaries")
    first = tri.find_simplex(x0)
    if (first == -1):
        raise(ValueError("First point outside domain"))
    n_l = _tools.c.c_int()
    n_xo = _tools.c.c_int()
    p_xo = _tools.c.POINTER(_tools.c.c_double)()
    p_l = _tools.c.POINTER(_tools.c.c_int)()
    ms = mesh_size(x, domain._projection)
    _tools.lib.gen_boundaries_from_points(
        _tools.c.c_int(x.shape[0]), _tools.np2c(x),
        _tools.np2c(tags, _tools.np.int32),
        _tools.c.c_int(tri.simplices.shape[0]),
        _tools.np2c(tri.simplices, _tools.np.int32),
        _tools.c.c_int(first), _tools.np2c(ms),
        _tools.c.byref(p_xo), _tools.c.byref(n_xo),
        _tools.c.byref(p_l), _tools.c.byref(n_l))
    # xptr = _tools.c.POINTER(n_xo.value*2*_tools.c.c_double)
    # xbuf = _tools.c.cast(p_xo, xptr).contents
    xtype = _tools.c.POINTER(n_xo.value*2*_tools.c.c_double)
    xbuf = _tools.c.cast(p_xo, xtype).contents
    xo = _tools.np.ctypeslib.frombuffer(xbuf, dtype=_tools.np.float64)
    xo = xo.reshape([-1, 2]).copy()
    _tools.lib.libcfree(p_xo)
    linesbuf = _tools.c.cast(p_l, _tools.c.POINTER(n_l.value*_tools.c.c_int))
    lines = _tools.np.ctypeslib.frombuffer(linesbuf.contents,
                                           dtype=_tools.np.int32)
    odomain = Domain(domain._projection)
    breaks = _tools.np.where(lines == -1)[0]
    tagi2str = dict((i, s) for (s, i) in str2tag.items())
    for i, b in enumerate(breaks):
        ifrom = 0 if i == 0 else breaks[i-1]+1
        ito = breaks[i]
        tag = tagi2str[lines[ifrom]]
        pts = lines[ifrom+1:ito]
        odomain.add_boundary_curve(xo[pts], tag,
                                   domain._projection,
                                   curve_type=CurveType.POLYLINE)
    _tools.lib.libcfree(p_l)
    return odomain
