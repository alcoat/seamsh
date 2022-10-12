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
from .geometry import Domain as _Domain
from .gmsh import _curve_sample


class Distance:
    """Callable evaluating the distance to a set of discretized curves.

    The curves are discretized as sets of points, then the distance to the
    closest point is computed.
    """

    def __init__(self, domain: _Domain, sampling: float,
                 tags: _tools.List[str] = None):
        """
        Args:
            domain: a Domain object containing the set of curves
            sampling: the interval between two consecutive sampling points.
            tags: List of physical tags specifying the curve from the domain.
                if None, all curves are taken into account.
        """
        _tools.log("Create distance field", True)
        points = []
        msg = "Sampling features for distance computation"
        progress = _tools.ProgressLog(msg)
        all_curves_iter = _tools.chain(domain._curves, domain._interior_curves)

        def size(x, proj):
            return _tools.np.full([x.shape[0]], sampling)

        icurve = 0
        for curve in all_curves_iter:
            if (tags is None) or (curve.tag in tags):
                points.append(_curve_sample(curve, size, domain._projection))
                icurve += 1
                progress.log("{} features sampled".format(icurve))
        for point in _tools.chain(domain._interior_points):
            if (tags is None) or (point.tag in tags):
                points.append(point.x)
        progress.end()
        points = _tools.np.vstack(points)
        _tools.log("Build KDTree with {} points".format(points.shape[0]))
        self._tree = _tools.cKDTree(points)
        self._projection = domain._projection

    def __call__(self, x: _tools.np.ndarray,
                 projection: _tools.osr.SpatialReference
                 ) -> _tools.np.ndarray:
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
        msg = "Create field from raster file \"{}\"".format(filename)
        _tools.log(msg, True)
        src_ds = _tools.gdal.Open(filename)
        self._geo_matrix = src_ds.GetGeoTransform()
        self._data = src_ds.GetRasterBand(1).ReadAsArray()
        self._projection = _tools.osr.SpatialReference()
        self._projection.ImportFromWkt(src_ds.GetProjection())
        if int(_tools.gdal.__version__.split(".")[0]) >= 3:
            order = _tools.osr.OAMS_TRADITIONAL_GIS_ORDER
            self._projection.SetAxisMappingStrategy(order)

    def __call__(self, x: _tools.np.ndarray,
                 projection: _tools.osr.SpatialReference
                 ) -> _tools.np.ndarray:
        """Evaluates the field value on each point of x.

        Keyword arguments:
            x: the points [n,2]
            projection: the coordinate system of the points
        Returns:
            The field value on points x. [n]
        """
        x = _tools.project_points(x, projection, self._projection)
        gm = self._geo_matrix
        lon, lat = x[:, 0], x[:, 1]
        det = gm[5]*gm[1]-gm[2]*gm[4]
        pixx = ((lon-gm[0])*gm[5]-(lat-gm[3])*gm[2])/det
        pixy = ((lat-gm[3])*gm[1]-(lon-gm[0])*gm[4])/det
        pixy = _tools.np.clip(pixy.astype(int),0,self._data.shape[0]-1)
        pixx = _tools.np.clip(pixx.astype(int),0,self._data.shape[1]-1)
        return self._data[pixy, pixx]

class Inpoly:
    """Callable evaluating the membership to a set of polygons

    """

    def __init__(self, domain: _Domain,
                 tags: _tools.List[str] = None):
        """
        Args:
            domain: a Domain object containing the set of curves
            tags: List of physical tags specifying the curve from the domain.
                if None, all curves are taken into account.
        """
        if not _tools.shapely_available:
            raise ValueError("The shapely python module is required to use Inpoly fields.")
        _tools.log("Initialisation Inpoly", True)
        self._projection = domain._projection

        domain._build_topology()
        self._area = []
        for ll in domain._curveloops:
            if (tags is None) or (ll.tag in tags):
                ptsid = []
                for lid, o in ll :
                    l = domain._curves[lid]
                    ptsid.extend(l.pointsid if o else np.flip(l.pointsid))
                assert(ptsid[-1] == ptsid[0])
                ptsid = ptsid[:-1]
                pts = domain._points[ptsid]
                self._area.append(_tools.shapely.geometry.Polygon(pts))

    def __call__(self, x: _tools.np.ndarray,
                 projection: _tools.osr.SpatialReference
                 ) -> _tools.np.ndarray:
        """Indicate for each point of x if it belongs to a set of polygons.

        Args:
            x: the points [n,2]
            projection: the coordinate system of the points, should be
                the same coordinate system as the domain, otherwise no
                conversion is done and an exception is raised.
        Returns:
            True if the point belongs to a polygon / False otherwise. [n]
        """
        x = _tools.project_points(x, projection, self._projection)

        xs = []
        inpoly = _tools.np.full(x.shape[0], False)
        xs = list(_tools.shapely.geometry.Point(xi[0], xi[1]) for xi in x)
        for polygon in self._area:
            for i, xsi in enumerate(xs):
                if not inpoly[i]:
                    inpoly[i] = xsi.within(polygon)

        return inpoly

