# %%
# Simple Geometry
# ===============
# This example illustrates how to insert different types of curves
# directly from numpy arrays and osr spatial references.

import seamsh
from seamsh.geometry import CurveType
import numpy as np
from osgeo import osr

# %%
# First,  let's create a domain object and its associated projection.
# Any projection supported by osgeo.osr can be used.

domain_srs = osr.SpatialReference()
domain_srs.ImportFromEPSG(4326)
domain = seamsh.geometry.Domain(domain_srs)

# %%
# Boundary curves are created from a list of control points, an osr projection,
# a physical tag and a curve type. If the projection does not match the
# domain's projection,  the points are re-projected.
#
# A STRICTPOLYLINE curve linearly interpolates control points and forces a mesh
# point on each control point.

domain.add_boundary_curve([[0, 0], [-0.2, 0.25], [-0.2, 0.75], [0, 1]],
                          "wall0", domain_srs,
                          curve_type=CurveType.STRICTPOLYLINE)

# %%
# A POLYLINE is similar to a STRICTPOLYLINE but no mesh point is forced on the
# control points,  the mesher can cut the corners.

x = [[0, 1], [0.25, 1.2], [0.75, 1.2], [1, 1]]
domain.add_boundary_curve(x, "wall1", domain_srs, CurveType.POLYLINE)

# %%
# A SPLINE uses a cubic spline interpolation through the control points

x = [[0, 0], [0.25, -0.2], [0.75, -0.2], [1, 0]]
domain.add_boundary_curve(x, "wall2", domain_srs, CurveType.SPLINE)

# %%
# A BSPLINE is smoother than a SPLINE but the control points are not on the
# curve.

x = [[1, 1], [1.2, 0.75], [1.2, 0.25], [1, 0]]
domain.add_boundary_curve(x, "wall3", domain_srs, CurveType.BSPLINE)

# %%
# If the last point is the same as the first one,  a periodic curve is
# created.

x = [[0.4, 0.6], [0.6, 0.6], [0.6, 0.4], [0.4, 0.4], [0.4, 0.6]]
domain.add_boundary_curve(x, "island", domain_srs, CurveType.BSPLINE)

# %%
# Interior curves are meshed on both sides,  they do not define the domain
# boundaries but forces edge alignment and mesh points. When an interior curve
# touch a boundary or another interior curve,  a control point should be
# present in both curves at the intersection.

x = [[-0.2, 0.75], [0.1, 0.6], [0.2, 0.4], [0.1, 0.4]]
domain.add_interior_curve(x, "interior1", domain_srs, CurveType.BSPLINE)
x = [[0.1, 0.4], [0.1, 0.6], [0.1, 0.8]]
domain.add_interior_curve(x, "interior2", domain_srs, CurveType.POLYLINE)

# %%
# Force several tagged mesh points

x = [[0.8, 0.13]]
domain.add_interior_points(x, "forcedpoint1", domain_srs)
x = [[0.8, 0.8], [0.8, 0.7]]
domain.add_interior_points(x, "forcedpoint2", domain_srs)


# %%
# A callback which takes an numpy array of points and a projection as argument
# defines the mesh element size. In this case,  a simple analytical function is
# used.

def mesh_size(x, projection):
    return (0.005+(1+np.sin(4*np.pi*x[:, 1]*np.pi*x[:, 0]))*0.01)


# %%
# Finally,  the seamesh.gmsh.mesh function generates the mesh.

seamsh.gmsh.mesh(domain, "basic_geometry.msh", mesh_size, version=2.0)
