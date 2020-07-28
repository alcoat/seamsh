# %%
# Boundary Coarsening
# ===================
# This example revisits the previous example (GIS Data) with an extra step to
# reduce the resolution of the input coastlines.
#
# This approach has two main advantages:
#
#   - Thanks to the fitting of the input data resolution to the final mesh
#     element size, the quality of the elements near the boundary is improved.
#   - The input curves do not have to form topologically valid closed loops.
#     The only requirement is that the external boundary should be closed 
#     but the lines defining the boundaries can intersect each others. 
# 
# But :
#
#   - Interior curves are not handled, they can be present as long as they
#     do not intersect the domain boundaries (external and islands) but they
#     won't be coarsened.
#
# This example starts as the previous one.

import os
import urllib.request
import tarfile
if not os.path.isdir("data") :
    urllib.request.urlretrieve("ftp://braque.mema.ucl.ac.be/seamsh/data-test-1.tar.gz",
                               "data-test-1.tar.gz")
    f = tarfile.open("data-test-1.tar.gz","r:*")
    f.extractall()

import seamsh
from seamsh.geometry import CurveType
import seamsh.geometry
import numpy as np
from osgeo import osr 

lonlat = osr.SpatialReference()
lonlat.ImportFromProj4("+proj=latlong +ellps=WGS84 +unit=degrees")
utm32 = osr.SpatialReference()
utm32.ImportFromProj4("+proj=utm +ellps=WGS84 +zone=32")
utm31 = osr.SpatialReference()
utm31.ImportFromProj4("+proj=utm +ellps=WGS84 +zone=31")

domain = seamsh.geometry.Domain(utm31)
domain.add_boundary_curves_shp("data/data_no_duplicate.shp","physical",CurveType.POLYLINE)

# %%
# Clip the domain by a circle in UTM (zone32) coordinates

alphas = np.linspace(0,2*np.pi,200)
circle_points = np.array([[272e3,4763e3]])+50e3*np.column_stack([np.sin(alphas),np.cos(alphas)])
domain.add_boundary_curve(circle_points,"open",utm32)

dist_coast_field = seamsh.field.Distance(domain,100,["coast","island"])
dist_porquerolles_field = seamsh.field.Distance(domain,20,["porquerolles"])

def mesh_size(x,projection) :
    s_coast = np.clip((dist_coast_field(x,projection)-400)*0.5,200,5000)
    s_porq = np.clip((dist_porquerolles_field(x,projection)-200)*0.5,50,5000)
    s_dist = np.minimum(s_coast,s_porq)
    return s_dist

# %%
# The geometry.coarsen_boundaries function requires in input the coordinates of one (any)
# point which is inside the domain. The last parameters should be at least two times smaller
# than the smallest value returned by the mesh_size callback, otherwise the resulting
# geometry will be invalid.

coarse = seamsh.geometry.coarsen_boundaries(domain,(6.21,42.95),lonlat,mesh_size,20)

# %%
# The resulting coarsened domain can be meshed normally.

seamsh.gmsh.mesh(coarse,"coarse_boundary_mesh.msh",mesh_size)

