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
# But for now, the coarsening algorithm presents several limitations:
#   - The physical tags are lost (this will be fixed soon).
#   - Interior curves hare not handled, they can be present as long as they
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

domain_srs = osr.SpatialReference()
domain_srs.ImportFromProj4("+proj=utm +ellps=WGS84 +zone=31")
domain = seamsh.geometry.Domain(domain_srs)
domain.add_boundary_curves_shp("data/data_no_duplicate.shp","physical",CurveType.POLYLINE)

dist_coast_field = seamsh.field.Distance(domain,100,["coast","island"])
dist_porquerolles_field = seamsh.field.Distance(domain,20,["porquerolles"])

def mesh_size(x,projection) :
    s_coast = np.clip((dist_coast_field(x,projection)-400)*0.5,200,5000)
    s_porq = np.clip((dist_porquerolles_field(x,projection)-200)*0.5,50,5000)
    s_dist = np.minimum(s_coast,s_porq)
    return s_dist

# %%
# The geometry.coarsen_boundaries function requires in input the coordinates of one (any)
# point which inside the domain. The last parameters should be at least two times smaller
# than the smallest value returned by the mesh_size callback, otherwise the resulting
# geometry will be invalid.

coarse = seamsh.geometry.coarsen_boundaries(domain,(8e5,4.68e6),domain_srs,mesh_size,20)

# %%
# The resulting coarsend domain can be meshed normally.

seamsh.gmsh.mesh(coarse,"coarse_boundary_mesh.msh",mesh_size)

