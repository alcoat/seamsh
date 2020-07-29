# %%
# Advanced Gmsh / Quad Mesh
# =========================
# Seamsh uses gmsh internally. The whole gmsh python api is exposed through the
# seamsh.gmsh.gmsh module (the repetition is not a typo). This example is a
# variation on the previous one using the gmsh api to generate a quadrangular
# instead of triangular mesh.
#
# We start as usual.

import os
import urllib.request
import tarfile
if not os.path.isdir("data"):
    url = "ftp://braque.mema.ucl.ac.be/seamsh/data-test-1.tar.gz"
    urllib.request.urlretrieve(url, "data-test-1.tar.gz")
    f = tarfile.open("data-test-1.tar.gz", "r:*")
    f.extractall()


import seamsh
from seamsh.geometry import CurveType
import seamsh.geometry
import numpy as np
from osgeo import osr

domain_srs = osr.SpatialReference()
domain_srs.ImportFromProj4("+proj=utm +ellps=WGS84 +zone=31")
domain = seamsh.geometry.Domain(domain_srs)
domain.add_boundary_curves_shp("data/data_no_duplicate.shp",
                               "physical", CurveType.POLYLINE)

dist_coast = seamsh.field.Distance(domain, 100, ["coast", "island"])
dist_porquerolles = seamsh.field.Distance(domain, 20, ["porquerolles"])


def mesh_size(x, projection):
    s_coast = np.clip((dist_coast(x, projection)-400)*0.5, 200, 5000)
    s_porq = np.clip((dist_porquerolles(x, projection)-200)*0.5, 50, 5000)
    s_dist = np.minimum(s_coast, s_porq)
    return s_dist


coarse = seamsh.geometry.coarsen_boundaries(domain, (8e5, 4.68e6),
                                            domain_srs, mesh_size)

# %%
# The gmsh options are set. The full list of gmsh mesh options can be found in
# `gmsh's documentation <https://gmsh.info/doc/texinfo/gmsh.html#Mesh-options-list>`_.

seamsh.gmsh.gmsh.option.setNumber("Mesh.RecombineAll", 1)
seamsh.gmsh.gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 1)
seamsh.gmsh.gmsh.option.setNumber("Mesh.Algorithm", 8)

# %%
# Eventually the mesh is generated.

seamsh.gmsh.mesh(coarse, "quad_mesh.msh", mesh_size)
