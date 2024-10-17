# %%
# GIS Data
# ========
# This example illustrates how the lirbary interacts with various GIS file
# formats directly from numpy arrays and osr spatial references.
#
# Let's download some data.

import os
import urllib.request
import tarfile
if not os.path.isdir("data2"):
    url = "https://nextcloud.cism.ucl.ac.be/s/TfY6tCzrqiKaTTX/download/data-test-2.tar.gz"
    urllib.request.urlretrieve(url, "data-test-2.tar.gz")
    f = tarfile.open("data-test-2.tar.gz", "r:*")
    f.extractall()

# %%
# The actual example starts here.

import seamsh
from seamsh.geometry import CurveType
import seamsh.geometry
import numpy as np
from osgeo import osr

# %%
# First, let's create a domain object and its associated projection.
# Any projection supported by osgeo.osr can be used.

domain_srs = osr.SpatialReference()
domain_srs.ImportFromProj4("+proj=utm +ellps=WGS84 +zone=31")
domain = seamsh.geometry.Domain(domain_srs)

# %%
# We load all curves from a given ESTRI shapefile as polylines.
# In the shapefile, a field named "physical" defines the physical tag of each
# curve. If a re-projection is required, it will be done automatically.

domain.add_boundary_curves_shp("data2/data_no_duplicate.shp",
                               "physical", CurveType.POLYLINE)

# %%
# Interior curves and interior points can be loaded in a similar way.
# Physical tags are optional

domain.add_interior_curves_shp("data2/interior.shp",
                               None, CurveType.STRICTPOLYLINE)
domain.add_interior_points_shp("data2/interior_points.shp", "physical")

# %%
# Seamsh provides helper classes to compute the element size field.
#
# field.Raster allows to load geotiff files (or any gdal raster file).

bath_field = seamsh.field.Raster("data2/medit.tiff")

# %%
# field.Distance can be used to compute the distance from given (tagged)
# boundaries. A first field returns the distance from boundaries with physical
# tag "coast", or "island". The curves are sampled at regular intervals and
# the computed distance is actually the distance from the closest sampling
# point. The second argument sets the length of the interval between the
# sampling points.

dist_coast = seamsh.field.Distance(domain, 100, ["coast", "island"])

# %%
# A second distance field returns the distance from the island tagged
# "porquerolles", in the shp file.

dist_porquerolles = seamsh.field.Distance(domain, 20, ["porquerolles"])


# %%
# The mesh size is defined based on those fields, smaller elements are
# prescribed around the Porquerolles island.

def mesh_size(x, projection):
    s_coast = np.clip((dist_coast(x, projection)-400)*0.5, 200, 5000)
    s_porq = np.clip((dist_porquerolles(x, projection)-200)*0.5, 50, 5000)
    s_dist = np.minimum(s_coast, s_porq)
    return s_dist


# %%
# Another option would be to take a mesh size proportional to the square root
# of the (clipped) bathymetry :

# def mesh_size(x,projection) :
#     bath = -bath_field(x,projection)
#     s_bath = np.sqrt(np.clip(bath,100,4000))*10
#     return  s_bath*10

# %%
# The "intermediate_file_name" option is used to save files containing
# intermediate meshes and mesh size fields. If this parameter takes the
# special value "-" an interactive gmsh graphical window will pop up
# after each meshing step.

output_srs = osr.SpatialReference()
output_srs.ImportFromProj4("+proj=latlon +ellps=WGS84")

seamsh.gmsh.mesh(domain, "gis_mesh.msh", mesh_size,
                 intermediate_file_name="debug", output_srs=output_srs)

seamsh.gmsh.reproject("gis_mesh.msh", None, "gis_mesh_lonlat.msh", output_srs)

# %%
# The gmsh.convert_to_gis function can be used to convert a gmsh .msh file
# into a shape file or into a geo package file.

seamsh.gmsh.convert_to_gis("gis_mesh.msh", None, "gis_mesh.gpkg")
