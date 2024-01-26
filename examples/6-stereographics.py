# %%
# Stereographic projection
# ========================
# This example illustrates the use of the stereographic projection to generate 
# a mesh covering most of the Earth's surface.

# %%
# Download the Natural Earth coastline data set and extract the shapefile.

import requests
import zipfile
import os

coastline_url = "https://naciscdn.org/naturalearth/10m/physical/ne_10m_coastline.zip"
shp_file = "ne_10m_coastline/ne_10m_coastline.shp"

if not os.path.isfile(shp_file):
    req = requests.get(coastline_url, allow_redirects=True)
    open("ne_10m_coastline.zip", "wb").write(req.content)
    zipfile.ZipFile("ne_10m_coastline.zip").extractall("ne_10m_coastline")

import seamsh
from seamsh.geometry import CurveType
import seamsh.geometry
import numpy as np
from osgeo import osr

# %%
# The mesh is generated in the polar stereographic projection.
# the singular point, at the opposite of the projection origin has to be outside the domain. Here
# the the North Pole is chosen as the origin of the stereographic projection so that the singular
# point (i.e. the South Pole) is outside the meshed domain.

domain_srs = osr.SpatialReference()
domain_srs.ImportFromProj4("+ellps=WGS84 +proj=stere +lat_0=90")
domain = seamsh.geometry.Domain(domain_srs)
domain.add_boundary_curves_shp(shp_file, "featurecla", CurveType.POLYLINE)

# %%
# "Cartesian" projection (i.e. no projection) is used to compute the distance to the coast.

cart_srs = osr.SpatialReference()
cart_srs.ImportFromProj4("+ellps=WGS84 +proj=cart +units=m +x_0=0 +y_0=0")


dist_coast = seamsh.field.Distance(domain, 10000, projection=cart_srs)


# %%
# it is necessary to convert the mesh size to stereographic coordinates

def mesh_size(x, projection):
    s_coast = np.clip((dist_coast(x, projection)-20000)*0.5, 20000, 100000)
    R = 6371000
    stereo_factor = 2/(1+x[:,0]**2/R**2+x[:,1]**2/R**2)
    return s_coast/stereo_factor


coarse = seamsh.geometry.coarsen_boundaries(domain, (0, 0), domain_srs, mesh_size)

# %%
# Eventually the mesh is generated.
# The mesh is saved both in stereographic and cartesian coordinates

seamsh.gmsh.mesh(coarse, "natural_earth.msh", mesh_size, output_srs=domain_srs)
seamsh.gmsh.reproject("natural_earth.msh", domain_srs, "natural_earth_cart.msh", cart_srs)
