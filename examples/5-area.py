# %%
# Area
# ========
# This example revisits the GIS example (GIS Data) with use of a specific area
# to adapt locally the meshsize.
# In fact we made a mix between
# - the distance from the coast and porquerolle
# - the distance computed from the bathymetry
#
# This example starts as the previous one.

import os
import urllib.request
import tarfile
if not os.path.isdir("data2"):
    url = "https://nextcloud.cism.ucl.ac.be/s/L62bAwF9kma3Ta3/download/data-test-2.tar.gz"
    urllib.request.urlretrieve(url, "data-test-2.tar.gz")
    f = tarfile.open("data-test-2.tar.gz", "r:*")
    f.extractall()

import seamsh
from seamsh.geometry import CurveType
import seamsh.geometry
import numpy as np
from osgeo import osr

domain_srs = osr.SpatialReference()
domain_srs.ImportFromProj4("+proj=utm +ellps=WGS84 +zone=31")
domain = seamsh.geometry.Domain(domain_srs)
domain_poly = seamsh.geometry.Domain(domain_srs)

domain.add_boundary_curves_shp("data2/data_no_duplicate.shp",
                               "physical", CurveType.POLYLINE)

domain.add_interior_curves_shp("data2/interior.shp",
                               None, CurveType.STRICTPOLYLINE)
domain.add_interior_points_shp("data2/interior_points.shp", "physical")

bath_field = seamsh.field.Raster("data2/medit.tiff")

domain_poly.add_boundary_curves_shp("data2/area_utm31.gpkg", "name", CurveType.POLYLINE)

dist_coast = seamsh.field.Distance(domain, 100, ["coast", "island"])
dist_porquerolles = seamsh.field.Distance(domain, 20, ["porquerolles"])
in_area = seamsh.field.Inpoly(domain_poly, ["noraster"])


def mesh_size(x, projection):
    s_coast = np.clip((dist_coast(x, projection)-400)*0.5, 200, 5000)
    s_porq = np.clip((dist_porquerolles(x, projection)-200)*0.5, 50, 5000)
    s_dist = np.minimum(s_coast, s_porq)
    bath = -bath_field(x, projection)
    s_bath = np.sqrt(np.clip(bath,100,4000))*100
    s_mix = np.where(in_area(x, projection), s_dist, s_bath)
    return s_mix


output_srs = osr.SpatialReference()
output_srs.ImportFromProj4("+proj=latlon +ellps=WGS84")

seamsh.gmsh.mesh(domain, "gis_mesh.msh", mesh_size,
                 intermediate_file_name="debug", output_srs=output_srs)


seamsh.gmsh.convert_to_gis("gis_mesh.msh", output_srs, "gis_mesh.gpkg")
