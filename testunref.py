import msea
import numpy as np
from osgeo import osr 

def mesh_size(x) :
    d = dist(x)
    return np.clip((d-0.01)*0.5,0.01,0.1)/2

domain_srs = osr.SpatialReference()
domain_srs.ImportFromEPSG(4326)

domain = msea.Domain(domain_srs)
#domain.add_shapefile("test/test_polygons.shp")
domain.add_shapefile("test/data_no_duplicate.shp","physical",curve_type=msea.POLYLINE)
dist = msea.DistanceField(domain,["coast","island"])
domain.unrefine_boundaries(42,6,mesh_size)
#
#domain.add_shapefile("test/interior.shp",None,True,curve_type=msea.POLYLINE)
#domain.add_interior_points("test/interior.shp")
domain.build_topology()
msea.mesh_gmsh(domain,"test.msh",mesh_size)
