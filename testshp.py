import msea
import numpy as np
from osgeo import osr 

def mesh_size(x) :
    s_coast = np.clip((dist_coast_field(x)-400)*0.5,200,5000)
    s_porquerolles = np.clip((dist_porquerolles_field(x)-200)*0.5,50,5000)
    return np.minimum(s_coast,s_porquerolles)

domain_srs = osr.SpatialReference()

domain_srs.ImportFromEPSG(32631)
#domain_srs.ImportFromProj4("+proj=utm +ellps=WGS84 +zone=31")

domain = msea.Domain(domain_srs)
domain.add_shapefile("test/data_no_duplicate.shp","physical",curve_type=msea.POLYLINE)
#domain.add_shapefile("test/interior.shp",None,True,curve_type=msea.POLYLINE)
#domain.add_interior_points("test/interior.shp")
dist_coast_field = msea.DistanceField(domain,100,["coast","island"])
dist_porquerolles_field = msea.DistanceField(domain,20,["porquerolles"])
domain.unrefine_boundaries(8e5,4.68e6,20,mesh_size)
domain.build_topology()
msea.mesh_gmsh(domain,"test.msh",mesh_size)
