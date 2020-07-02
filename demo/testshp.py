import msea
from msea.geometry import CurveType
import numpy as np
from osgeo import osr 

def mesh_size(x,projection) :
    s_coast = np.clip((dist_coast_field(x,projection)-400)*0.5,200,5000)
    s_porquerolles = np.clip((dist_porquerolles_field(x,projection)-200)*0.5,50,5000)
    return np.minimum(s_coast,s_porquerolles)

domain_srs = osr.SpatialReference()

domain_srs.ImportFromEPSG(32631)
#domain_srs.ImportFromProj4("+proj=utm +ellps=WGS84 +zone=31")

domain = msea.geometry.Domain(domain_srs)
domain.add_boundary_curves_shp("test/data_no_duplicate.shp","physical",CurveType.POLYLINE)
#domain.add_interior_curves_shp("test/interior.shp",None,CurveType.STRICTPOLYLINE)
#domain.add_interior_points_shp("test/interior.shp")
bath_field = msea.field.Raster("medit.tiff")
dist_coast_field = msea.field.Distance(domain,100,["coast","island"])
dist_porquerolles_field = msea.field.Distance(domain,20,["porquerolles"])
#coarse = msea.geometry.coarsen_boundaries(domain,(8e5,4.68e6),domain_srs,mesh_size,20)
msea.gmsh.mesh(domain,"test.msh",mesh_size,intermediate_file_name="log")
msea.gmsh.convert_to_gis("test.msh",domain_srs,"test.gpkg")
