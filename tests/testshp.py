# seamsh - Copyright (C) <2010-2020>
# <Universite catholique de Louvain (UCL), Belgium
# 	
# List of the contributors to the development of seamsh: see AUTHORS file.
# Description and complete License: see LICENSE file.
# 	
# This program (seamsh) is free software: 
# you can redistribute it and/or modify it under the terms of the GNU Lesser General 
# Public License as published by the Free Software Foundation, either version
# 3 of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License
# along with this program (see COPYING file).  If not, 
# see <http://www.gnu.org/licenses/>.

import os
import urllib.request
import tarfile
if not os.path.isdir("data") :
    urllib.request.urlretrieve("ftp://braque.mema.ucl.ac.be/seamsh/data-test-1.tar.gz", "data-test-1.tar.gz")
    f = tarfile.open("data-test-1.tar.gz","r:*")
    f.extractall()

import seamsh
from seamsh.geometry import CurveType
import seamsh.geometry
import numpy as np
from osgeo import osr 

def mesh_size(x,projection) :
    s_coast = np.clip((dist_coast_field(x,projection)-400)*0.5,200,5000)
    s_porquerolles = np.clip((dist_porquerolles_field(x,projection)-200)*0.5,50,5000)
    return np.minimum(s_coast,s_porquerolles)

domain_srs = osr.SpatialReference()

domain_srs.ImportFromEPSG(32631)
#domain_srs.ImportFromProj4("+proj=utm +ellps=WGS84 +zone=31")

domain = seamsh.geometry.Domain(domain_srs)
domain.add_boundary_curves_shp("data/data_no_duplicate.shp","physical",CurveType.POLYLINE)
#domain.add_interior_curves_shp("data/interior.shp",None,CurveType.STRICTPOLYLINE)
#domain.add_interior_points_shp("data/interior.shp")
bath_field = seamsh.field.Raster("data/medit.tiff")
dist_coast_field = seamsh.field.Distance(domain,100,["coast","island"])
dist_porquerolles_field = seamsh.field.Distance(domain,20,["porquerolles"])
#coarse = seamsh.geometry.coarsen_boundaries(domain,(8e5,4.68e6),domain_srs,mesh_size,20)
seamsh.gmsh.mesh(domain,"test.msh",mesh_size,intermediate_file_name="log")
seamsh.gmsh.convert_to_gis("test.msh",domain_srs,"test.gpkg")

seamsh.gmsh.gmsh.model.add("test")
seamsh.gmsh.gmsh.open("test.msh")
tag,nodes = seamsh.gmsh.gmsh.model.mesh.getElementsByType(2)
ntri = len(tag)
print("ntri",ntri)
assert(ntri>40000 and ntri<41000)
