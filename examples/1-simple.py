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

# %%
# my super example simple
# =======================

import seamsh
import numpy as np
from osgeo import osr 

def mesh_size(x) :
    return (0.0025+(1+np.sin(4*np.pi*x[:,1]*np.pi*x[:,0]))*0.01)*2

domain_srs = osr.SpatialReference()
domain_srs.ImportFromEPSG(4326)

domain = seamsh.Domain(domain_srs)
domain.add_curve([[0,0],[-0.2,0.25],[-0.2,0.75],[0,1]],
        "wall",domain_srs,curve_type=seamsh.SPLINE)
domain.add_curve([[0,1],[0.25,1.2],[0.75,1.2],[1,1]],
        "wall1",domain_srs,curve_type=seamsh.BSPLINE)
domain.add_curve([[0,0],[0.25,-0.2],[0.75,-0.2],[1,0]],
        "wall2",domain_srs,curve_type=seamsh.STRICTPOLYLINE)
domain.add_curve([[1,1],[1.2,0.75],[1.2,0.25],[1,0]],
        "wall3",domain_srs,curve_type=seamsh.STRICTPOLYLINE)
domain.add_curve([[0.4,0.6],[0.6,0.6],[0.6,0.4],[0.4,0.4],[0.4,0.6]],
        "wallin",domain_srs,curve_type=seamsh.BSPLINE)
domain.build_topology()

seamsh.mesh_gmsh(domain,"test.msh",mesh_size)
