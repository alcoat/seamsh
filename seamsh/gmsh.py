# seamsh - Copyright (C) <2010-2020>
# <Universite catholique de Louvain (UCL), Belgium
#
# List of the contributors to the development of seamsh: see AUTHORS file.
# Description and complete License: see LICENSE file.
#
# This program (seamsh) is free software:
# you can redistribute it and/or modify it under the terms of the GNU
# Lesser General Public License as published by the Free Software Foundation,
# either version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program (see COPYING file).  If not,
# see <http://www.gnu.org/licenses/>.

import gmsh
from . import geometry as _geometry
from . import _tools
__all__ = ["mesh", "convert_to_gis"]

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)
gmsh.option.setNumber("General.Verbosity", 2)


def _gmsh_curve_geo(curve_type: _geometry.CurveType, pointsid):
    pts = list(i+1 for i in pointsid)
    if curve_type == _geometry.CurveType.SPLINE:
        return [gmsh.model.geo.addSpline(pts)]
    elif curve_type == _geometry.CurveType.BSPLINE:
        return [gmsh.model.geo.addBSpline(pts)]
    elif curve_type == _geometry.CurveType.STRICTPOLYLINE:
        pairs = zip(pts[:-1], pts[1:])
        return list(gmsh.model.geo.addLine(p0, p1) for p0, p1 in pairs)
    elif curve_type == _geometry.CurveType.POLYLINE:
        return [gmsh.model.geo.addPolyline(pts)]
    else:
        raise ValueError("unknown curve_type '%s'" % curve_type)


def _curve_sample_gmsh_tag(tag, lc, projection):
    bounds = gmsh.model.getParametrizationBounds(1, tag)
    xi = _tools.np.linspace(bounds[0][0], bounds[1][0], 5)
    count = 0
    x = gmsh.model.getValue(1, tag, xi).reshape([-1, 3])
    size = lc(x, projection)
    while True:
        count += 1
        dist = _tools.np.linalg.norm(x[1:, :]-x[:-1, :], axis=1)
        target = _tools.np.minimum(size[1:], size[:-1])
        if not _tools.np.any(target < dist):
            break
        exi = (xi[1:][target < dist]+xi[:-1][target < dist])/2
        ex = gmsh.model.getValue(1, tag, exi).reshape([-1, 3])
        esize = lc(ex, projection)
        xi = _tools.np.hstack([xi, exi])
        s = _tools.np.argsort(xi)
        xi = xi[s]
        x = _tools.np.row_stack([x, ex])[s, :]
        size = _tools.np.hstack([size, esize])[s]
    return x[:, :2], xi, size


def _curve_sample(curve, lc, projection):
    gmsh.model.add(str(_tools.uuid.uuid4()))
    tags = list([gmsh.model.geo.addPoint(*x, 0) -
                 1 for x in curve.points[:-1, :]])
    if _tools.np.linalg.norm(curve.points[0, :]-curve.points[-1, :]) < 1e-8:
        tags.append(tags[0])
    else:
        tags.append(gmsh.model.geo.addPoint(*curve.points[-1, :], 0)-1)
    ltag = _gmsh_curve_geo(curve.curve_type, tags)
    gmsh.model.geo.synchronize()
    r = _curve_sample_gmsh_tag(ltag[0], lc, projection)[0]
    gmsh.model.remove()
    return r

def _create_gmsh_geometry(domain: _geometry.Domain) :
    _tools.log("Build gmsh model")
    gmsh.option.setNumber("Mesh.CharacteristicLengthFromPoints", 0)
    gmsh.option.setNumber("Mesh.LcIntegrationPrecision", 1e-5)
    gmsh.option.setNumber("Mesh.CharacteristicLengthFactor", 1)
    gmsh.option.setNumber("Mesh.CharacteristicLengthFromCurvature", 0)
    gmsh.option.setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.CharacteristicLengthFromParametricPoints",
                          1)
    ctags = []
    cltag = []
    physicals = {}
    for i, x in enumerate(domain._points):
        gmsh.model.geo.addPoint(*x, 0, tag=i+1)
    for cl in domain._curveloops:
        ctag = []
        for i, (l, o) in enumerate(cl):
            curve = domain._curves[l]
            tags = _gmsh_curve_geo(curve.curve_type, curve.pointsid)
            physicals.setdefault(curve.tag, []).extend(tags)
            ctag.extend(tags if o else list(-t for t in reversed(tags)))
        cltag.append(gmsh.model.geo.addCurveLoop(ctag))
        ctags += ctag
    stag = gmsh.model.geo.addPlaneSurface(cltag)
    embeded = []
    for curve in domain._interior_curves:
        tags = _gmsh_curve_geo(curve.curve_type, curve.pointsid)
        embeded.extend(tags)
        physicals.setdefault(curve.tag, []).extend(tags)
    physicalspt = {}
    for point in domain._interior_points:
        physicalspt.setdefault(point.tag, []).append(point.pointid+1)
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.embed(1, embeded, 2, stag)
    gmsh.model.mesh.embed(
        0, [p.pointid+1 for p in domain._interior_points], 2, stag)
    for name, tags in physicals.items():
        tag = gmsh.model.addPhysicalGroup(1, tags)
        gmsh.model.setPhysicalName(1, tag, name)
    for name, tags in physicalspt.items():
        tag = gmsh.model.addPhysicalGroup(0, tags)
        gmsh.model.setPhysicalName(0, tag, name)
    tag = gmsh.model.addPhysicalGroup(2, [stag])
    gmsh.model.setPhysicalName(2, stag, "domain")


def _mesh_bgrid(domain: _geometry.Domain, mesh_size: _geometry.MeshSizeCallback, smoothness:float):
    _tools.log("Build smooth mesh size field")
    np = _tools.np
    x0 = np.min(domain._points,axis=0)
    x1 = np.max(domain._points,axis=0)
    d = np.max(x1-x0)/100
    x = np.arange(x0[0],x1[0]+d,d)
    y = np.arange(x0[1],x1[1]+d,d)
    xx,yy = np.meshgrid(x,y)
    xy = np.stack((xx.T,yy.T),axis=2).reshape(-1,2)
    v = mesh_size(xy,domain._projection)
    poly = _tools.PolyMesh(x0, x1)
    poly.add_points(xy)
    while True:
        xyz = np.c_[xy.reshape(-1,2),np.zeros((xy.reshape(-1,2).shape[0],1))]
        tri = poly.get_triangles()
        edges = np.r_[tri[:,[0,1]],tri[:,[1,2]],tri[:,[2,0]]]
        edges.sort(axis=1)
        shift = 2**32
        ekey = np.unique(edges[:,0]*shift+edges[:,1])
        edges = np.c_[ekey//shift,ekey%shift]
        dedges = xy[edges[:,0]]-xy[edges[:,1]]
        elength = np.hypot(dedges[:,0], dedges[:,1])
        while True:
            vo = v.copy()
            cond = v[edges[:,0]]>v[edges[:,1]]+elength*smoothness
            v[edges[cond,0]] = v[edges[cond,1]]+elength[cond]*smoothness
            cond = v[edges[:,1]]>v[edges[:,0]]+elength*smoothness
            v[edges[cond,1]] = v[edges[cond,0]]+elength[cond]*smoothness
            if np.max(np.abs(v-vo)) < 1:
                break
        eminsize = np.minimum(v[edges[:,0]],v[edges[:,1]])
        tocut = np.where(eminsize*2<elength)[0]
        if tocut.size ==0:
            break
        newp = (xy[edges[tocut,0]]+xy[edges[tocut,1]])/2
        xy = np.r_[xy, newp]
        v = np.r_[v, mesh_size(newp, domain._projection)]
        poly.add_points(newp)
    xtri = xyz[tri]
    data = np.c_[xtri.swapaxes(1,2).reshape(-1,9),v.flatten()[tri]]
    view = gmsh.view.add("tri")
    gmsh.view.add_list_data(view, "ST", tri.shape[0], data.flatten())
    gmsh.fltk.run()
    field = gmsh.model.mesh.field.add("PostView")
    gmsh.model.mesh.field.setNumber(field, "ViewTag",view)
    gmsh.model.mesh.field.setAsBackgroundMesh(field)
    gmsh.option.setNumber("Mesh.CharacteristicLengthFromPoints", 0)
    gmsh.option.setNumber("Mesh.LcIntegrationPrecision", 1e-5)
    gmsh.option.setNumber("Mesh.CharacteristicLengthFactor", 1)
    gmsh.option.setNumber("Mesh.CharacteristicLengthFromCurvature", 0)
    gmsh.option.setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.CharacteristicLengthFromParametricPoints", 0)
    _tools.log("Mesh with gmsh")
    gmsh.model.mesh.generate(2)


def _mesh_successive(domain: _geometry.Domain, mesh_size: _geometry.MeshSizeCallback, intermediate_file_name: str=None):
    for curve in domain._curves:
        curve.mesh_size = mesh_size(curve.points, domain._projection)

    # 1D mesh
    progress = _tools.ProgressLog("Sample curves for mesh size")
    for icurve, (dim, tag) in enumerate(gmsh.model.getEntities(1)):
        _, xi, size = _curve_sample_gmsh_tag(tag,
                                             lambda x, p: mesh_size(x, p)/2,
                                             domain._projection)
        gmsh.model.mesh.setSizeAtParametricPoints(dim, tag, xi, size*2)
        progress.log("{} curve sampled".format(icurve+1))
    progress.end()
    _tools.log("Generate 1D mesh")
    gmsh.model.mesh.generate(1)
    if intermediate_file_name is not None:
        if intermediate_file_name == "-":
            gmsh.fltk.run()
        else:
            gmsh.write(intermediate_file_name+"_1d.msh")
    # 2D mesh ##
    _tools.log("Generate initial 2D mesh")
    _, x, u = gmsh.model.mesh.getNodes()
    x = x.reshape([-1, 3])
    domain_size = _tools.np.max(x, axis=0)-_tools.np.min(x, axis=0)
    initlc = _tools.np.linalg.norm(domain_size)/100
    gmsh.option.setNumber("Mesh.CharacteristicLengthFromParametricPoints", 0)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", initlc)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", initlc)
    gmsh.model.mesh.generate(2)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 1e22)
    bg_field = gmsh.model.mesh.field.add("PostView")
    gmsh.model.mesh.field.setAsBackgroundMesh(bg_field)
    nadapt = 3
    for i in range(nadapt):
        _tools.log("Generate refined 2D mesh (pass {}/{})".format(i+1, nadapt))
        node_tags, node_x, _ = gmsh.model.mesh.getNodes(-1, -1)
        node_x = node_x.reshape([-1, 3])
        nodes_map = dict({tag: i for i, tag in enumerate(node_tags)})
        sf_view = gmsh.view.add("mesh size field")
        node_lc = mesh_size(node_x, domain._projection)
        etypes, etags, enodes = gmsh.model.mesh.getElements(2, -1)
        for etype, enode in zip(etypes, enodes):
            nn = 3 if etype == 2 else 4
            enode = list(nodes_map[t] for t in enode)
            enode = _tools.np.array(enode).reshape([-1, nn])
            xnode = node_x[enode, :].swapaxes(1, 2).reshape([-1, nn*3])
            data = _tools.np.column_stack([xnode, node_lc[enode]])
            gmsh.view.addListData(sf_view, "ST"if etype == 2 else "SQ",
                                  enode.shape[0], data.reshape([-1]))
        gmsh.model.mesh.field.setNumber(bg_field, "ViewTag", sf_view)
        if intermediate_file_name is not None:
            if intermediate_file_name == "-":
                gmsh.fltk.run()
            else:
                gmsh.view.write(sf_view,
                                intermediate_file_name+"_2d_"+str(i)+".pos")
        gmsh.model.mesh.generate(2)
        gmsh.view.remove(sf_view)
    gmsh.model.mesh.field.remove(bg_field)

def mesh(domain: _geometry.Domain, filename: str,
         mesh_size: _geometry.MeshSizeCallback,
         version: float = 4.0, intermediate_file_name: str=None, smoothness=-1) -> None:
    """ Calls gmsh to generate a mesh from a geometry and a mesh size callback

    Args:
        domain: the input geometry
        filename: output mesh file (.msh)
        mesh_size: callbable prescribing the mesh element size
        version: msh file version (typically 2.0 or 4.0)
        smoothness: if positive, controls the maximum gradation of the mesh size, 0.25 is
            a good choice. When used, the mesh size is recursively evaluated on a regular grid
            which may slow down the meshing process.
        intermediate_file_name: if not None, save intermediate meshes to those
            files for debugging purpose (suffixes and extensions will be
            appended), if == "-", an interactive gmsh graphical window will pop
            up after each meshing step.
    """
    _tools.log("Generate mesh", True)
    domain._build_topology()
    _create_gmsh_geometry(domain)
    if smoothness>0:
        _mesh_bgrid(domain, mesh_size, smoothness)
    else:
        _mesh_successive(domain, mesh_size, intermediate_file_name)

    # remove nodes that do not touch any triangle
    _,keepnodes = gmsh.model.mesh.getElementsByType(2)
    keepnodes = set(keepnodes)
    rm = []
    for e in gmsh.model.getEntities(1) :
        enodes,_,_ = gmsh.model.mesh.getNodes(*e)
        if (len(enodes) == 0) or (not (enodes[0] in keepnodes)) :
            rm.append(e)
    if len(rm) != 0 :
        gmsh.model.removeEntities(rm,True)
    _tools.log("Write \"{}\" (msh version {})".format(filename, version))
    gmsh.option.setNumber("Mesh.MshFileVersion", version)
    gmsh.write(filename)
    gmsh.model.remove()


def convert_to_gis(input_filename: str,
                   projection: _tools.osr.SpatialReference,
                   output_filename: str) -> None:
    """ Converts a triangular gmsh mesh into shapefile or geopackage.

    Args:
        input_filename : any mesh file readable by gmsh (typically
            a .msh file)
        projection: the projection assigned to the output layer, mesh
            files do not store any projection so neither checks nor
            re-projections are performed.
        output_filename : shape file (.shp) or geopackage (.gpkg) file
    """
    _tools.log("Convert \"{}\" into \"{}\"".format(input_filename,
               output_filename), True)
    gmsh.model.add(str(_tools.uuid.uuid4()))
    gmsh.open(input_filename)
    if output_filename.endswith(".shp"):
        shpdriver = _tools.ogr.GetDriverByName('ESRI Shapefile')
    elif output_filename.endswith(".gpkg"):
        shpdriver = _tools.ogr.GetDriverByName('GPKG')
    else:
        raise ValueError("Unknown file extension '" + output_filename+"'")
    if _tools.os.path.exists(output_filename):
        shpdriver.DeleteDataSource(output_filename)
    out_data_source = shpdriver.CreateDataSource(output_filename)
    out_layer = out_data_source.CreateLayer(output_filename,
                                            geom_type=_tools.ogr.wkbPolygon,
                                            srs=projection)
    out_layer.CreateField(_tools.ogr.FieldDefn("id", _tools.ogr.OFTInteger))
    id_field = out_layer.FindFieldIndex("id", 1)

    tri_i, tri_n = gmsh.model.mesh.getElementsByType(2)
    node_tags, x, _ = gmsh.model.mesh.getNodesByElementType(2)
    x = x.reshape([-1, 3])
    nodes_map = dict({tag: i for i, tag in enumerate(node_tags)})
    tri_n = _tools.np.array(list(nodes_map[t] for t in tri_n)).reshape([-1, 3])

    tri_n = _tools.np.column_stack([tri_n, tri_n[:, [0]]])
    header = _tools.struct.pack("<biii", 1, 3, 1, 4)
    out_data_source.StartTransaction()
    for i, tn in zip(tri_i, tri_n):
        wkb = header+x[tn, :2].astype('<f8').tobytes()
        feat = _tools.ogr.Feature(out_layer.GetLayerDefn())
        feat.SetGeometryDirectly(_tools.ogr.CreateGeometryFromWkb(wkb))
        feat.SetField(id_field, int(i))
        out_layer.CreateFeature(feat)
        feat = None
    out_data_source.CommitTransaction()
    gmsh.model.remove()


@_tools.atexit.register
def _finalize():
    gmsh.finalize()
