import gmsh
import numpy as np
import uuid
import typing
from .geometry import Domain, CurveType, MeshSizeCallback
from osgeo import ogr, osr
import os
import struct
import atexit

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)
GmshOptionValue = typing.Union[float, str, typing.Tuple[float, float, float]]


def _gmsh_curve_geo(curve_type: CurveType, pointsid):
    pts = list(i+1 for i in pointsid)
    if curve_type == CurveType.SPLINE:
        return [gmsh.model.geo.addSpline(pts)]
    elif curve_type == CurveType.BSPLINE:
        return [gmsh.model.geo.addBSpline(pts)]
    elif curve_type == CurveType.STRICTPOLYLINE:
        print("STRICT POLYLINE !!!")
        pairs = zip(pts[:-1], pts[1:])
        return list(gmsh.model.geo.addLine(p0, p1) for p0, p1 in pairs)
    elif curve_type == CurveType.POLYLINE:
        return [gmsh.model.geo.addPolyline(pts)]
    else:
        raise ValueError("unknown curve_type '%s'" % curve_type)


def _curve_sample(curve, lc):
    gmsh.model.add(str(uuid.uuid4()))
    tags = list([gmsh.model.geo.addPoint(*x, 0) -
                 1 for x in curve.points[:-1, :]])
    if np.linalg.norm(curve.points[0, :]-curve.points[-1, :]) < 1e-8:
        tags.append(tags[0])
    else:
        tags.append(gmsh.model.geo.addPoint(*curve.points[-1, :], 0)-1)
    ltag = _gmsh_curve_geo(curve.curve_type, tags)
    gmsh.model.geo.synchronize()
    gmsh.option.setNumber("Mesh.LcIntegrationPrecision", 1e-3)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", lc)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", lc)
    gmsh.model.mesh.generate(1)
    nodes = gmsh.model.mesh.getNodes(1, ltag[0], includeBoundary=True)[
        1].reshape([-1, 3])
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 1e22)
    gmsh.model.remove()
    return nodes[:, :2]


def mesh(domain: Domain, filename: str, mesh_size: MeshSizeCallback,
         version: float = 4.0,
         gmsh_options: typing.Dict[str, GmshOptionValue] = {}) -> None:
    """ Use gmsh to generate a mesh from a geometry and a mesh size callback

    Args:
        domain: the input geometry
        filename: output mesh file (.msh)
        mesh_size: callbable prescribing the mesh element size
        version: msh file version (typically 2.0 or 4.0)
        gmsh_options: additional options passed to gmsh
    """

    domain._build_topology()
    nadapt = 3
    nadapt1d = 4
    for curve in domain._curves:
        curve.mesh_size = mesh_size(curve.points, domain._projection)
    gmsh.option.setNumber("Mesh.CharacteristicLengthFromPoints", 0)
    gmsh.option.setNumber("Mesh.LcIntegrationPrecision", 1e-3)
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
            ctag.extend(tags)
        cltag.append(gmsh.model.geo.addCurveLoop(ctag))
        ctags += ctag
    stag = gmsh.model.geo.addPlaneSurface(cltag)
    embeded = []
    for curve in domain._interior_curves:
        tags = _gmsh_curve_geo(curve.curve_type, curve.pointsid)
        embeded.extend(tags)
        physicals.setdefault(curve.tag, []).extend(tags)
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.embed(1, embeded, 2, stag)
    gmsh.model.mesh.embed(
        0, [i+1 for i in domain._interior_points_id], 2, stag)
    it = 0
    for cl in domain._curveloops:
        for i, (l, o) in enumerate(cl):
            curve = domain._curves[l]
            sizes = curve.mesh_size if o else np.flip(curve.mesh_size)
            u0, u1 = gmsh.model.getParametrizationBounds(1, ctags[it])
            u = np.linspace(u0, u1, curve.points.shape[0])
            gmsh.model.mesh.setSizeAtParametricPoints(
                1, ctags[it], u, sizes)
            it += 1
            """
            elif curve.curve_type == STRICTPOLYLINE:
                for j in range(curve.points.shape[0]-1):
                    gmsh.model.mesh.setSizeAtParametricPoints(
                        1, ctags[it], [0, 1], [sizes[j], sizes[j+1]])
                    it += 1
            """
    for name, tags in physicals.items():
        tag = gmsh.model.addPhysicalGroup(1, tags)
        gmsh.model.setPhysicalName(1, tag, name)
    tag = gmsh.model.addPhysicalGroup(2, [stag])
    gmsh.model.setPhysicalName(2, stag, "domain")
    # 1D mesh
    gmsh.model.mesh.generate(1)
    for i in range(nadapt1d):
        for (dim, tag) in gmsh.model.getEntities(1):
            _, x, u = gmsh.model.mesh.getNodes(dim, tag, True)
            size = mesh_size(x.reshape([-1, 3]), domain._projection)
            gmsh.model.mesh.setSizeAtParametricPoints(dim, tag, u, size)
        print("pass ", i, nadapt)
        gmsh.model.mesh.clear()
        gmsh.model.mesh.generate(1)
    # gmsh.fltk.run()
    # 2D mesh ##
    initial_sampling = 0.1
    gmsh.model.mesh.generate(1)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", initial_sampling)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", initial_sampling)
    gmsh.option.setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.CharacteristicLengthFromParametricPoints",
                          0)
    gmsh.model.mesh.generate(2)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 1e22)
    bg_field = gmsh.model.mesh.field.add("PostView")
    gmsh.model.mesh.field.setAsBackgroundMesh(bg_field)
    for i in range(nadapt):
        node_tags, node_x, _ = gmsh.model.mesh.getNodes(-1, -1)
        node_x = node_x.reshape([-1, 3])
        nodes_map = dict({tag: i for i, tag in enumerate(node_tags)})
        sf_view = gmsh.view.add("mesh size field")
        node_lc = mesh_size(node_x, domain._projection)
        tri = gmsh.model.mesh.getElements(2, stag)[2][0]
        tri = np.array(list(nodes_map[t] for t in tri)).reshape([-1, 3])
        xnode = node_x[tri, :].swapaxes(1, 2).reshape([-1, 9])
        data = np.column_stack([xnode, node_lc[tri]]).reshape([-1])
        gmsh.view.addListData(sf_view, "ST", tri.shape[0], data)
        gmsh.model.mesh.field.setNumber(bg_field, "ViewTag", sf_view)
        gmsh.model.mesh.generate(2)
        gmsh.view.remove(sf_view)
    gmsh.model.mesh.field.remove(bg_field)
    gmsh.fltk.run()
    gmsh.option.setNumber("Mesh.MshFileVersion", version)
    gmsh.write(filename)


def convert_to_gis(input_filename: str, projection: osr.SpatialReference,
        output_filename: str) ->None:
    """ Convert a triangular gmsh mesh into shapefile or geopackage.

    Args:
        input_filename : any mesh file readable by gmsh (typically
            a .msh file)
        projection: the projection assigned to the output layer, mesh
            files do not store any projection so neither checks nor
            re-projection are performed.
        output_filename : shape file (.shp) or geopackage (.gpkg) file
    """
    gmsh.model.add(str(uuid.uuid4()))
    gmsh.open(input_filename)
    if output_filename.endswith(".shp"):
        shpdriver = ogr.GetDriverByName('ESRI Shapefile')
    elif output_filename.endswith(".gpkg"):
        shpdriver = ogr.GetDriverByName('GPKG')
    else:
        raise ValueError("Unknown file extension '" + output_filename+"'")
    if os.path.exists(output_filename):
        shpdriver.DeleteDataSource(output_filename)
    out_data_source = shpdriver.CreateDataSource(output_filename)
    out_layer = out_data_source.CreateLayer(output_filename,
                                            geom_type=ogr.wkbPolygon,
                                            srs=projection)
    out_layer.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))
    id_field = out_layer.FindFieldIndex("id", 1)

    tri_i, tri_n = gmsh.model.mesh.getElementsByType(2)
    node_tags, x, _ = gmsh.model.mesh.getNodesByElementType(2)
    x = x.reshape([-1, 3])
    nodes_map = dict({tag: i for i, tag in enumerate(node_tags)})
    tri_n = np.array(list(nodes_map[t] for t in tri_n)).reshape([-1, 3])

    tri_n = np.column_stack([tri_n, tri_n[:, [0]]])
    header = struct.pack("<biii", 1, 3, 1, 4)
    out_data_source.StartTransaction()
    for i, tn in zip(tri_i, tri_n):
        wkb = header+x[tn, :2].astype('<f8').tobytes()
        feat = ogr.Feature(out_layer.GetLayerDefn())
        feat.SetGeometryDirectly(ogr.CreateGeometryFromWkb(wkb))
        feat.SetField(id_field, int(i))
        out_layer.CreateFeature(feat)
        feat = None
    out_data_source.CommitTransaction()
    gmsh.model.remove()


@atexit.register
def finalize():
    print("finalizing")
    gmsh.finalize()
