from osgeo import ogr,osr
import numpy as np
from scipy.spatial import cKDTree
import logging
import gmsh

POLYLINE = 0
BSPLINE = 1
SPLINE = 2
STRICTPOLYLINE = 3

class _Curve :

    def __init__(self,points,tag,curve_type):
        self.points = points
        self.mesh_size = None
        self.tag = tag
        self.curve_type = curve_type


class Domain :

    def __init__(self, projection) :
        self._projection = projection
        self._interior_curves = []
        self._curves = []

    def add_curve(self, points, tag, projection, interior=False, curve_type=BSPLINE) :
        points = np.asarray(points)
        assert(len(points.shape) == 2)
        if not projection.IsSame(self._projection) :
            logging.info("reprojecting coordinates")
            points = osr.CoordinateTransformation(projection,self._projection).TransformPoints(points)
            points = np.asarray(points)
        if points.shape[1] == 2 :
            points = np.column_stack([points,np.zeros((points.shape[0],1))])
        curve = _Curve(points,tag,curve_type)
        if interior :
            self._interior_curves.append(curve)
        else :
            self._curves.append(curve)

    def build_topology(self) :
        allpoints = np.row_stack([c.points for c in self._curves])
        bbmin = np.min(allpoints,axis=0)
        bbmax = np.max(allpoints,axis=0)
        eps = np.linalg.norm(bbmax-bbmin)*1e-8
        epoints = np.row_stack([np.array([c.points[0],c.points[-1]]) for c in self._curves])
        nc = len(self._curves)
        curve_to_p = np.ndarray((nc,2),np.int32)
        touched = np.full((nc),False,np.bool)
        p_to_curve = []
        self._curveloops = []
        for ip, (p0,p1) in enumerate(cKDTree(epoints).query_pairs(eps)) :
            p_to_curve.append(((p0//2,p0%2),(p1//2,p1%2)))
            curve_to_p[p0//2,p0%2] = ip
            curve_to_p[p1//2,p1%2] = ip
        try :
            for i,curve in enumerate(self._curves) :
                print(i,len(self._curves))
                if touched[i] : continue
                self._curveloops.append([])
                j = 0
                while (not self._curveloops[-1]) or self._curveloops[-1][0][0] != i:
                    self._curveloops[-1].append((i,j==0))
                    p = curve_to_p[i,(j+1)%2]
                    i,j = p_to_curve[p][0] if p_to_curve[p][0][0] != i else p_to_curve[p][1]
                    touched[i] = True
        except :
            raise ValueError("Invalid topology")
        loopbboxarea = np.zeros([len(self._curveloops)])
        for i,loop in enumerate(self._curveloops) :
            lpts = np.row_stack([self._curves[j].points for j,o in loop])
            bbox = np.max(lpts,axis=0)-np.min(lpts,axis=0)
            loopbboxarea[i] = bbox[0]*bbox[1]
        self._curveloops.insert(0,self._curveloops.pop(np.argmax(loopbboxarea)))

    def _add_geometry(self,geometry,tag,projection,curve_type) :
        if geometry.GetGeometryCount() != 0 :
            for subg in geometry :
                self._add_geometry(subg,tag,projection,curve_type)
            return
        assert(geometry.GetGeometryType() == 2)
        self.add_curve(geometry.GetPoints(),tag,projection,curve_type)

    def add_shapefile(self, filename, interior=False, curve_type=POLYLINE) :
        proj = osr.SpatialReference()
        driver = ogr.GetDriverByName('ESRI Shapefile')
        data = driver.Open(filename,0)
        layer = data.GetLayer()
        layerproj = layer.GetSpatialRef()
        for i in layer :
            g = i.geometry()
            self._add_geometry(g,1,layerproj,curve_type)

    def clip_circle(self, center, radius, resolution) :
        pass

    def unrefine_geometry(self) :
        pass


def _gen_mesh_gmsh(sf_callback,adapt_step=2) :
    if xo.shape[1] == 2 :
        xo = np.column_stack([xo,np.zeros([xo.shape[0],1])])
    lcpoint = sf_callback(xo)
    gmsh.option.setNumber("Mesh.CharacteristicLengthFromPoints",0)
    gmsh.option.setNumber("Mesh.CharacteristicLengthFactor",1)
    gmsh.option.setNumber("Mesh.CharacteristicLengthFromCurvature",0)
    gmsh.option.setNumber("Mesh.CharacteristicLengthExtendFromBoundary",0)
    gmsh.option.setNumber("Mesh.CharacteristicLengthFromParametricPoints",1)
    gmsh.model.mesh.generate(1)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin",initial_sampling)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax",initial_sampling)
    gmsh.model.mesh.generate(2)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin",0)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax",1e22)
    bg_field = gmsh.model.mesh.field.add("PostView");
    gmsh.model.mesh.field.setAsBackgroundMesh(bg_field);
    for i in range(adapt_step) :
        node_tags, node_x, _  = gmsh.model.mesh.getNodes()
        node_x = node_x.reshape([-1,3])
        nodes_map = dict({tag:i for i,tag in enumerate(node_tags)})
        sf_view = gmsh.view.add("mesh size field")
        node_lc = sf_callback(node_x)
        tri = gmsh.model.mesh.getElementsByType(2)[1].reshape([-1,3])
        tri = np.vectorize(nodes_map.get)(tri)
        data = np.column_stack([node_x[tri,:].swapaxes(1,2).reshape([-1,9]),node_lc[tri]])
        gmsh.view.addListData(sf_view,"ST",tri.shape[0],data.reshape([-1]))
        gmsh.model.mesh.field.setNumber(bg_field, "ViewTag", sf_view);
        gmsh.model.mesh.generate(2)
        gmsh.view.remove(sf_view)
    gmsh.model.mesh.field.remove(bg_field)
    if not (outfilename is None):
        gmsh.write(outfilename)

def mesh_gmsh(domain,filename,mesh_size,version=4.0,gmsh_options={}) :
    nadapt = 3
    for curve in domain._curves :
        curve.mesh_size = mesh_size(curve.points)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("Mesh.CharacteristicLengthFromPoints",0)
    gmsh.option.setNumber("Mesh.CharacteristicLengthFactor",1)
    gmsh.option.setNumber("Mesh.CharacteristicLengthFromCurvature",0)
    gmsh.option.setNumber("Mesh.CharacteristicLengthExtendFromBoundary",0)
    gmsh.option.setNumber("Mesh.CharacteristicLengthFromParametricPoints",1)
    ## occ geometry mesh ##
    allctag = []
    cltag = []
    physicals = {}
    for cl in domain._curveloops :
        clpts = list([gmsh.model.occ.addPoint(
            *domain._curves[l].points[0 if o else -1,:]) for l,o in cl])
        ctag = []
        for i,(l,o) in enumerate(cl) :
            curve = domain._curves[l]
            ipts = list([gmsh.model.occ.addPoint(*p) for p in curve.points[1:-1]])
            if not o :
                ipts.reverse()
            pts = [clpts[i]]+ipts+[clpts[(i+1)%len(cl)]]
            if curve.curve_type == SPLINE :
                tag = gmsh.model.occ.addSpline(pts)
                physicals.setdefault(curve.tag,[]).append(tag)
                ctag.append(tag)
            elif curve.curve_type == BSPLINE :
                tag = gmsh.model.occ.addBSpline(pts)
                physicals.setdefault(curve.tag,[]).append(tag)
                ctag.append(tag)
                u = np.linspace(0,1,curve.points.shape[0])
                gmsh.model.mesh.setSizeAtParametricPoints(1,ctag[-1],u,curve.mesh_size)
            elif curve.curve_type == STRICTPOLYLINE :
                for p0,p1 in zip(pts[:-1],pts[1:]) :
                    tag = gmsh.model.occ.addLine(p0,p1)
                    physicals.setdefault(curve.tag,[]).append(tag)
                    ctag.append(tag)
            elif curve.curve_type == POLYLINE :
                tag = gmsh.model.occ.addBSpline(pts,degree=1)
                physicals.setdefault(curve.tag,[]).append(tag)
                ctag.append(tag)
        cltag.append(gmsh.model.occ.addWire(ctag))
        allctag += ctag
    gmsh.model.occ.addPlaneSurface(cltag)
    gmsh.model.occ.synchronize()
    it = 0
    for cl in domain._curveloops :
        for i,(l,o) in enumerate(cl) :
            curve = domain._curves[l]
            sizes = curve.mesh_size if o else np.flip(curve.mesh_size)
            if curve.curve_type in [BSPLINE,SPLINE,POLYLINE] :
                u0,u1 = gmsh.model.getParametrizationBounds(1,allctag[it])
                u = np.linspace(u0,u1,curve.points.shape[0])
                gmsh.model.mesh.setSizeAtParametricPoints(1,allctag[it],u,sizes)
                it += 1
            elif curve.curve_type == STRICTPOLYLINE :
                for j in range(curve.points.shape[0]-1) :
                    gmsh.model.mesh.setSizeAtParametricPoints(1,allctag[it],[0,1],[sizes[j],sizes[j+1]])
                    it += 1
    for name,tags in physicals.items() :
        tag = gmsh.model.addPhysicalGroup(1,tags)
        gmsh.model.setPhysicalName(1,tag,name)
    ## 1D mesh ##
    gmsh.model.mesh.generate(1)
    for i in range(nadapt) :
        for (dim,tag) in gmsh.model.getEntities(1) :
            _,x,u = gmsh.model.mesh.getNodes(dim,tag,includeBoundary=True)
            size = mesh_size(x.reshape([-1,3]))
            gmsh.model.mesh.setSizeAtParametricPoints(dim,tag,u,size)
        gmsh.model.mesh.generate(1)
    ## 2D mesh ##
    gmsh.model.mesh.generate(2)
    bg_field = gmsh.model.mesh.field.add("PostView");
    gmsh.model.mesh.field.setAsBackgroundMesh(bg_field);
    for i in range(nadapt) :
        node_tags, node_x, _  = gmsh.model.mesh.getNodes()
        node_x = node_x.reshape([-1,3])
        nodes_map = dict({tag:i for i,tag in enumerate(node_tags)})
        sf_view = gmsh.view.add("mesh size field")
        node_lc = mesh_size(node_x)
        tri = gmsh.model.mesh.getElementsByType(2)[1].reshape([-1,3])
        tri = np.vectorize(nodes_map.get)(tri)
        data = np.column_stack([node_x[tri,:].swapaxes(1,2).reshape([-1,9]),node_lc[tri]])
        gmsh.view.addListData(sf_view,"ST",tri.shape[0],data.reshape([-1]))
        gmsh.model.mesh.field.setNumber(bg_field, "ViewTag", sf_view);
        gmsh.model.mesh.generate(2)
        gmsh.view.remove(sf_view)
    gmsh.model.mesh.field.remove(bg_field)
    gmsh.fltk.run()
    gmsh.option.setNumber("Mesh.MshFileVersion",version)
    gmsh.write(filename)
    gmsh.finalize()

def convert_msh_shp() :
    pass

