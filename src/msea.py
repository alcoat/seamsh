from osgeo import ogr,osr,gdal
import platform 
import numpy as np
from scipy.spatial import cKDTree,Delaunay
from scipy import interpolate
import logging
import gmsh
import os
import struct
import uuid
import itertools
import ctypes as c

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)


libdir = os.path.dirname(os.path.realpath(__file__))
if platform.system() == "Windows":
    libpath = os.path.join(libdir, "unref.dll")
elif platform.system() == "Darwin":
    libpath = os.path.join(libdir, "libunref.dylib")
else:
    libpath = os.path.join(libdir, "libmsea.so")

lib = c.CDLL(libpath)


POLYLINE = 0
BSPLINE = 1
SPLINE = 2
STRICTPOLYLINE = 3

def _gmsh_geo(curve_type,pointsid) :
    pts = list(i+1 for i in pointsid)
    if curve_type == SPLINE :
        return [gmsh.model.geo.addSpline(pts)]
    elif curve_type == BSPLINE :
        return [gmsh.model.geo.addBSpline(pts)]
    elif curve_type == STRICTPOLYLINE :
        return list(gmsh.model.geo.addLine(p0,p1)for p0,p1 in zip(pts[:-1],pts[1:]))
    elif curve_type == POLYLINE :
        return [gmsh.model.geo.addPolyline(pts)]


def _generate_unique_points(x) :
    tree = cKDTree(x)
    unique_id = np.full([x.shape[0]],-1,np.int32)
    bbmin = np.min(x,axis=0)
    bbmax = np.max(x,axis=0)
    eps = np.linalg.norm(bbmax-bbmin)*1e-12
    eps = 1e-12
    cid = 0
    for p0,p1 in tree.query_pairs(eps) :
        uid = max(unique_id[p0],unique_id[p1])
        if uid == -1 :
            uid = cid
            cid += 1
        unique_id[p0] = uid
        unique_id[p1] = uid
    nbreakpoints = cid
    for i,u in enumerate(unique_id) :
        if unique_id[i] == -1 :
            unique_id[i] = cid
            cid += 1
    xo  = np.ndarray([cid,2])
    xo[unique_id,:] = x
    return xo,unique_id,nbreakpoints

class _Curve :

    def __init__(self,points,tag,curve_type):
        self.points = points[:,:2]
        self.mesh_size = None
        self.tag = tag
        self.curve_type = curve_type
    
    def sample(self,lc) :
        gmsh.model.add(str(uuid.uuid4()))
        tags = list([gmsh.model.geo.addPoint(*x,0)-1 for x in self.points[:-1,:]])
        if np.linalg.norm(self.points[0,:]-self.points[-1,:]) < 1e-8 :
            tags.append(tags[0])
        else :
            tags.append(gmsh.model.geo.addPoint(*self.points[-1,:],0)-1)
        ltag = _gmsh_geo(self.curve_type,tags)
        gmsh.model.geo.synchronize()
        gmsh.option.setNumber("Mesh.LcIntegrationPrecision",1e-3)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin",lc)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax",lc)
        gmsh.model.mesh.generate(1)
        nodes = gmsh.model.mesh.getNodes(1,ltag[0],includeBoundary=True)[1].reshape([-1,3])
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin",0)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax",1e22)
        gmsh.model.remove()
        return nodes[:,:2]


class DistanceField :

    def __init__(self, domain, sampling_size, tags=None) :
        points = []
        for curve in domain._curves :
            if (tags is None) or (curve.tag in tags) :
                points.append(curve.sample(sampling_size))
        points = np.vstack(points)
        self._tree = cKDTree(points)

    def __call__(self,x) :
        x = x[:,:2]
        return self._tree.query(x)[0]


class GeoTiffField :

    def __init__(self,filename) :
        src_ds = gdal.Open(filename)
        self.geo_matrix = src_ds.GetGeoTransform()
        self.data = src_ds.GetRasterBand(1).ReadAsArray()
        assert(self.geo_matrix[2] == 0.)
        assert(self.geo_matrix[4] == 0.)

    def __call__(self,x) :
        gm = self.geo_matrix
        xi = (x[:,0]-gm[3])/gm[5]
        eta = (x[:,1]-gm[0])/gm[1]
        return self.data[xi.astype(int),eta.astype(int)]


class Domain :

    def __init__(self, projection) :
        self._projection = projection
        self._interior_curves = []
        self._interior_points = []
        self._curves = []

    def add_points(self,points,projection) :
        points = np.asarray(points)[:,:2]
        if not projection.IsSame(self._projection) :
            logging.info("reprojecting coordinates")
            points = osr.CoordinateTransformation(projection,self._projection).TransformPoints(points)
            points = np.asarray(points)
        if points.shape[1] == 2 :
            points = np.column_stack([points,np.zeros((points.shape[0],1))])
        self._interior_points.append(points)

    def add_curve(self, points, tag, projection, interior=False, curve_type=BSPLINE) :
        points = np.asarray(points)[:,:2]
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
        curvesiter = itertools.chain(self._curves,self._interior_curves)
        allpoints = np.row_stack(list(itertools.chain((c.points for c in curvesiter),self._interior_points)))
        self._points,unique_id,nbreakpoints = _generate_unique_points(allpoints)
        #assign points id in curves
        cid = 0
        for c in itertools.chain(self._curves,self._interior_curves) :
            n = c.points.shape[0]
            c.pointsid = unique_id[cid:cid+n]
            cid += n
        # assign points id for interior points
        self._interior_points_id = unique_id[cid:]
        # break curves on repeated points
        def split_curves(curves,nbk) :
            newcurves = []
            for curve in curves:
                breaks = np.where(curve.pointsid[1:-1] < nbk)[0]
                if len(breaks) == 0:
                    newcurves.append(curve)
                else :
                    breaks = [0] + list(breaks+1) + [len(curve.points)-1]
                    for i,j in zip(breaks[:-1],breaks[1:]) :
                        ncurve = _Curve(curve.points[i:j+1],curve.tag,curve.curve_type)
                        ncurve.pointsid = curve.pointsid[i:j+1]
                        newcurves.append(ncurve)
            return newcurves
        self._curves = split_curves(self._curves,nbreakpoints)
        self._interior_curves = split_curves(self._interior_curves,nbreakpoints)
        p2curve = np.full([self._points.shape[0],2,2],-1,int)
        curve2p = np.ndarray([len(self._curves),2],int)
        for icurve,c in enumerate(self._curves) :
            print("curve",icurve, "from", c.pointsid[0],"to",c.pointsid[-1])
            curve2p[icurve,0] = c.pointsid[0]
            curve2p[icurve,1] = c.pointsid[-1]
            p0 = p2curve[c.pointsid[0]]
            i0 = 0 if p0[0][0] == -1 else 1
            p0[i0][0] = icurve
            p0[i0][1] = 0
            p1 = p2curve[c.pointsid[-1]]
            i1 = 0 if p1[0][0] == -1 else 1
            p1[i1][0] = icurve
            p1[i1][1] = 1
        touched = np.full((len(self._curves)),False,np.bool)
        self._curveloops = []
        count = 0
        try :
            for i,_ in enumerate(self._curves) :
                if touched[i] : continue
                self._curveloops.append([])
                j = 0
                while (not self._curveloops[-1]) or self._curveloops[-1][0][0] != i:
                    self._curveloops[-1].append((i,j==0))
                    p = curve2p[i,(j+1)%2]
                    i,j = p2curve[p][0] if p2curve[p][0][0] != i else p2curve[p][1]
                    touched[i] = True
                    assert(i!=-1)
                    count += 1
        except :
            raise ValueError("Invalid topology")
        loopbboxarea = np.zeros([len(self._curveloops)])
        for i,loop in enumerate(self._curveloops) :
            lpts = np.row_stack([self._curves[j].points for j,o in loop])
            bbox = np.max(lpts,axis=0)-np.min(lpts,axis=0)
            loopbboxarea[i] = bbox[0]*bbox[1]
        self._curveloops.insert(0,self._curveloops.pop(np.argmax(loopbboxarea)))

    def _add_geometry(self,geometry,tag,projection,curve_type,interior,onlypoints=False) :
        if geometry.GetGeometryCount() != 0 :
            for subg in geometry :
                self._add_geometry(subg,tag,projection,curve_type,interior,onlypoints)
            return
        if onlypoints :
            self.add_points(geometry.GetPoints(),projection)
        else :
            assert(geometry.GetGeometryType() == 2)
            self.add_curve(geometry.GetPoints(),tag,projection,interior,curve_type)

    def add_interior_points(self, filename, physical_name_field=None) :
        proj = osr.SpatialReference()
        driver = ogr.GetDriverByName('ESRI Shapefile')
        data = driver.Open(filename,0)
        layer = data.GetLayer()
        layerdef = layer.GetLayerDefn()
        layerproj = layer.GetSpatialRef()
        for i in layer :
            self._add_geometry(i.geometry(),None,layerproj,None,None,True)

    def add_shapefile(self, filename, physical_name_field=None, interior=False, curve_type=POLYLINE) :
        proj = osr.SpatialReference()
        driver = ogr.GetDriverByName('ESRI Shapefile')
        data = driver.Open(filename,0)
        layer = data.GetLayer()
        layerdef = layer.GetLayerDefn()
        physfield = None
        if not physical_name_field is None :
            for i in range(layerdef.GetFieldCount()):
                field_name = layerdef.GetFieldDefn(i).GetName()
                if field_name == physical_name_field :
                    physfield = i
            if physfield is None :
                raise ValueError("field '"+physical_name_field+"' not found in shapefile")

        layerproj = layer.GetSpatialRef()
        for i in layer :
            phys = i.GetField(physfield) if physfield else "boundary"
            self._add_geometry(i.geometry(),phys,layerproj,curve_type,interior)

    def clip_circle(self, center, radius, resolution) :
        pass

    def unrefine_boundaries(self,x0,y0,sampling,mesh_size_cb) :
        #x = np.vstack(list(c.points for c in self._curves))
        x = np.vstack(list(c.sample(sampling) for c in self._curves))
        x,_,_ = _generate_unique_points(x)
        x = np.copy(x[:,:2])
        #avoid cocircle points
        eps = (np.max(x,axis=0,keepdims=True)-np.min(x,axis=0,keepdims=True))*1e-8
        x = x + np.random.random(x.shape)*eps
        #tag = np.vstack(np.full(c.points.shape[0],c.tag,np.int32) for c in self._curves)
        tag = np.hstack(list(np.full(c.points.shape[0],i,np.int32) for i,c in enumerate(self._curves)))
        def np2c(a,dtype=np.float64) :
            tmp = np.require(a,dtype,"C")
            r = c.c_void_p(tmp.ctypes.data)
            r.tmp = tmp
            return r
        tri = Delaunay(x).simplices
        n_ll = c.c_int()
        n_l = c.c_int()
        n_xo = c.c_int()
        p_xo = c.POINTER(c.c_double)()
        p_l = c.POINTER(c.c_int)()
        p_ll = c.POINTER(c.c_int)()
        mesh_size = mesh_size_cb(x)
        lib.gen_boundaries_from_points(
                c.c_int(x.shape[0]),np2c(x),np2c(tag,np.int32),
                c.c_int(tri.shape[0]),np2c(tri,np.int32),
                c.c_double(x0),c.c_double(y0),
                np2c(mesh_size),
                c.byref(p_xo),c.byref(n_xo),
                c.byref(p_l),c.byref(n_l),
                c.byref(p_ll),c.byref(n_ll))
        xo = np.ctypeslib.frombuffer(c.cast(p_xo,c.POINTER(n_xo.value*2*c.c_double)).contents,dtype=np.float64).reshape([-1,2]).copy()
        lib.libcfree(p_xo)
        l = np.ctypeslib.frombuffer(c.cast(p_l,c.POINTER(n_l.value*3*c.c_int)).contents,dtype=np.int32).reshape([-1,3]).copy()
        lib.libcfree(p_l)
        #ll = np.ctypeslib.frombuffer(c.cast(p_ll,c.POINTER(n_ll.value*2*c.c_int)).contents,dtype=np.int32).reshape([-1,2]).copy()
        #lib.libcfree(p_ll)
        self._curves = []
        for i,(lfrom,lto,llast) in enumerate(l):
            pts = list(range(lfrom,lto+1))+[llast]
            self.add_curve(xo[pts],"boundaries",self._projection,curve_type=POLYLINE)


def mesh_gmsh(domain,filename,mesh_size,version=4.0,gmsh_options={}) :
    nadapt = 3
    nadapt1d = 4
    for curve in domain._curves :
        curve.mesh_size = mesh_size(curve.points)
    gmsh.option.setNumber("Mesh.CharacteristicLengthFromPoints",0)
    gmsh.option.setNumber("Mesh.LcIntegrationPrecision",1e-3)
    gmsh.option.setNumber("Mesh.CharacteristicLengthFactor",1)
    gmsh.option.setNumber("Mesh.CharacteristicLengthFromCurvature",0)
    gmsh.option.setNumber("Mesh.CharacteristicLengthExtendFromBoundary",0)
    gmsh.option.setNumber("Mesh.CharacteristicLengthFromParametricPoints",1)
    allctag = []
    cltag = []
    physicals = {}
    for i,x in enumerate(domain._points) :
        gmsh.model.geo.addPoint(*x,0,tag=i+1)
    for cl in domain._curveloops :
        ctag = []
        for i,(l,o) in enumerate(cl) :
            curve = domain._curves[l]
            tags = _gmsh_geo(curve.curve_type, curve.pointsid)
            physicals.setdefault(curve.tag,[]).extend(tags)
            ctag.extend(tags)
        cltag.append(gmsh.model.geo.addCurveLoop(ctag))
        allctag += ctag
    stag = gmsh.model.geo.addPlaneSurface(cltag)
    embeded = []
    for curve in domain._interior_curves :
        tags = _gmsh_geo(curve.curve_type, curve.pointsid)
        embeded.extend(tags)
        physicals.setdefault(curve.tag,[]).extend(tags)
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.embed(1,embeded,2,stag)
    gmsh.model.mesh.embed(0,[i+1 for i in domain._interior_points_id],2,stag)
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
    tag = gmsh.model.addPhysicalGroup(2,[stag])
    gmsh.model.setPhysicalName(2,stag,"domain")
    ## 1D mesh ##
    gmsh.model.mesh.generate(1)
    for i in range(nadapt1d) :
        for (dim,tag) in gmsh.model.getEntities(1) :
            _,x,u = gmsh.model.mesh.getNodes(dim,tag,includeBoundary=True)
            size = mesh_size(x.reshape([-1,3]))
            gmsh.model.mesh.setSizeAtParametricPoints(dim,tag,u,size)
        print("pass ",i,nadapt)
        gmsh.model.mesh.clear()
        gmsh.model.mesh.generate(1)
    #gmsh.fltk.run()
    ## 2D mesh ##
    initial_sampling=0.1
    gmsh.model.mesh.generate(1)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin",initial_sampling)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax",initial_sampling)
    gmsh.option.setNumber("Mesh.CharacteristicLengthExtendFromBoundary",0)
    gmsh.option.setNumber("Mesh.CharacteristicLengthFromParametricPoints",0)
    gmsh.model.mesh.generate(2)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin",0)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax",1e22)
    bg_field = gmsh.model.mesh.field.add("PostView");
    gmsh.model.mesh.field.setAsBackgroundMesh(bg_field);
    for i in range(nadapt) :
        node_tags, node_x, _  = gmsh.model.mesh.getNodes(-1,-1)
        node_x = node_x.reshape([-1,3])
        nodes_map = dict({tag:i for i,tag in enumerate(node_tags)})
        sf_view = gmsh.view.add("mesh size field")
        node_lc = mesh_size(node_x)
        tri = gmsh.model.mesh.getElements(2,stag)[2][0]
        tri = np.array(list(nodes_map[t] for t in tri)).reshape([-1,3])
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


def convert_msh_gis(inputfn, outputfn) :
    gmsh.model.add(str(uuid.uuid4()))
    gmsh.open(inputfn)
    if outputfn.endswith(".shp") :
        shpdriver = ogr.GetDriverByName('ESRI Shapefile')
    elif outputfn.endswith(".gpkg") :
        shpdriver = ogr.GetDriverByName('GPKG')
    else :
        raise ValueError("Unknown file extension '" + outputfn+"'")
    if os.path.exists(outputfn):
        shpdriver.DeleteDataSource(outputfn)
    out_data_source = shpdriver.CreateDataSource(outputfn)
    out_layer = out_data_source.CreateLayer(outputfn, geom_type=ogr.wkbPolygon)
    out_layer.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))
    id_field = out_layer.FindFieldIndex("id",1)

    tri_i,tri_n = gmsh.model.mesh.getElementsByType(2)
    node_tags, x, _ = gmsh.model.mesh.getNodesByElementType(2)
    x = x.reshape([-1,3])
    nodes_map = dict({tag:i for i,tag in enumerate(node_tags)})
    tri_n = np.array(list(nodes_map[t] for t in tri_n)).reshape([-1,3])

    tri_n = np.column_stack([tri_n,tri_n[:,[0]]])
    header = struct.pack("<biii",1,3,1,4)
    out_data_source.StartTransaction()
    for i,tn in zip(tri_i,tri_n) :
        wkb = header+x[tn,:2].astype('<f8').tobytes()
        feat = ogr.Feature(out_layer.GetLayerDefn())
        feat.SetGeometryDirectly(ogr.CreateGeometryFromWkb(wkb))
        feat.SetField(id_field,int(i))
        out_layer.CreateFeature(feat)
        feat = None
    out_data_source.CommitTransaction()
    gmsh.model.remove()


