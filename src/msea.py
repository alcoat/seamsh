from osgeo import ogr,osr,gdal
import numpy as np
from scipy.spatial import cKDTree
from scipy import interpolate
import logging
import gmsh
import os
import struct
import uuid
import itertools

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)

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
    
    def _gmsh_geo(self) :
        pts = list(i+1 for i in self._pointsid)
        if self.curve_type == SPLINE :
            return [gmsh.model.geo.addSpline(pts)]
        elif self.curve_type == BSPLINE :
            return [gmsh.model.geo.addBSpline(pts)]
        elif self.curve_type == STRICTPOLYLINE :
            return list(gmsh.model.geo.addLine(p0,p1)for p0,p1 in zip(pts[:-1],pts[1:]))
        elif self.curve_type == POLYLINE :
            return [gmsh.model.geo.addPolyline(pts)]

    def sample(self,lc) :
        gmsh.model.add(str(uuid.uuid4()))
        gmsh.model.geo.addPoint(*self.points[0,:],tag=self._pointsid[0]+1)
        for i,x in zip(self._pointsid[1:],self.points[1:]) :
            if i != self._pointsid[0] :
                gmsh.model.geo.addPoint(*x,tag=i+1)
        self._gmsh_geo()
        gmsh.model.geo.synchronize()
        gmsh.option.setNumber("Mesh.LcIntegrationPrecision",1e-3)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin",lc)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax",lc)
        gmsh.model.mesh.generate(1)
        nodes = gmsh.model.mesh.getNodes()[1].reshape([-1,3])
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin",0)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax",1e22)
        gmsh.model.remove()
        return nodes


class DistanceField :

    def __init__(self, domain, tags) :
        points = []
        for curve in domain._curves :
            if curve.tag in tags :
                points.append(curve.sample(0.01))
        points = np.vstack(points)
        self._tree = cKDTree(points)

    def __call__(self,x) :
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
        # generate unique points
        curvesiter = itertools.chain(self._curves,self._interior_curves)
        allpoints = np.row_stack([c.points for c in curvesiter])
        tree = cKDTree(allpoints)
        unique_id = np.full([allpoints.shape[0]],-1,np.int32)
        epoints = np.row_stack([np.array([c.points[0],c.points[-1]]) for c in self._curves])
        bbmin = np.min(epoints,axis=0)
        bbmax = np.max(epoints,axis=0)
        eps = np.linalg.norm(bbmax-bbmin)*1e-8
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
        self._points = np.ndarray([cid,3])
        self._points[unique_id,:] = allpoints
        
        # assign points id in curves
        isendpoint = np.full([cid],False,bool)
        cid = 0
        for c in itertools.chain(self._curves,self._interior_curves) :
            n = c.points.shape[0]
            c._pointsid = unique_id[cid:cid+n]
            cid += n

        # break curves on repeated points
        def split_curves(curves,nbk) :
            newcurves = []
            for curve in curves:
                breaks = np.where(curve._pointsid[1:-1] < nbk)[0]
                if len(breaks) == 0:
                    newcurves.append(curve)
                else :
                    breaks = [0] + list(breaks+1) + [len(curve.points)-1]
                    for i,j in zip(breaks[:-1],breaks[1:]) :
                        print("break",i,j)
                        ncurve = _Curve(curve.points[i:j+1],curve.tag,curve.curve_type)
                        ncurve._pointsid = curve._pointsid[i:j+1]
                        newcurves.append(ncurve)
            return newcurves
        self._curves = split_curves(self._curves,nbreakpoints)
        self._interior_curves = split_curves(self._interior_curves,nbreakpoints)

        p2curve = np.full([self._points.shape[0],2,2],-1,int)
        curve2p = np.ndarray([len(self._curves),2],int)
        for icurve,c in enumerate(self._curves) :
            curve2p[icurve,0] = c._pointsid[0]
            curve2p[icurve,1] = c._pointsid[-1]
            p0 = p2curve[c._pointsid[0]]
            i0 = 0 if p0[0][0] == -1 else 1
            p0[i0][0] = icurve
            p0[i0][1] = 0
            p1 = p2curve[c._pointsid[-1]]
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
                    count += 1
        except :
            raise ValueError("Invalid topology")
        loopbboxarea = np.zeros([len(self._curveloops)])
        for i,loop in enumerate(self._curveloops) :
            lpts = np.row_stack([self._curves[j].points for j,o in loop])
            bbox = np.max(lpts,axis=0)-np.min(lpts,axis=0)
            loopbboxarea[i] = bbox[0]*bbox[1]
        self._curveloops.insert(0,self._curveloops.pop(np.argmax(loopbboxarea)))


    def _add_geometry(self,geometry,tag,projection,curve_type,interior) :
        if geometry.GetGeometryCount() != 0 :
            for subg in geometry :
                self._add_geometry(subg,tag,projection,curve_type,interior)
            return
        assert(geometry.GetGeometryType() == 2)
        self.add_curve(geometry.GetPoints(),tag,projection,interior,curve_type)

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

    def unrefine_geometry(self) :
        pass


def mesh_gmsh(domain,filename,mesh_size,version=4.0,gmsh_options={}) :
    nadapt = 3
    for curve in domain._curves :
        curve.mesh_size = mesh_size(curve.points)
    gmsh.option.setNumber("Mesh.CharacteristicLengthFromPoints",0)
    gmsh.option.setNumber("Mesh.LcIntegrationPrecision",1e-1)
    gmsh.option.setNumber("Mesh.CharacteristicLengthFactor",1)
    gmsh.option.setNumber("Mesh.CharacteristicLengthFromCurvature",0)
    gmsh.option.setNumber("Mesh.CharacteristicLengthExtendFromBoundary",0)
    gmsh.option.setNumber("Mesh.CharacteristicLengthFromParametricPoints",1)
    allctag = []
    cltag = []
    physicals = {}
    for i,x in enumerate(domain._points) :
        gmsh.model.geo.addPoint(*x,tag=i+1)
    for cl in domain._curveloops :
        ctag = []
        for i,(l,o) in enumerate(cl) :
            curve = domain._curves[l]
            tags = curve._gmsh_geo()
            physicals.setdefault(curve.tag,[]).extend(tags)
            ctag.extend(tags)
        cltag.append(gmsh.model.geo.addCurveLoop(ctag))
        allctag += ctag
    stag = gmsh.model.geo.addPlaneSurface(cltag)
    embeded = []
    for curve in domain._interior_curves :
        tags = curve._gmsh_geo()
        embeded.extend(tags)
        physicals.setdefault(curve.tag,[]).extend(tags)
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.embed(1,embeded,2,stag)
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
    for i in range(nadapt) :
        for (dim,tag) in gmsh.model.getEntities(1) :
            _,x,u = gmsh.model.mesh.getNodes(dim,tag,includeBoundary=True)
            size = mesh_size(x.reshape([-1,3]))
            gmsh.model.mesh.setSizeAtParametricPoints(dim,tag,u,size)
        print("pass ",i,nadapt)
        gmsh.model.mesh.clear()
        gmsh.model.mesh.generate(1)
    gmsh.fltk.run()
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
